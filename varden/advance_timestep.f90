module advance_module 

  use BoxLib
  use f2kcli
  use bl_constants_module
  use list_box_module
  use ml_boxarray_module
  use ml_multifab_module
  use layout_module
  use multifab_module
  use bc_module
  use define_bc_module
  use pre_advance_module       , only : advance_premac
  use velocity_advance_module  , only : velocity_advance
  use scalar_advance_module    , only : scalar_advance
  use macproject_module        , only : macproject
  use hgproject_module         , only : hgproject
  use rhohalf_module           , only : make_at_halftime
  use explicit_diffusive_module, only : get_explicit_diffusive_term
  use probin_module, only : nscal, visc_coef, diffusion_type, stencil_order, &
                            verbose, mg_verbose, cg_verbose
  use proj_parameters

contains

  subroutine advance_timestep(istep,mla,sold,uold,snew,unew,gp,p,ext_vel_force,ext_scal_force,&
                              the_bc_tower,dt,time,dx,press_comp,proj_type)

    implicit none

    integer        , intent(in   ) :: istep
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: unew(:)
    type(multifab) , intent(inout) ::   gp(:)
    type(multifab) , intent(inout) ::    p(:)
    type(multifab) , intent(inout) :: ext_vel_force(:)
    type(multifab) , intent(inout) :: ext_scal_force(:)
    real(dp_t)     , intent(in   ) :: dt,time,dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp
    integer        , intent(in   ) :: proj_type

    type(multifab) :: lapu(mla%nlevel)
    type(multifab) :: umac(mla%nlevel, mla%dim)
    type(multifab) :: rhohalf(mla%nlevel)
    type(multifab) :: mac_rhs(mla%nlevel)

    real(kind=dp_t) ::  sa_time,  sa_time_start,  sa_time_max
    real(kind=dp_t) ::  va_time,  va_time_start,  va_time_max
    real(kind=dp_t) :: mac_time, mac_time_start, mac_time_max
    real(kind=dp_t) ::  hg_time,  hg_time_start,  hg_time_max

    integer    :: i,n,comp,dm,nlevs

    dm    = mla%dim
    nlevs = mla%nlevel

    do n = 1,nlevs
       ! We must build mac_rhs with one ghost cell because if viscous, 
       !    the RHS includes a term that looks like grad(divu)
       call multifab_build(mac_rhs(n), mla%la(n),     1, 1)
       call multifab_build(   lapu(n), mla%la(n),    dm, 0)
       call multifab_build(rhohalf(n), mla%la(n),    dm, 1)

       call setval(mac_rhs(n),0.d0,all=.true.)

       do i = 1,dm
         call multifab_build_edge( umac(n,i), mla%la(n),1,1,i)
         call setval( umac(n,i),1.d20, all=.true.)
       end do

    end do

    if ( verbose .ge. 1 ) call print_old(uold,proj_type,time,istep)

    ! compute lapu
    if (visc_coef .gt. ZERO) then
       do comp = 1, dm
         call get_explicit_diffusive_term(mla,lapu,uold,comp,comp,dx,the_bc_tower)
       end do
    else
       do n = 1, nlevs
         call setval(lapu(n),ZERO)
       enddo
    endif

    call advance_premac(mla,uold,sold,lapu,umac,gp,ext_vel_force,dx,dt,the_bc_tower)
  
    mac_time_start = parallel_wtime()
 
    call macproject(mla,umac,sold,mac_rhs,dx,the_bc_tower,press_comp)

    call parallel_barrier()
    mac_time = parallel_wtime() - mac_time_start

    sa_time_start = parallel_wtime()
    call scalar_advance(mla,sold,snew,umac,ext_scal_force, &
                        dx,dt,the_bc_tower)
    call parallel_barrier()
    sa_time = parallel_wtime() - sa_time_start
    
    call make_at_halftime(mla,rhohalf,sold,snew,1,1,the_bc_tower%bc_tower_array)

    if (diffusion_type .eq. 2) then
       do n = 1, nlevs
          call setval(lapu(n),ZERO)
       enddo
    end if

    va_time_start = parallel_wtime()
    call velocity_advance(mla,uold,unew,sold,lapu,rhohalf,umac,gp, &
                          ext_vel_force,mac_rhs,dx,dt,the_bc_tower)
    call parallel_barrier()
    va_time = parallel_wtime() - va_time_start

    ! Project the new velocity field.
    hg_time_start = parallel_wtime()
    call hgproject(proj_type,mla,unew,uold,rhohalf,p,gp,dx,dt, &
                   the_bc_tower,press_comp)
    call parallel_barrier()
    hg_time = parallel_wtime() - hg_time_start

    if ( verbose .ge. 1 ) call print_new(unew,proj_type,time,dt,istep)

    do n = 1,nlevs
       call multifab_destroy(lapu(n))
       call multifab_destroy(rhohalf(n))
       call multifab_destroy(mac_rhs(n))
       do i = 1,dm
         call multifab_destroy(umac(n,i))
       end do
    end do

    call parallel_reduce(mac_time_max, mac_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())
    call parallel_reduce( hg_time_max,  hg_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())
    call parallel_reduce( sa_time_max,  sa_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())
    call parallel_reduce( va_time_max,  va_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    if (parallel_IOProcessor()) then
       print *, 'Timing summary:'
       print *, 'Scalar   update: ',  sa_time_max, ' seconds'
       print *, 'Velocity update: ',  va_time_max, ' seconds'
       print *, ' MAC Projection: ', mac_time_max, ' seconds'
       print *, '  HG Projection: ',  hg_time_max, ' seconds'
       print *, ' '
    endif

  end subroutine advance_timestep

  subroutine print_old(u,proj_type,time,istep)

     use probin_module, only : nlevs

     type(multifab), intent(in) :: u(:)
     integer       , intent(in) :: proj_type
     real(dp_t)    , intent(in) :: time
     integer       , intent(in) :: istep

     real(dp_t) :: nrm1, nrm2, nrm3
     integer    :: n,dm

     dm = u(1)%dim

     do n = 1,nlevs
         nrm1 = norm_inf(u(n),1,1)
         nrm2 = norm_inf(u(n),2,1)
         if (dm > 2) nrm3 = norm_inf(u(n),3,1)
         if ( parallel_IOProcessor() .and. dm .eq. 2) then
            if ( proj_type .eq. pressure_iters) then
              write(6,1000) n,istep,nrm1,nrm2
            else if ( proj_type .eq. regular_timestep) then
              write(6,1001) n,time,nrm1,nrm2
            else 
              call bl_error('UNKNOWN PROJ_TYPE IN ADVANCE ')
            end if
         else if ( parallel_IOProcessor() .and. dm .eq. 3) then
            if ( proj_type .eq. pressure_iters) then
              write(6,2000) n,istep,nrm1,nrm2,nrm3
            else if ( proj_type .eq. regular_timestep) then
              write(6,2001) n,time,nrm1,nrm2,nrm3
            else 
              call bl_error('UNKNOWN PROJ_TYPE IN ADVANCE ')
            end if
         end if
     end do
     if ( parallel_IOProcessor() ) print *,' '

1000 format('LEVEL: ',i3,' ITER: ',   i3,' UOLD/VOLD MAX: ',e15.9,1x,e15.9)
1001 format('LEVEL: ',i3,' TIME: ',f14.8,' UOLD/VOLD MAX: ',e15.9,1x,e15.9)
2000 format('LEVEL: ',i3,' ITER: ',   i3,' UOLD/VOLD/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)
2001 format('LEVEL: ',i3,' TIME: ',f14.8,' UOLD/VOLD/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)

  end subroutine print_old

  subroutine print_new(u,proj_type,time,dt,istep)

     use probin_module, only : nlevs

     type(multifab), intent(in) :: u(:)
     integer       , intent(in) :: proj_type
     real(dp_t)    , intent(in) :: time,dt
     integer       , intent(in) :: istep

     real(dp_t) :: nrm1, nrm2, nrm3
     integer    :: n,dm

     dm = u(1)%dim

      do n = 1,nlevs
         nrm1 = norm_inf(u(n),1,1)
         nrm2 = norm_inf(u(n),2,1)
         if (dm > 2) nrm3 = norm_inf(u(n),3,1)
         if ( parallel_IOProcessor() .and. dm .eq. 2) then
            if ( proj_type .eq. pressure_iters) then
              write(6,1002) n,istep,nrm1,nrm2
            else if ( proj_type .eq. regular_timestep) then
              write(6,1003) n,time+dt,nrm1,nrm2
            else 
              call bl_error('UNKNOWN PROJ_TYPE IN ADVANCE ')
            end if
         else if ( parallel_IOProcessor() .and. dm .eq. 3) then
            if ( proj_type .eq. pressure_iters) then
              write(6,2002) n,istep,nrm1,nrm2,nrm3
            else if ( proj_type .eq. regular_timestep) then
              write(6,2003) n,time+dt,nrm1,nrm2,nrm3
       else 
              call bl_error('UNKNOWN PROJ_TYPE IN ADVANCE ')
            end if
         end if
      end do
      if ( parallel_IOProcessor() ) print *,' '

1002 format('LEVEL: ',i3,' ITER: ',   i3,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)
1003 format('LEVEL: ',i3,' TIME: ',f14.8,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)
2002 format('LEVEL: ',i3,' ITER: ',   i3,' UNEW/VNEW/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)
2003 format('LEVEL: ',i3,' TIME: ',f14.8,' UNEW/VNEW/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)

  end subroutine print_new

end module advance_module 
