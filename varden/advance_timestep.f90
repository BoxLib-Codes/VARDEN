module advance_module 

  use BoxLib
  use omp_module
  use f2kcli
  use bl_constants_module
  use list_box_module
  use ml_boxarray_module
  use ml_multifab_module
  use layout_module
  use multifab_module
  use bc_module
  use define_bc_module
  use pre_advance_module
  use velocity_advance_module
  use scalar_advance_module
  use macproject_module
  use hgproject_module
  use rhohalf_module
  use explicit_diffusive_module
  use viscous_module
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

    integer    :: i,n,comp,dm,nlevs
    real(dp_t) :: nrm1,nrm2,nrm3
    logical    :: umac_nodal_flag(mla%dim)

    dm    = mla%dim
    nlevs = mla%nlevel

    do n = 1,nlevs
       call multifab_build(   lapu(n), mla%la(n),    dm, 0)
       call multifab_build(rhohalf(n), mla%la(n),    dm, 1)

       do i = 1,dm
         umac_nodal_flag(:) = .false.
         umac_nodal_flag(i) = .true.
         call multifab_build( umac(n,i), mla%la(n),1,1,nodal=umac_nodal_flag)
         call setval( umac(n,i),1.d20, all=.true.)
       end do

    end do

    if ( verbose .ge. 1 ) then
           do n = 1,nlevs
              nrm1 = norm_inf(uold(n),1,1)
              nrm2 = norm_inf(uold(n),2,1)
              if (dm > 2) nrm3 = norm_inf(uold(n),3,1)
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
    end if

    ! compute lapu
    if (visc_coef .gt. ZERO .and. diffusion_type .eq. 1) then
       do comp = 1, dm
         call get_explicit_diffusive_term(mla,lapu,uold,comp,comp,dx,the_bc_tower)
       end do
    else
       do n = 1, nlevs
         call setval(lapu(n),ZERO)
       enddo
    endif

    call advance_premac(mla,uold,sold,lapu,umac,gp,ext_vel_force,dx,dt, &
                        the_bc_tower%bc_tower_array)

    call macproject(mla,umac,sold,dx,the_bc_tower,press_comp)

    call scalar_advance(mla,uold,sold,snew,umac,ext_scal_force, &
                        dx,dt,the_bc_tower)
    
    call make_at_halftime(mla,rhohalf,sold,snew,1,1,the_bc_tower%bc_tower_array)

    call velocity_advance(mla,uold,unew,sold,lapu,rhohalf,umac,gp, &
                              ext_vel_force,dx,dt,the_bc_tower)

    ! Project the new velocity field.
    call hgproject(proj_type,mla,unew,uold,rhohalf,p,gp,dx,dt, &
                   the_bc_tower,press_comp)

    if ( verbose .ge. 1 ) then
           do n = 1,nlevs
              nrm1 = norm_inf(unew(n),1,1)
              nrm2 = norm_inf(unew(n),2,1)
              if (dm > 2) nrm3 = norm_inf(unew(n),3,1)
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
    end if

    do n = 1,nlevs
       call multifab_destroy(lapu(n))
       call multifab_destroy(rhohalf(n))
       do i = 1,dm
         call multifab_destroy(umac(n,i))
       end do
    end do

1000 format('LEVEL: ',i3,' ITER: ',   i3,' UOLD/VOLD MAX: ',e15.9,1x,e15.9)
1001 format('LEVEL: ',i3,' TIME: ',f14.8,' UOLD/VOLD MAX: ',e15.9,1x,e15.9)
1002 format('LEVEL: ',i3,' ITER: ',   i3,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)
1003 format('LEVEL: ',i3,' TIME: ',f14.8,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)

2000 format('LEVEL: ',i3,' ITER: ',   i3,' UOLD/VOLD/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)
2001 format('LEVEL: ',i3,' TIME: ',f14.8,' UOLD/VOLD/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)
2002 format('LEVEL: ',i3,' ITER: ',   i3,' UNEW/VNEW/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)
2003 format('LEVEL: ',i3,' TIME: ',f14.8,' UNEW/VNEW/WOLD MAX: ',e15.9,1x,e15.9,1x,e15.9)

  end subroutine advance_timestep

end module advance_module 
