module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use explicit_diffusive_module
  use viscous_module
  use fabio_module
  use macproject_module
  use probin_module, only : nscal, diff_coef, diffusion_type, stencil_order, &
                            verbose, mg_verbose, cg_verbose, mass_fractions

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(mla,uold,sold,snew,umac, &
                            ext_scal_force,dx,dt,t,the_bc_tower)

    use mkflux_module
    use mkforce_module
    use update_module
    use bl_constants_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,t
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    type(multifab) :: scal_force(mla%nlevel)
    type(multifab) :: divu(mla%nlevel)
    type(multifab) :: laps(mla%nlevel)
    type(multifab) :: sflux(mla%nlevel,mla%dim)
    type(multifab) :: sedge(mla%nlevel,mla%dim)
    ! stores del dot (us)
    ! filled but not used 
    type(multifab) :: adv_s(mla%nlevel)

    logical :: is_conservative(nscal)
    logical :: umac_nodal_flag(mla%dim)

    integer         :: i,n
    integer         :: comp,bc_comp,nlevs,dm
    logical         :: is_vel
    real(kind=dp_t) :: diff_fac,visc_mu
    real(kind=dp_t) :: smin,smax

    character(len=20)         :: plot_names(nscal)
    character(len=20)         :: sd_name
    integer, save             :: iter = 0


    plot_names(1) = 'Density'
    plot_names(2) = 'SpeciesA'
    plot_names(3) = 'SpeciesB'
    plot_names(4) = 'SpeciesC'

    nlevs = mla%nlevel
    dm    = mla%dim

    
    if (mass_fractions) then
       ! using mass fractions
       is_conservative(:) = .true.
    else
       is_conservative(1) = .true.
       ! using concentrations
       do comp = 2, nscal
          is_conservative(comp) = .false.
       end do
    end if

    is_vel  = .false.

    do n = 1, nlevs
       call multifab_build( scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build( divu(n),mla%la(n),    1,1)
       call multifab_build( laps(n),mla%la(n),nscal,0)
!how many ghost cells?
       call multifab_build( adv_s(n),mla%la(n),nscal-1,0)
       do i = 1,dm
         umac_nodal_flag(:) = .false.
         umac_nodal_flag(i) = .true.
         call multifab_build(sflux(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call multifab_build(sedge(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call setval(sflux(n,i),ZERO,all=.true.)
         call setval(sedge(n,i),ZERO,all=.true.)
       end do

       call setval(scal_force(n),0.0_dp_t,all=.true.)
       call setval(divu(n),0.0_dp_t,all=.true.)
       call setval(laps(n),0.0_dp_t,all=.true.)
    enddo

    !***********************************
    ! Compute laps for passive scalar only
    !***********************************

call write_plotfile(8000+iter,nscal,sold)
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,laps,sold,comp,bc_comp,dx,the_bc_tower,is_vel=.false.)
       enddo
    else
      do n = 1, nlevs
         call setval(laps(n),ZERO)
      enddo
    endif
! REMOVE ME
call write_plotfile(9000+iter,nscal,laps)
    !***********************************
    ! Create scalar force at time n.
    !***********************************

    diff_fac = ONE
    call mkscalforce(nlevs,scal_force,ext_scal_force,laps,diff_fac)
call write_plotfile(9100+iter,nscal,scal_force)
    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************

    call mkflux(mla,sold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_tower%bc_tower_array,is_vel,is_conservative)

    !***********************************
    ! Create scalar force at time n+1/2.
    !***********************************

    ! The laps term will be added to the rhs in diff_scalar_solve
    ! for Crank-Nicolson
    diff_fac = ZERO
    call mkscalforce(nlevs,scal_force,ext_scal_force,laps,diff_fac)
call write_plotfile(9200+iter,nscal,scal_force)
    !***********************************
    ! Update the scalars with conservative or convective differencing.
    !***********************************

    call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s,dx,dt,t,is_vel, &
                is_conservative,the_bc_tower%bc_tower_array)
call write_plotfile(9300+iter,nscal,snew)

    if (verbose .ge. 1) then
       do n = 1, nlevs
          do comp = 1, nscal
             smin = multifab_min_c(snew(n),comp) 
             smax = multifab_max_c(snew(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2000) n,smin,smax
             else if (comp .eq. 2) then
                if (parallel_IOProcessor()) write(6,2001) n,smin,smax
             end if
          end do
       end do
    end if

    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(divu(n))
       do i = 1,dm
         call multifab_destroy(sflux(n,i))
         call multifab_destroy(sedge(n,i))
       end do
    enddo

    if (diff_coef > ZERO) then
       do comp = 2, nscal
           bc_comp = dm+comp
 
           ! Crank-Nicolson
           if (diffusion_type .eq. 1) then
              visc_mu = HALF*dt*diff_coef
 
           ! backward Euler
           else if (diffusion_type .eq. 2) then
              visc_mu = dt*diff_coef
           else
             call bl_error('BAD DIFFUSION TYPE ')
          end if

          call diff_scalar_solve(mla,snew,laps,dx,t,visc_mu,&
                                 the_bc_tower,comp,bc_comp)

        enddo
     end if

!--- for debugging -- FIX ME 
!       do comp = 2, nscal
!          bc_comp = dm+comp
!          call get_explicit_diffusive_term(mla,laps,snew,comp,bc_comp,dx,the_bc_tower,is_vel=.false.)
!       enddo
!call write_plotfile(100+iter,nscal,laps)
call write_plotfile(9400+iter,nscal,snew)
!--------

    do n = 1,nlevs
       call multifab_destroy(laps(n))
       call multifab_destroy(adv_s(n))
    end do

    iter = iter + 1

2000 format('... level ', i2,' new min/max : density           ',e17.10,2x,e17.10)
2001 format('... level ', i2,' new min/max :  tracer           ',e17.10,2x,e17.10)


!-- For debugging -- Fix me
  contains
  
   subroutine write_plotfile(istep_to_write, n_plot_comps, mf)

    integer,         intent(in   ) :: istep_to_write
    integer,         intent(in   ) :: n_plot_comps
    type(multifab),  intent(in   ) :: mf(:)


    type(multifab), allocatable :: plotdata(:) 
    integer        :: n

 
    allocate(plotdata(nlevs))
 
    do n = 1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy(plotdata(n), mf(n))

    end do
        
    write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), t, dx(1,:))

    do n = 1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)
    
  end subroutine write_plotfile

  end subroutine scalar_advance

end module scalar_advance_module
