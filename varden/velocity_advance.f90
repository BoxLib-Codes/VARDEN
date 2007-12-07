module velocity_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module

  implicit none

contains

  subroutine velocity_advance(nlevs,mla,uold,unew,sold,lapu,rhohalf,umac,uedge,flux,gp,p, &
                              ext_vel_force,dx,time,dt,the_bc_level,visc_coef,verbose, &
                              use_godunov_debug,use_minion)

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: lapu(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: uedge(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(multifab) , intent(inout) :: unew(:)
    type(multifab) , intent(inout) :: rhohalf(:)
    type(multifab) , intent(inout) :: gp(:)
    type(multifab) , intent(inout) :: p(:)
    type(multifab) , intent(inout) :: ext_vel_force(:)
    real(kind=dp_t), intent(inout) :: dx(:,:),time,dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: visc_coef
    integer        , intent(in   ) :: verbose
    logical        , intent(in)    :: use_godunov_debug
    logical        , intent(in)    :: use_minion

    ! local
    type(multifab), allocatable :: vel_force(:)
    type(multifab), allocatable :: divu(:)

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: unp(:,:,:,:)
    real(kind=dp_t), pointer:: ump(:,:,:,:)
    real(kind=dp_t), pointer:: vmp(:,:,:,:)
    real(kind=dp_t), pointer:: wmp(:,:,:,:)
    real(kind=dp_t), pointer:: gpp(:,:,:,:)
    real(kind=dp_t), pointer::  fp(:,:,:,:)
    real(kind=dp_t), pointer::  ep(:,:,:,:)
    real(kind=dp_t), pointer:: uepx(:,:,:,:)
    real(kind=dp_t), pointer:: uepy(:,:,:,:)
    real(kind=dp_t), pointer:: uepz(:,:,:,:)
    real(kind=dp_t), pointer:: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer:: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer:: fluxpz(:,:,:,:)

    real(kind=dp_t), pointer:: lapup(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: snp(:,:,:,:)
    real(kind=dp_t), pointer::  rp(:,:,:,:) 
    real(kind=dp_t), pointer::  dp(:,:,:,:) 
    real(kind=dp_t), pointer:: sepx(:,:,:,:)
    real(kind=dp_t), pointer:: sepy(:,:,:,:)

    integer :: lo(uold(1)%dim),hi(uold(1)%dim)
    integer :: i,n,dm,d,comp,ng_cell,ng_rho
    logical :: is_vel,is_conservative(uold(1)%dim)
    real(kind=dp_t) :: visc_fac,visc_mu,half_dt
    real(kind=dp_t) :: umin,umax

    allocate(vel_force(nlevs),divu(nlevs))

    half_dt = HALF * dt

    ng_cell = uold(1)%ng
    ng_rho  = rhohalf(1)%ng
    dm      = uold(1)%dim

    is_conservative = .false.

    do n = 1, nlevs
       call multifab_build(vel_force(n),ext_vel_force(n)%la,dm,1)
       call multifab_build(divu(n),vel_force(n)%la,1,1)
       call setval(divu(n),0.0_dp_t,all=.true.)
    enddo


    !********************************************************
    ! Create the velocity forcing term at time n using rho and the full viscous term.
    !********************************************************
    visc_fac = ONE
    call mkvelforce(nlevs,vel_force,ext_vel_force,sold,gp,uold,lapu,dx,visc_coef,visc_fac)

    do n = 1, nlevs


       !********************************************************
       ! Create the edge states of velocity using the MAC velocity 
       !********************************************************
       is_vel = .true.
       do i = 1, uold(n)%nboxes
          if ( multifab_remote(uold(n), i) ) cycle
          uop  => dataptr(uold(n), i)
          uepx => dataptr(uedge(n,1), i)
          uepy => dataptr(uedge(n,2), i)
          fluxpx => dataptr(flux(n,1), i)
          fluxpy => dataptr(flux(n,2), i)
          ump  => dataptr(umac(n,1), i)
          vmp  => dataptr(umac(n,2), i)
          fp   => dataptr(vel_force(n), i)
          dp   => dataptr(divu(n), i)
          lo = lwb(get_box(uold(n), i))
          hi = upb(get_box(uold(n), i))
          select case (dm)
          case (2)
             if(use_godunov_debug) then
                call mkflux_debug_2d(uop(:,:,1,:), uop(:,:,1,:), &
                                     uepx(:,:,1,:), uepy(:,:,1,:), &
                                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                     ump(:,:,1,1), vmp(:,:,1,1), &
                                     fp(:,:,1,:), dp(:,:,1,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                                     ng_cell, use_minion, is_conservative)
             else
                call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                               uepx(:,:,1,:), uepy(:,:,1,:), &
                               fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               fp(:,:,1,:), dp(:,:,1,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               ng_cell, use_minion, is_conservative)
             endif
          case (3)
             uepz => dataptr(uedge(n,3), i)
             fluxpz => dataptr(flux(n,3), i)
             wmp  => dataptr(umac(n,3), i)
             if(use_godunov_debug) then
                call mkflux_debug_3d(uop(:,:,:,:), uop(:,:,:,:), &
                                     uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                                     fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                     ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                     fp(:,:,:,:), dp(:,:,:,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                                     ng_cell, use_minion, is_conservative)
             else
                call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                               uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                               fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               fp(:,:,:,:), dp(:,:,:,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               ng_cell, use_minion, is_conservative)
             endif
          end select
       end do

    enddo ! do n = 1, nlevs

    ! synchronize fluxes
    do n = nlevs, 2, -1
       do comp = 1, dm
          if(is_conservative(comp)) then
             do d = 1, dm
                call ml_edge_restriction_c(flux(n-1,d),comp,flux(n,d),comp, &
                                           mla%mba%rr(n-1,:),d,1)
             enddo
          endif
       enddo
    enddo

    !********************************************************
    ! Now create vel_force at half-time using rhohalf and half the viscous term.
    !********************************************************
    visc_fac = HALF
    call mkvelforce(nlevs,vel_force,ext_vel_force,rhohalf,gp,uold,lapu,dx, &
                    visc_coef,visc_fac)

    !********************************************************
    ! Update the velocity with convective differencing
    !********************************************************
    call update(nlevs,uold,umac,uedge,flux,vel_force,unew,rhohalf,dx,dt,is_vel, &
                is_conservative,the_bc_level,mla)

    do n = 1, nlevs
       call multifab_destroy(vel_force(n))
       call multifab_destroy(divu(n))
    enddo

    if (verbose .ge. 1) then
       do n = 1, nlevs
          do comp = 1, dm
             umin = multifab_min_c(unew(n),comp) 
             umax = multifab_max_c(unew(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2001) n,umin,umax
             else if (comp .eq. 2) then
                if (parallel_IOProcessor()) write(6,2002) n,umin,umax
             else if (comp .eq. 3) then
                if (parallel_IOProcessor()) write(6,2003) n,umin,umax
             end if
          end do
       end do
    end if

2001 format('... level ', i2,' new min/max : x-vel           ',e17.10,2x,e17.10)
2002 format('... level ', i2,' new min/max : y-vel           ',e17.10,2x,e17.10)
2003 format('... level ', i2,' new min/max : z-vel           ',e17.10,2x,e17.10)

  end subroutine velocity_advance

end module velocity_advance_module
