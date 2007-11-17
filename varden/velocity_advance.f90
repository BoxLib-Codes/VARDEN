module velocity_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module

  implicit none

contains

   subroutine velocity_advance(nlevs,uold,unew,sold,rhohalf,&
                               umac,uedge,flux,gp,p, &
                               ext_vel_force,dx,time,dt, &
                               the_bc_level, &
                               visc_coef,verbose,use_godunov_debug, &
                               use_minion)
 
      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: uold(:)
      type(multifab) , intent(inout) :: sold(:)
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

      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  dp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)

      integer :: lo(uold(1)%dim),hi(uold(1)%dim)
      integer :: ir(uold(1)%dim)
      integer :: i,n,comp,dm,ng_cell,ng_rho,ilev
      logical :: is_vel,is_conservative(uold(1)%dim)
      real(kind=dp_t) :: visc_fac,visc_mu,half_dt

      allocate(vel_force(nlevs),divu(nlevs))

      ! refinement ratio
      ir = 2

      half_dt = HALF * dt

      ng_cell = uold(1)%ng
      ng_rho  = rhohalf(1)%ng
      dm      = uold(1)%dim

      is_conservative = .false.
      
      do ilev=1,nlevs
         call multifab_build(vel_force(ilev),ext_vel_force(ilev)%la,dm,1)
         call multifab_build(divu(ilev),vel_force(ilev)%la,1,1)
         call setval(divu(ilev),0.0_dp_t,all=.true.)
      enddo

      do ilev=1,nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the velocity forcing term at time n using rho and the full viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      visc_fac = ONE
      do i = 1, uold(ilev)%nboxes
         if ( multifab_remote(uold(ilev), i) ) cycle
         fp  => dataptr(vel_force(ilev), i)
         ep  => dataptr(ext_vel_force(ilev), i)
         gpp => dataptr(gp(ilev), i)
         rp  => dataptr(sold(ilev), i)
         uop => dataptr(uold(ilev), i)
         lo = lwb(get_box(uold(ilev), i))
         hi = upb(get_box(uold(ilev), i))
         select case (dm)
         case (2)
            call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                               gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                               ng_cell, ng_cell, dx(ilev,:), &
                               the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         case (3)
            call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                               gpp(:,:,:,:), rp(:,:,:,1), uop(:,:,:,:), &
                               ng_cell, ng_cell, dx(ilev,:), &
                               the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(vel_force(ilev))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of velocity using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      is_vel = .true.
      do i = 1, uold(ilev)%nboxes
         if ( multifab_remote(uold(ilev), i) ) cycle
         uop  => dataptr(uold(ilev), i)
         uepx => dataptr(uedge(ilev,1), i)
         uepy => dataptr(uedge(ilev,2), i)
         fluxpx => dataptr(flux(ilev,1), i)
         fluxpy => dataptr(flux(ilev,2), i)
         ump  => dataptr(umac(ilev,1), i)
         vmp  => dataptr(umac(ilev,2), i)
         fp   => dataptr(vel_force(ilev), i)
         dp   => dataptr(divu(ilev), i)
         lo = lwb(get_box(uold(ilev), i))
         hi = upb(get_box(uold(ilev), i))
         select case (dm)
         case (2)
            if(use_godunov_debug) then
               call mkflux_debug_2d(uop(:,:,1,:), uop(:,:,1,:), &
                                    uepx(:,:,1,:), uepy(:,:,1,:), &
                                    fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                    ump(:,:,1,1), vmp(:,:,1,1), &
                                    fp(:,:,1,:), dp(:,:,1,1), &
                                    lo, dx(ilev,:), dt, is_vel, &
                                    the_bc_level(ilev)%phys_bc_level_array(i,:,:), &
                                    the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                              uepx(:,:,1,:), uepy(:,:,1,:), &
                              fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                              ump(:,:,1,1), vmp(:,:,1,1), &
                              fp(:,:,1,:), dp(:,:,1,1), &
                              lo, dx(ilev,:), dt, is_vel, &
                              the_bc_level(ilev)%phys_bc_level_array(i,:,:), &
                              the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                              ng_cell, use_minion, is_conservative)
            endif
         case (3)
            uepz => dataptr(uedge(ilev,3), i)
            fluxpz => dataptr(flux(ilev,3), i)
            wmp  => dataptr(umac(ilev,3), i)
            if(use_godunov_debug) then
               call mkflux_debug_3d(uop(:,:,:,:), uop(:,:,:,:), &
                                    uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                                    fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    fp(:,:,:,:), dp(:,:,:,1), &
                                    lo, dx(ilev,:), dt, is_vel, &
                                    the_bc_level(ilev)%phys_bc_level_array(i,:,:), &
                                    the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                              uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                              fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                              ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                              fp(:,:,:,:), dp(:,:,:,1), &
                              lo, dx(ilev,:), dt, is_vel, &
                              the_bc_level(ilev)%phys_bc_level_array(i,:,:), &
                              the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                              ng_cell, use_minion, is_conservative)
            endif
         end select
      end do

      enddo ! do ilev=1,nlevs

      ! synchronize fluxes
      ! for now we don't call ml_edge restriction because currently the velocity is
      ! updated with convective differencing
!      do ilev=2,nlevs
!         do n=1,dm
!            call ml_edge_restriction(flux(ilev-1,n),flux(ilev,n),ir,n)
!         enddo
!      enddo

      do ilev=1,nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Now create vel_force at half-time using rhohalf and half the viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      visc_fac = HALF
      do i = 1, uold(ilev)%nboxes
         if ( multifab_remote(uold(ilev), i) ) cycle
         fp  => dataptr(vel_force(ilev), i)
         ep  => dataptr(ext_vel_force(ilev), i)
         gpp => dataptr(gp(ilev), i)
         rp  => dataptr(rhohalf(ilev), i)
         uop => dataptr(uold(ilev), i)
         lo = lwb(get_box(uold(ilev), i))
         hi = upb(get_box(uold(ilev), i))
         select case (dm)
         case (2)
            call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                               gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                               ng_cell, ng_rho, dx(ilev,:), &
                               the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         case (3)
            call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                               gpp(:,:,:,:), rp(:,:,:,1), uop(:,:,:,:), &
                               ng_cell, ng_rho, dx(ilev,:), &
                               the_bc_level(ilev)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(vel_force(ilev))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the velocity with convective differencing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, uold(ilev)%nboxes
         if ( multifab_remote(uold(ilev), i) ) cycle
         uop  => dataptr(uold(ilev), i)
         ump  => dataptr(umac(ilev,1), i)
         vmp  => dataptr(umac(ilev,2), i)
         uepx => dataptr(uedge(ilev,1), i)
         uepy => dataptr(uedge(ilev,2), i)
         fluxpx => dataptr(flux(ilev,1), i)
         fluxpy => dataptr(flux(ilev,2), i)
         unp  => dataptr(unew(ilev), i)
         fp   => dataptr(vel_force(ilev), i)
         rp   => dataptr(rhohalf(ilev), i)
         lo = lwb(get_box(uold(ilev), i))
         hi = upb(get_box(uold(ilev), i))
         select case (dm)
         case (2)
            call update_2d(ilev,uop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                           uepx(:,:,1,:), uepy(:,:,1,:), &
                           fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                           fp(:,:,1,:), unp(:,:,1,:), &
                           rp(:,:,1,1), &
                           lo, hi, ng_cell,dx(ilev,:),dt,is_vel,is_conservative,verbose)
         case (3)
            wmp => dataptr(umac(ilev,3), i)
            uepz => dataptr(uedge(ilev,3), i)
            fluxpz => dataptr(flux(ilev,3), i)
            call update_3d(ilev,uop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                           uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                           fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                           fp(:,:,:,:), unp(:,:,:,:), &
                           rp(:,:,:,1), &
                           lo, hi, ng_cell,dx(ilev,:),dt,is_vel,is_conservative,verbose)
         end select
      end do

      enddo ! do ilev=1,nlevs

      ! use restriction so coarse cells are the average
      ! of the corresponding fine cells
      do ilev=2,nlevs
         call ml_cc_restriction(unew(ilev-1),unew(ilev),ir)
      enddo

      do ilev=1,nlevs
         call multifab_destroy(vel_force(ilev))
         call multifab_destroy(divu(ilev))
      enddo

   end subroutine velocity_advance

end module velocity_advance_module
