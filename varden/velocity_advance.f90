module velocity_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module

  implicit none

contains

   subroutine velocity_advance(nlevs,mla,uold,unew,sold,lapu,rhohalf,&
                               umac,uedge,flux,gp,p, &
                               ext_vel_force,dx,time,dt, &
                               the_bc_level, &
                               visc_coef,verbose,use_godunov_debug, &
                               use_minion)
 
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
      type(box) :: fine_domain

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

      do n = 1, nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the velocity forcing term at time n using rho and the full viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      visc_fac = ONE
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         fp    => dataptr(vel_force(n), i)
         ep    => dataptr(ext_vel_force(n), i)
         gpp   => dataptr(gp(n), i)
         rp    => dataptr(sold(n), i)
         uop   => dataptr(uold(n), i)
         lapup => dataptr(lapu(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                               gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                               lapup(:,:,1,:), &
                               ng_cell, ng_cell, dx(n,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         case (3)
            call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                               gpp(:,:,:,:), rp(:,:,:,1), uop(:,:,:,:), &
                               lapup(:,:,:,:), &
                               ng_cell, ng_cell, dx(n,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(vel_force(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of velocity using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      do n = 1, nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Now create vel_force at half-time using rhohalf and half the viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      visc_fac = HALF
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         fp    => dataptr(vel_force(n), i)
         ep    => dataptr(ext_vel_force(n), i)
         gpp   => dataptr(gp(n), i)
         rp    => dataptr(rhohalf(n), i)
         uop   => dataptr(uold(n), i)
         lapup => dataptr(lapu(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                               gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                               lapup(:,:,1,:), &
                               ng_cell, ng_rho, dx(n,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         case (3)
            call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                               gpp(:,:,:,:), rp(:,:,:,1), uop(:,:,:,:), &
                               lapup(:,:,:,:), &
                               ng_cell, ng_rho, dx(n,:), &
                               the_bc_level(n)%ell_bc_level_array(i,:,:,1:dm), &
                               visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(vel_force(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the velocity with convective differencing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         uop  => dataptr(uold(n), i)
         ump  => dataptr(umac(n,1), i)
         vmp  => dataptr(umac(n,2), i)
         uepx => dataptr(uedge(n,1), i)
         uepy => dataptr(uedge(n,2), i)
         fluxpx => dataptr(flux(n,1), i)
         fluxpy => dataptr(flux(n,2), i)
         unp  => dataptr(unew(n), i)
         fp   => dataptr(vel_force(n), i)
         rp   => dataptr(rhohalf(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            call update_2d(n,uop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                           uepx(:,:,1,:), uepy(:,:,1,:), &
                           fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                           fp(:,:,1,:), unp(:,:,1,:), &
                           rp(:,:,1,1), &
                           lo, hi, ng_cell,dx(n,:),dt,is_vel,is_conservative,verbose)
         case (3)
            wmp => dataptr(umac(n,3), i)
            uepz => dataptr(uedge(n,3), i)
            fluxpz => dataptr(flux(n,3), i)
            call update_3d(n,uop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                           uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                           fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                           fp(:,:,:,:), unp(:,:,:,:), &
                           rp(:,:,:,1), &
                           lo, hi, ng_cell,dx(n,:),dt,is_vel,is_conservative,verbose)
         end select
      end do

      call multifab_fill_boundary(unew(n))

      do i = 1, unew(n)%nboxes
         if ( multifab_remote(unew(n), i) ) cycle
         unp => dataptr(unew(n), i)
         lo = lwb(get_box(unew(n), i))
         select case (dm)
         case (2)
            do comp = 1, dm
               call setbc_2d(unp(:,:,1,comp), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,comp), &
                             dx(n,:),comp)
            end do
         case (3)
            do comp = 1, dm
               call setbc_3d(unp(:,:,:,comp), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,comp), &
                             dx(n,:),comp)
            end do
         end select
      end do

      enddo ! do n = 1, nlevs

      ! use restriction so coarse cells are the average
      ! of the corresponding fine cells
      do n = nlevs, 2, -1
         call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:))
      enddo

      do n = 2, nlevs
         fine_domain = layout_get_pd(mla%la(n))
         call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                        ng_cell,mla%mba%rr(n-1,:), &
                                        the_bc_level(n-1), the_bc_level(n), &
                                        1,1,dm)
      end do

      do n = 1, nlevs
         call multifab_destroy(vel_force(n))
         call multifab_destroy(divu(n))
      enddo

   end subroutine velocity_advance

end module velocity_advance_module
