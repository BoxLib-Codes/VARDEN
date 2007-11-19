module scalar_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module
  use setbc_module
  use layout_module

  implicit none

contains

   subroutine scalar_advance (nlevs,mla,uold,sold,snew,rhohalf, &
                              umac,sedge,flux,ext_scal_force,&
                              dx,time,dt, &
                              the_bc_level, &
                              diff_coef,&
                              verbose,use_godunov_debug, &
                              use_minion)
 
      integer        , intent(in   ) :: nlevs
      type(ml_layout), intent(inout) :: mla
      type(multifab) , intent(inout) :: uold(:)
      type(multifab) , intent(inout) :: sold(:)
      type(multifab) , intent(inout) :: snew(:)
      type(multifab) , intent(inout) :: rhohalf(:)
      type(multifab) , intent(inout) :: umac(:,:)
      type(multifab) , intent(inout) :: sedge(:,:)
      type(multifab) , intent(inout) :: flux(:,:)
      type(multifab) , intent(inout) :: ext_scal_force(:)

      real(kind=dp_t), intent(inout) :: dx(:,:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level(:)
      real(kind=dp_t), intent(in   ) :: diff_coef
      integer        , intent(in   ) :: verbose
      logical        , intent(in)    :: use_godunov_debug
      logical        , intent(in)    :: use_minion
 
      ! local variables
      type(multifab), allocatable :: scal_force(:), divu(:)
      logical, allocatable        :: is_conservative(:)
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)

      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  dp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)
      real(kind=dp_t), pointer:: fluxpx(:,:,:,:)
      real(kind=dp_t), pointer:: fluxpy(:,:,:,:)
      real(kind=dp_t), pointer:: fluxpz(:,:,:,:)

      integer :: nscal,dm,d,n
      integer :: lo(uold(1)%dim),hi(uold(1)%dim)
      integer :: i,ng_cell,ng_rho
      logical :: is_vel,make_divu
      real(kind=dp_t) :: diff_fac
      real(kind=dp_t) :: half_dt
      type(box) :: fine_domain

      nscal = ncomp(sold(1))

      allocate(scal_force(nlevs),divu(nlevs))
      allocate(is_conservative(nscal))
      is_conservative(1) = .true.
      is_conservative(2) = .false.
      half_dt = HALF * dt

      ng_cell = uold(1)%ng
      ng_rho  = rhohalf(1)%ng
      dm      = uold(1)%dim

      is_vel  = .false.


      do n = 1, nlevs
         call multifab_build(scal_force(n),ext_scal_force(n)%la,nscal,1)
         call multifab_build(divu(n),scal_force(n)%la,1,1)

         call setval(scal_force(n),0.0_dp_t,all=.true.)
         call setval(divu(n),0.0_dp_t,all=.true.)
      enddo

      do n = 1, nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff_fac = ONE
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         fp  => dataptr(scal_force(n), i)
         ep  => dataptr(ext_scal_force(n), i)
         sop => dataptr(sold(n), i)
         lo = lwb(get_box(sold(n), i))
         hi = upb(get_box(sold(n), i))
         select case (dm)
         case (2)
            call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sop(:,:,1,:), &
                                ng_cell, dx(n,:), &
                                the_bc_level(n)%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                diff_coef, diff_fac)
         case (3)
            call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sop(:,:,:,:), &
                                ng_cell, dx(n,:), &
                                the_bc_level(n)%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                diff_coef, diff_fac)
         end select
      end do

      call multifab_fill_boundary(scal_force(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of scalar using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1, sold(n)%nboxes
         if ( multifab_remote(sold(n), i) ) cycle
         sop    => dataptr(sold(n), i)
         uop    => dataptr(uold(n), i)
         sepx   => dataptr(sedge(n,1), i)
         sepy   => dataptr(sedge(n,2), i)
         fluxpx => dataptr(flux(n,1), i)
         fluxpy => dataptr(flux(n,2), i)
         ump    => dataptr(umac(n,1), i)
         vmp    => dataptr(umac(n,2), i)
         fp     => dataptr(scal_force(n) , i)
         dp     => dataptr(divu(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            if(use_godunov_debug) then
               call mkflux_debug_2d(sop(:,:,1,:), uop(:,:,1,:), &
                                    sepx(:,:,1,:), sepy(:,:,1,:), &
                                    fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                    ump(:,:,1,1), vmp(:,:,1,1), &
                                    fp(:,:,1,:), dp(:,:,1,1), &
                                    lo, dx(n,:), dt, is_vel, &
                                    the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                    the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                              sepx(:,:,1,:), sepy(:,:,1,:), &
                              fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                              ump(:,:,1,1), vmp(:,:,1,1), &
                              fp(:,:,1,:), dp(:,:,1,1), &
                              lo, dx(n,:), dt, is_vel, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:), &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                              ng_cell, use_minion, is_conservative)
            endif
         case (3)
            sepz   => dataptr(sedge(n,3), i)
            fluxpz => dataptr(flux(n,3), i)
            wmp  => dataptr(umac(n,3), i)
            if(use_godunov_debug) then
               call mkflux_debug_3d(sop(:,:,:,:), uop(:,:,:,:), &
                                    sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                    fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    fp(:,:,:,:), dp(:,:,:,1), &
                                    lo, dx(n,:), dt, is_vel, &
                                    the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                    the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                              sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                              fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                              ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                              fp(:,:,:,:), dp(:,:,:,1), &
                              lo, dx(n,:), dt, is_vel, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:), &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                              ng_cell, use_minion, is_conservative)
            endif
         end select
      end do

      enddo ! do n = 1, nlevs

      ! sychronize fluxes
      do n = 2, nlevs
         do d = 1, dm
            call ml_edge_restriction(flux(n-1,d),flux(n,d),mla%mba%rr(n-1,:),d)
         enddo
      enddo

      do n = 1, nlevs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n+1/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff_fac = HALF
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         fp  => dataptr(scal_force(n), i)
         ep  => dataptr(ext_scal_force(n), i)
         sop => dataptr(sold(n) , i)
         lo =  lwb(get_box(sold(n), i))
         hi =  upb(get_box(sold(n), i))
         select case (dm)
         case (2)
            call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sop(:,:,1,:), &
                                ng_cell, dx(n,:), &
                                the_bc_level(n)%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                diff_coef, diff_fac)
         case (3)
            call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sop(:,:,:,:), &
                                ng_cell, dx(n,:), &
                                the_bc_level(n)%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                diff_coef, diff_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force(n))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the scalars with conservative or convective differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, sold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         sop  => dataptr(sold(n), i)
         ump  => dataptr(umac(n,1), i)
         vmp  => dataptr(umac(n,2), i)
         sepx => dataptr(sedge(n,1), i)
         sepy => dataptr(sedge(n,2), i)
         fluxpx => dataptr(flux(n,1), i)
         fluxpy => dataptr(flux(n,2), i)
         snp  => dataptr(snew(n), i)
          rp  => dataptr(rhohalf(n), i)
          fp  => dataptr(scal_force(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            call update_2d(n,sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                           sepx(:,:,1,:), sepy(:,:,1,:), &
                           fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                           fp(:,:,1,:) , snp(:,:,1,:), &
                           rp(:,:,1,1) , &
                           lo, hi, ng_cell,dx(n,:),dt,is_vel,is_conservative,verbose)
         case (3)
            wmp => dataptr(umac(n,3), i)
            sepz => dataptr(sedge(n,3), i)
            fluxpz => dataptr(flux(n,3), i)
            call update_3d(n,sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                           sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                           fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                           fp(:,:,:,:) , snp(:,:,:,:), &
                           rp(:,:,:,1) , &
                           lo, hi, ng_cell,dx(n,:),dt,is_vel,is_conservative,verbose)
         end select
      end do
      
      call multifab_fill_boundary(rhohalf(n))
      call multifab_fill_boundary(snew(n))
      
      do i = 1, sold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         snp => dataptr(snew(n), i)
         rp  => dataptr(rhohalf(n), i)
         lo = lwb(get_box(uold(n), i))
         select case (dm)
         case (2)
            do d = 1,nscal
               call setbc_2d(snp(:,:,1,d), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,dm+d), &
                             dx(n,:),dm+d)
            end do
            call setbc_2d(rp(:,:,1,1), lo, ng_rho, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1), &
                          dx(n,:),dm+1)
         case (3)
            do d = 1,nscal
               call setbc_3d(snp(:,:,:,d), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,dm+d), &
                             dx(n,:),dm+d)
            end do
            call setbc_3d(rp(:,:,:,1), lo, ng_rho, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1), &
                          dx(n,:),dm+1)
         end select
      end do
      
      enddo ! do n = 1, nlevs

      ! use restriction so coarse cells are the average
      ! of the corresponding fine cells
      do n = 2, nlevs
         call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))
      enddo

      do n = 2, nlevs
         fine_domain = layout_get_pd(mla%la(n))
         call multifab_fill_ghost_cells(snew(n),snew(n-1),fine_domain, &
              ng_cell,mla%mba%rr(n-1,:), &
              the_bc_level(n-1)%adv_bc_level_array(0,:,:,:), &
              1,dm+1,nscal)
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(is_conservative)

      do n = 1,nlevs
         call multifab_destroy(scal_force(n))
         call multifab_destroy(divu(n))
      enddo

   end subroutine scalar_advance

  subroutine modify_force(a, targ, b, src_b, c, src_c, nc, mult, all)
    integer, intent(in)           :: targ, src_b, src_c
    integer, intent(in), optional :: nc
    logical, intent(in), optional :: all
    type(multifab), intent(inout) :: a
    type(multifab), intent(in   ) :: b
    type(multifab), intent(in   ) :: c
    real(dp_t)    , intent(in   ) :: mult
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)
    integer :: i
    logical :: lall
    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src_b, nc)
          cp => dataptr(c, i, src_c, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src_b, nc)
          cp => dataptr(c, i, get_ibox(c, i), src_c, nc)
       end if
       ap = ap + mult*bp*cp
    end do
    !$OMP END PARALLEL DO
  end subroutine modify_force

end module scalar_advance_module
