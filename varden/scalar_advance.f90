module scalar_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module
  use setbc_module

  implicit none

contains

   subroutine scalar_advance (lev,uold,sold,snew,rhohalf,&
                              umac,sedge, &
                              ext_scal_force,&
                              dx,time,dt, &
                              the_bc_level, &
                              diff_coef,&
                              verbose,use_godunov_debug, &
                              use_minion)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: sedge(:)
      type(multifab) , intent(inout) :: ext_scal_force
!
      real(kind=dp_t), intent(inout) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: diff_coef
      logical        , intent(in)    :: use_godunov_debug
      logical        , intent(in)    :: use_minion
! 
      integer        , intent(in   ) :: lev,verbose
! 
      type(multifab) :: scal_force
! 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  dp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)
!
      type(multifab) :: divu
      integer :: nscal
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,comp,dm,ng_cell,ng_rho
      logical :: is_vel, make_divu
      logical, allocatable :: is_conservative(:)
      real(kind=dp_t) :: diff_fac
      real(kind=dp_t) :: half_dt

      half_dt = HALF * dt

      ng_cell = uold%ng
      ng_rho  = rhohalf%ng
      dm      = uold%dim

      nscal   = ncomp(sold)
      is_vel  = .false.

      allocate(is_conservative(nscal))
      is_conservative(1) = .true.
      is_conservative(2) = .false.

      call multifab_build(scal_force,ext_scal_force%la,nscal,1)
      call multifab_build(divu,scal_force%la,1,1)

      call setval(divu,0.0_dp_t,all=.true.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff_fac = ONE
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(scal_force, i)
          ep => dataptr(ext_scal_force, i)
         sop => dataptr(sold , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sop(:,:,1,:), &
                                  ng_cell, dx, &
                                  the_bc_level%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                  diff_coef, diff_fac)
            case (3)
              call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sop(:,:,:,:), &
                                  ng_cell, dx, &
                                  the_bc_level%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                  diff_coef, diff_fac)
         end select
      end do

      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of scalar using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         uop  => dataptr(uold, i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
          fp  => dataptr(scal_force , i)
          dp  => dataptr(divu, i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
         case (2)
            if(use_godunov_debug) then
               call mkflux_debug_2d(sop(:,:,1,:), uop(:,:,1,:), &
                                    sepx(:,:,1,:), sepy(:,:,1,:), &
                                    ump(:,:,1,1), vmp(:,:,1,1), &
                                    fp(:,:,1,:), dp(:,:,1,1), &
                                    lo, dx, dt, is_vel, &
                                    the_bc_level%phys_bc_level_array(i,:,:), &
                                    the_bc_level%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                              sepx(:,:,1,:), sepy(:,:,1,:), &
                              ump(:,:,1,1), vmp(:,:,1,1), &
                              fp(:,:,1,:), dp(:,:,1,1), &
                              lo, dx, dt, is_vel, &
                              the_bc_level%phys_bc_level_array(i,:,:), &
                              the_bc_level%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                              ng_cell, use_minion, is_conservative)
            endif
         case (3)
            sepz => dataptr(sedge(3), i)
            wmp  => dataptr(umac(3), i)
            if(use_godunov_debug) then
               call mkflux_debug_3d(sop(:,:,:,:), uop(:,:,:,:), &
                                    sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                    ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                    fp(:,:,:,:), dp(:,:,:,1), &
                                    lo, dx, dt, is_vel, &
                                    the_bc_level%phys_bc_level_array(i,:,:), &
                                    the_bc_level%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                    ng_cell, use_minion, is_conservative)
            else
               call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                              sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                              ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                              fp(:,:,:,:), dp(:,:,:,1), &
                              lo, dx, dt, is_vel, &
                              the_bc_level%phys_bc_level_array(i,:,:), &
                              the_bc_level%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                              ng_cell, use_minion, is_conservative)
            endif
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n+1/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff_fac = HALF
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(scal_force, i)
          ep => dataptr(ext_scal_force, i)
          sop => dataptr(sold , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sop(:,:,1,:), &
                                  ng_cell, dx, the_bc_level%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                  diff_coef, diff_fac)
            case (3)
              call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sop(:,:,:,:), &
                                  ng_cell, dx, the_bc_level%ell_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                  diff_coef, diff_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the scalars with conservative or convective differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         sop => dataptr(sold, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         snp => dataptr(snew, i)
          rp => dataptr(rhohalf, i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_2d(lev,sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             fp(:,:,1,:) , snp(:,:,1,:), &
                             rp(:,:,1,1) , &
                             lo, hi, ng_cell,dx,dt,is_vel,is_conservative,verbose)
            case (3)
               wmp => dataptr(umac(3), i)
               sepz => dataptr(sedge(3), i)
              call update_3d(lev,sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                             sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                             fp(:,:,:,:) , snp(:,:,:,:), &
                             rp(:,:,:,1) , &
                             lo, hi, ng_cell,dx,dt,is_vel,is_conservative,verbose)
         end select
      end do

      call multifab_fill_boundary(rhohalf)
      call multifab_fill_boundary(snew)

      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         snp => dataptr(snew, i)
          rp => dataptr(rhohalf, i)
         lo =  lwb(get_box(uold, i))
         select case (dm)
            case (2)
              do n = 1,nscal
                call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
              call setbc_2d(rp(:,:,1,1), lo, ng_rho, &
                            the_bc_level%adv_bc_level_array(i,:,:,dm+1),dx,dm+1)
            case (3)
              do n = 1,nscal
                call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
              call setbc_3d(rp(:,:,:,1), lo, ng_rho, &
                            the_bc_level%adv_bc_level_array(i,:,:,dm+1),dx,dm+1)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(is_conservative)
      call multifab_destroy(scal_force)
      call multifab_destroy(divu)

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
