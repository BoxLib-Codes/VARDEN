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

   subroutine scalar_advance (uold,sold,snew,rhohalf,&
                              umac,sedgex,sedgey, &
                              utrans, &
                              ext_scal_force,&
                              dx,time,dt, &
                              the_bc_level, &
                              diff_coef,&
                              verbose,mg_verbose)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: umac
      type(multifab) , intent(inout) :: sedgex,sedgey
      type(multifab) , intent(inout) :: utrans
      type(multifab) , intent(inout) :: ext_scal_force
!
      real(kind=dp_t), intent(inout) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
!     integer        , intent(in   ) :: domain_scal_bc(:,:)
!     type(bc_level) , intent(in   ) :: scal_bc
      real(kind=dp_t), intent(in   ) :: diff_coef
! 
      integer        , intent(in   ) :: verbose, mg_verbose
! 
      type(multifab) :: force,scal_force
! 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
!
      integer :: nscal,velpred
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,comp,dm,ng_cell,ng_edge,ng_rho
      logical :: is_vel,is_conservative
      real(kind=dp_t) :: visc_fac, visc_mu
      real(kind=dp_t) :: half_dt

      half_dt = HALF * dt

      ng_cell = uold%ng
      ng_rho  = rhohalf%ng
      ng_edge = umac%ng
      dm      = uold%dim

      nscal   = 2
      is_vel = .false.

      call multifab_build(scal_force,ext_scal_force%la,nscal,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      visc_fac = ONE
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
                                  diff_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of scalar using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      velpred = 0
      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         uop  => dataptr(uold, i)
         sepx => dataptr(sedgex, i)
         sepy => dataptr(sedgey, i)
         ump  => dataptr(umac, i)
         utp  => dataptr(utrans, i)
          fp  => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             ump(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                             velpred, ng_cell, ng_edge)
            case (3)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n+1/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      visc_fac = HALF
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
                                  diff_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the scalars with conservative differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      is_conservative = .true.
      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         sop => dataptr(sold, i)
         ump => dataptr(umac, i)
         sepx => dataptr(sedgex, i)
         sepy => dataptr(sedgey, i)
         snp => dataptr(snew, i)
          rp => dataptr(rhohalf, i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_2d(sop(:,:,1,:), ump(:,:,1,:), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             fp(:,:,1,:) , snp(:,:,1,:), &
                             rp(:,:,1,1) , &
                             lo, hi, ng_cell, ng_edge, dx, time, dt, is_conservative)
              do n = 1,nscal
                call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
              call setbc_2d(rp(:,:,1,1), lo, ng_rho, &
                            the_bc_level%adv_bc_level_array(i,:,:,dm+1),dx,dm+1)
         end select
      end do
      call multifab_fill_boundary(rhohalf)
      call multifab_fill_boundary(snew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call multifab_destroy(scal_force)

   end subroutine scalar_advance

end module scalar_advance_module
