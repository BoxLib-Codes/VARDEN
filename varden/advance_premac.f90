module pre_advance_module

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use mkflux_module
  use mkutrans_module
  use mkforce_module
  use setbc_module

  implicit none

contains

   subroutine advance_premac(uold,sold,umac,uedge,utrans, &
                             gp,p,ext_force, &
                             dx,time,dt, &
                             the_bc_level, &
                             visc_coef)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: uedge(:)
      type(multifab) , intent(inout) :: utrans(:)
      type(multifab) , intent(in   ) :: gp
      type(multifab) , intent(in   ) :: p
      type(multifab) , intent(in   ) :: ext_force
      real(kind=dp_t), intent(in   ) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: visc_coef
 
      type(multifab) :: force
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: vtp(:,:,:,:)
      real(kind=dp_t), pointer:: wtp(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer:: uepx(:,:,:,:)
      real(kind=dp_t), pointer:: uepy(:,:,:,:)
      real(kind=dp_t), pointer:: uepz(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  up(:,:,:,:)
      real(kind=dp_t), pointer::  sp(:,:,:,:)
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
!
      integer :: edge_based
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: velpred
      integer :: dm,ng_cell
      integer :: i,n,comp
      logical :: is_conservative(uold%dim)
      logical :: is_vel
      real(kind=dp_t) :: visc_fac, visc_mu

      real(kind=dp_t) :: half_dt

      half_dt = HALF * dt

      is_conservative = .false.

      ng_cell = uold%ng
      dm      = uold%dim

      call multifab_build(     force,ext_force%la,dm,1)

!     Impose boundary conditions on sold and uold.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold , i)
         sop => dataptr(sold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              do n = 1,dm
                call setbc_2d(uop(:,:,1,n), lo, ng_cell, the_bc_level%adv_bc_level_array(i,:,:,n),dx,n)
              end do
              n = dm+1
                call setbc_2d(sop(:,:,1,1), lo, ng_cell, the_bc_level%adv_bc_level_array(i,:,:,n),dx,n)
            case (3)
              do n = 1,dm
                call setbc_3d(uop(:,:,:,n), lo, ng_cell, the_bc_level%adv_bc_level_array(i,:,:,n),dx,n)
              end do
              n = dm+1
                call setbc_3d(sop(:,:,:,1), lo, ng_cell, the_bc_level%adv_bc_level_array(i,:,:,n),dx,n)
         end select
      end do
      call multifab_fill_boundary(uold)
      call multifab_fill_boundary(sold)

!     Create force.
      visc_fac = ONE
      edge_based = 0
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(force, i)
          ep => dataptr(ext_force, i)
         gpp => dataptr(gp   , i)
          rp => dataptr(sold , i)
          up => dataptr(uold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         uop     => dataptr( uold, i)
         select case (dm)
            case (2)
              call mkforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                              gpp(:,:,1,:), rp(:,:,1,1), up(:,:,1,:), &
                              ng_cell, ng_cell, dx, &
                              the_bc_level%ell_bc_level_array(i,:,:,:), &
                              visc_coef, visc_fac)
            case (3)
              call mkforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                              gpp(:,:,:,:), rp(:,:,:,1), up(:,:,:,:), &
                              ng_cell, ng_cell, dx, &
                              the_bc_level%ell_bc_level_array(i,:,:,:), &
                              visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(force)

!     Create utrans.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold  , i)
         utp => dataptr(utrans(1), i)
         vtp => dataptr(utrans(2), i)
          fp => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkutrans_2d(uop(:,:,1,:), utp(:,:,1,1), vtp(:,:,1,1), &
                               fp(:,:,1,:), &
                               lo,dx,dt,ng_cell,&
                               the_bc_level%adv_bc_level_array(i,:,:,:), &
                               the_bc_level%phys_bc_level_array(i,:,:))
            case (3)
               wtp => dataptr(utrans(3), i)
              call mkutrans_3d(uop(:,:,:,:), utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), &
                               fp(:,:,:,:), &
                               lo,dx,dt,ng_cell,&
                               the_bc_level%adv_bc_level_array(i,:,:,:), &
                               the_bc_level%phys_bc_level_array(i,:,:))
         end select
      end do

!     Create the edge states to be used for the MAC velocity 
      velpred = 1
      is_vel = .true.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold  , i)
         uepx => dataptr(uedge(1), i)
         uepy => dataptr(uedge(2), i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         utp  => dataptr(utrans(1), i)
         vtp  => dataptr(utrans(2), i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             ump(:,:,1,1),  vmp(:,:,1,1), &
                             utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, is_conservative, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,:), &
                             velpred, ng_cell)
            case (3)
               uepz => dataptr(uedge(3), i)
               wmp  => dataptr(umac(3), i)
               wtp  => dataptr(utrans(3), i)
              call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                             uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                             ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                             utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                             lo, dx, dt, is_vel, is_conservative, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,:), &
                             velpred, ng_cell)
         end select
      end do

      call multifab_destroy(force)

   end subroutine advance_premac

end module pre_advance_module
