module pre_advance_module

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use velpred_module
  use mkforce_module
  use setbc_module

  implicit none

contains

   subroutine advance_premac(uold,sold,umac, &
                             gp,p,ext_force, &
                             dx,time,dt, &
                             the_bc_level, &
                             visc_coef,use_godunov_debug, &
                             use_minion)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(in   ) :: gp
      type(multifab) , intent(in   ) :: p
      type(multifab) , intent(in   ) :: ext_force
      real(kind=dp_t), intent(in   ) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: visc_coef
      logical        , intent(in)    :: use_godunov_debug
      logical        , intent(in)    :: use_minion
 
      type(multifab) :: force
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
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
      integer :: dm,ng_cell
      integer :: i,n,comp
      logical :: is_conservative(uold%dim)
      real(kind=dp_t) :: visc_fac, visc_mu

      real(kind=dp_t) :: half_dt

      half_dt = HALF * dt

      is_conservative = .false.

      ng_cell = uold%ng
      dm      = uold%dim

      call multifab_build(     force,ext_force%la,dm,1)

      call multifab_fill_boundary(uold)
      call multifab_fill_boundary(sold)

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

!     Create the edge states to be used for the MAC velocity 
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold  , i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
         case (2)
            if(use_godunov_debug) then
               call velpred_debug_2d(uop(:,:,1,:), &
                                     ump(:,:,1,1),  vmp(:,:,1,1), &
                                     fp(:,:,1,:), &
                                     lo, dx, dt, &
                                     the_bc_level%phys_bc_level_array(i,:,:), &
                                     the_bc_level%adv_bc_level_array(i,:,:,:), &
                                     ng_cell, use_minion)
            else
               call velpred_2d(uop(:,:,1,:), &
                               ump(:,:,1,1),  vmp(:,:,1,1), &
                               fp(:,:,1,:), &
                               lo, dx, dt, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,:), &
                               ng_cell, use_minion)
            endif
         case (3)
            wmp  => dataptr(umac(3), i)
            if(use_godunov_debug) then
               call velpred_debug_3d(uop(:,:,:,:), &
                                     ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                                     fp(:,:,:,:), &
                                     lo, dx, dt, &
                                     the_bc_level%phys_bc_level_array(i,:,:), &
                                     the_bc_level%adv_bc_level_array(i,:,:,:), &
                                     ng_cell, use_minion)
            else
               call velpred_3d(uop(:,:,:,:), &
                               ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                               fp(:,:,:,:), &
                               lo, dx, dt, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,:), &
                               ng_cell, use_minion)
            endif
         end select
      end do

      call multifab_destroy(force)

   end subroutine advance_premac

end module pre_advance_module
