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

   subroutine advance_premac(nlevs,uold,sold,umac, &
                             gp,ext_vel_force, &
                             dx,time,dt, &
                             the_bc_level, &
                             visc_coef,use_godunov_debug, &
                             use_minion)
 
      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: uold(:)
      type(multifab) , intent(inout) :: sold(:)
      type(multifab) , intent(inout) :: umac(:,:)
      type(multifab) , intent(in   ) :: gp(:)
      type(multifab) , intent(in   ) :: ext_vel_force(:)
      real(kind=dp_t), intent(in   ) :: dx(:,:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level(:)
      real(kind=dp_t), intent(in   ) :: visc_coef
      logical        , intent(in)    :: use_godunov_debug
      logical        , intent(in)    :: use_minion
 
      type(multifab), allocatable :: vel_force(:)
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)

      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  up(:,:,:,:)
      real(kind=dp_t), pointer::  sp(:,:,:,:)
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)

      integer :: lo(uold(1)%dim),hi(uold(1)%dim)
      integer :: dm,ng_cell
      integer :: i,n,d
      logical :: is_conservative(uold(1)%dim)
      real(kind=dp_t) :: visc_fac, visc_mu
      real(kind=dp_t) :: half_dt

      allocate(vel_force(nlevs))

      half_dt = HALF * dt

      is_conservative = .false.

      ng_cell = uold(1)%ng
      dm      = uold(1)%dim

      do n = 1, nlevs
         call multifab_build(vel_force(n),ext_vel_force(n)%la,dm,1)

         call setval(vel_force(n),0.0_dp_t,all=.true.)

         call multifab_fill_boundary(uold(n))
         call multifab_fill_boundary(sold(n))
      enddo

      do n = 1, nlevs

      ! Impose boundary conditions on sold and uold.
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         uop => dataptr(uold(n), i)
         sop => dataptr(sold(n), i)
         lo =  lwb(get_box(uold(n), i))
         hi =  upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            do d = 1, dm
               call setbc_2d(uop(:,:,1,d), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,d),dx(n,:),d)
            end do
            d = dm+1
            call setbc_2d(sop(:,:,1,1), lo, ng_cell, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,d),dx(n,:),d)
         case (3)
            do d = 1, dm
               call setbc_3d(uop(:,:,:,d), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,d),dx(n,:),d)
            end do
            d = dm+1
            call setbc_3d(sop(:,:,:,1), lo, ng_cell, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,d),dx(n,:),d)
         end select
      end do
      
      ! Create vel_force
      visc_fac = ONE
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         fp  => dataptr(vel_force(n), i)
         ep  => dataptr(ext_vel_force(n), i)
         gpp => dataptr(gp(n), i)
         rp  => dataptr(sold(n), i)
         up  => dataptr(uold(n), i)
         lo = lwb(get_box(uold(n), i))
         hi = upb(get_box(uold(n), i))
         select case (dm)
            case (2)
              call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), &
                                 gpp(:,:,1,:), rp(:,:,1,1), up(:,:,1,:), &
                                 ng_cell, ng_cell, dx(n,:), &
                                 the_bc_level(n)%ell_bc_level_array(i,:,:,:), &
                                 visc_coef, visc_fac)
            case (3)
              call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), &
                                 gpp(:,:,:,:), rp(:,:,:,1), up(:,:,:,:), &
                                 ng_cell, ng_cell, dx(n,:), &
                                 the_bc_level(n)%ell_bc_level_array(i,:,:,:), &
                                 visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(vel_force(n))

      ! Create the edge states to be used for the MAC velocity 
      do i = 1, uold(n)%nboxes
         if ( multifab_remote(uold(n), i) ) cycle
         uop => dataptr(uold(n), i)
         ump => dataptr(umac(n,1), i)
         vmp => dataptr(umac(n,2), i)
         fp => dataptr(vel_force(n) , i)
         lo =  lwb(get_box(uold(n), i))
         hi =  upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            if(use_godunov_debug) then
               call velpred_debug_2d(uop(:,:,1,:), &
                                     ump(:,:,1,1),  vmp(:,:,1,1), &
                                     fp(:,:,1,:), &
                                     lo, dx(n,:), dt, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                     ng_cell, use_minion)
            else
               call velpred_2d(uop(:,:,1,:), &
                               ump(:,:,1,1),  vmp(:,:,1,1), &
                               fp(:,:,1,:), &
                               lo, dx(n,:), dt, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                               ng_cell, use_minion)
            endif
         case (3)
            wmp  => dataptr(umac(n,3), i)
            if(use_godunov_debug) then
               call velpred_debug_3d(uop(:,:,:,:), &
                                     ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                                     fp(:,:,:,:), &
                                     lo, dx(n,:), dt, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                     ng_cell, use_minion)
            else
               call velpred_3d(uop(:,:,:,:), &
                               ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                               fp(:,:,:,:), &
                               lo, dx(n,:), dt, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                               ng_cell, use_minion)
            endif
         end select
      end do
      
      do d = 1, dm
         call multifab_fill_boundary(umac(n,d))
      enddo
      
      enddo ! do n = 1, nlevs

      do n = 1, nlevs
         call multifab_destroy(vel_force(n))
      enddo
      
    end subroutine advance_premac

end module pre_advance_module
