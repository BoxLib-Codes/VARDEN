module velocity_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module

  implicit none

contains

   subroutine velocity_advance(uold,unew,sold,rhohalf,&
                              umac,uedgex,uedgey, &
                              utrans,gp,p, &
                              ext_force, &
                              dx,time,dt, &
                              the_bc_level, &
                              visc_coef,verbose,mg_verbose)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: umac
      type(multifab) , intent(inout) :: unew
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: uedgex,uedgey,utrans
      type(multifab) , intent(inout) :: gp
      type(multifab) , intent(inout) :: p
      type(multifab) , intent(inout) :: ext_force

      real(kind=dp_t), intent(inout) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: visc_coef
 
      integer        , intent(in   ) :: verbose, mg_verbose
 
      type(multifab) :: force
 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: unp(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer:: uepx(:,:,:,:)
      real(kind=dp_t), pointer:: uepy(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
!
      integer :: irz,velpred,edge_based
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,comp,dm,ng_cell,ng_edge,ng_rho
      logical :: is_vel,is_conservative(uold%dim)
      real(kind=dp_t) :: visc_fac, visc_mu
      real(kind=dp_t) :: half_dt

      half_dt = HALF * dt

      ng_cell = uold%ng
      ng_rho  = rhohalf%ng
      ng_edge = umac%ng
      dm      = uold%dim

      is_conservative = .false.

      irz = 0

      call multifab_build(     force,ext_force%la,   dm,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the velocity forcing term at time n using rho and the full viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      visc_fac = ONE
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(force, i)
          ep => dataptr(ext_force, i)
         gpp => dataptr(gp   , i)
          rp => dataptr(sold , i)
         uop => dataptr(uold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         uop     => dataptr( uold, i)
         select case (dm)
            case (2)
              call mkforce_2d( fp(:,:,1,:), ep(:,:,1,:), &
                              gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                              ng_cell, ng_cell, dx, &
                              the_bc_level%ell_bc_level_array(i,:,:,1:dm), &
                              visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(force)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of velocity using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      velpred = 0
      is_vel = .true.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold, i)
         uepx => dataptr(uedgex, i)
         uepy => dataptr(uedgey, i)
         ump  => dataptr(umac, i)
         utp  => dataptr(utrans, i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             ump(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, is_conservative, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%ell_bc_level_array(i,:,:,1:dm), &
                             velpred, ng_cell, ng_edge)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Now create the force at half-time using rhohalf and half the viscous term.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      visc_fac = HALF
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(force, i)
          ep => dataptr(ext_force, i)
         gpp => dataptr(gp   , i)
          rp => dataptr(rhohalf , i)
         uop => dataptr(uold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkforce_2d( fp(:,:,1,:), ep(:,:,1,:), &
                              gpp(:,:,1,:), rp(:,:,1,1), uop(:,:,1,:), &
                              ng_cell, ng_rho, dx, &
                              the_bc_level%ell_bc_level_array(i,:,:,1:dm), &
                              visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update the velocity with convective differencing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold, i)
         ump => dataptr(umac, i)
         uepx => dataptr(uedgex, i)
         uepy => dataptr(uedgey, i)
         unp => dataptr(unew, i)
          fp => dataptr(force, i)
          rp => dataptr(rhohalf, i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_2d(uop(:,:,1,:), ump(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             fp(:,:,1,:), unp(:,:,1,:), &
                             rp(:,:,1,1), &
                             lo, hi, ng_cell, ng_edge, dx, time, dt, is_vel, is_conservative)
         end select
      end do

      call multifab_destroy(force)

   end subroutine velocity_advance

end module velocity_advance_module
