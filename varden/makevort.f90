module vort_module

  use bl_types
  use bc_module
  use define_bc_module
  use multifab_module

  implicit none

contains

   subroutine make_vorticity (vort,u,dx,bc)

      type(multifab) , intent(inout) :: vort
      type(multifab) , intent(in   ) :: u
      real(kind=dp_t), intent(in   ) :: dx(:)
      type(bc_level) , intent(in   ) :: bc

      real(kind=dp_t), pointer:: up(:,:,:,:)
      real(kind=dp_t), pointer:: vp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      integer :: i

      ng = u%ng
      dm = u%dim

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         up => dataptr(u, i)
         vp => dataptr(vort, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call makevort_2d(vp(:,:,1,1),up(:,:,1,:), lo, hi, ng, dx, bc%phys_bc_level_array(i,:,:))
            case (3)
               call makevort_3d(vp(:,:,:,1),up(:,:,:,:), lo, hi, ng, dx, bc%phys_bc_level_array(i,:,:))
               call bl_warn('3d vorticitiy not implemented')
         end select
      end do

   end subroutine make_vorticity

   subroutine makevort_2d (vort,u,lo,hi,ng,dx,bc)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: vort(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      integer           , intent(in   ) :: bc(:,:)

!     Local variables
      integer :: i, j, n
      real (kind = dp_t) :: vx,uy

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        enddo
      enddo

      if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. &
          bc(1,1) .eq. NO_SLIP_WALL) then
        i = lo(1)
        do j = lo(2), hi(2)
           vx = (u(i+1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i-1,j,2)) / dx(1)
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. &
          bc(1,2) .eq. NO_SLIP_WALL) then
        i = hi(1)
        do j = lo(2), hi(2)
           vx = -(u(i-1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i+1,j,2)) / dx(1)
           uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. &
          bc(2,1) .eq. NO_SLIP_WALL) then
        j = lo(2)
        do i = lo(1), hi(1)
           vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = (u(i,j+1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j-1,1)) / dx(2)
           vort(i,j) = vx - uy
        end do
      end if

      if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. &
          bc(2,2) .eq. NO_SLIP_WALL) then
        j = hi(2)
        do i = lo(1), hi(1)
           vx =  (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
           uy = -(u(i,j-1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j+1,1)) / dx(2)
           vort(i,j) = vx - uy
        end do
      end if

   end subroutine makevort_2d

   subroutine makevort_3d (vort,u,lo,hi,ng,dx,bc)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: vort(lo(1):,lo(2):,lo(3):)
      real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      integer           , intent(in   ) :: bc(:,:)

      vort = 0.0_dp_t

   end subroutine makevort_3d


end module vort_module
