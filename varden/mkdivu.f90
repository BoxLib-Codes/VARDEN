module mkdivu_module

  use bl_types
  implicit none

contains

      subroutine mkdivu_2d(divu,umac,vmac,dx)

      real(kind=dp_t), intent(in   ) :: umac(0:,0:)
      real(kind=dp_t), intent(in   ) :: vmac(0:,0:)
      real(kind=dp_t), intent(  out) :: divu(:,:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j

      do j = 1,size(divu,dim=2)
      do i = 1,size(divu,dim=1)
        divu(i,j) = (umac(i+1,j) - umac(i,j)) / dx(1) &
                   +(vmac(i,j+1) - vmac(i,j)) / dx(2)
      end do
      end do

      end subroutine mkdivu_2d

      subroutine mkdivu_3d(divu,umac,vmac,wmac,dx)

      real(kind=dp_t), intent(in   ) :: umac(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: vmac(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: wmac(0:,0:,0:)
      real(kind=dp_t), intent(  out) :: divu(:,:,:)
      real(kind=dp_t), intent(in   ) ::   dx(:)


      integer :: i,j,k

      do k = 1,size(divu,dim=3)
      do j = 1,size(divu,dim=2)
      do i = 1,size(divu,dim=1)
        divu(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                     +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                     +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
      end do
      end do
      end do

      end subroutine mkdivu_3d

end module mkdivu_module
