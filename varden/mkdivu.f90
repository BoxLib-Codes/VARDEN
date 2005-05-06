module mkdivu_module

  use bl_types
  implicit none

contains

      subroutine mkdivu_2d(divu,umac,dx)

      real(kind=dp_t), intent(in   ) :: umac(:,:,:)
      real(kind=dp_t), intent(  out) :: divu(:,:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j

      do j = 1,size(divu,dim=2)
      do i = 1,size(divu,dim=1)
        divu(i,j) = (umac(i+1,j,1) - umac(i,j,1)) / dx(1) &
                   +(umac(i,j+1,2) - umac(i,j,2)) / dx(2)
      end do
      end do

      end subroutine mkdivu_2d

      subroutine mkdivu_3d(divu,umac,dx)

      real(kind=dp_t), intent(in   ) :: umac(:,:,:,:)
      real(kind=dp_t), intent(  out) :: divu(:,:,:)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j,k

      do j = 1,size(divu,dim=2)
      do i = 1,size(divu,dim=1)
        divu(i,j,k) = (umac(i+1,j,k,1) - umac(i,j,k,1)) / dx(1) &
                     +(umac(i,j+1,k,2) - umac(i,j,k,2)) / dx(2)
                     +(umac(i,j,k+1,3) - umac(i,j,k,2)) / dx(3)
      end do
      end do

      end subroutine mkdivu_3d

end module mkdivu_module
