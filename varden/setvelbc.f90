module setvelbc_module

  use bl_types
  use bc_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: EIGHTH = 0.125_dp_t

contains

      subroutine setvelbc_2d(vel,lo,ng,bc,irz)

      integer        , intent(in   ) :: lo(2),ng
      real(kind=dp_t), intent(inout) :: vel(lo(1)-ng:, lo(2)-ng:,:)
      integer, intent(in) ::  bc(2,2)
      integer, intent(in) ::  irz

!     Local variables
      integer :: i,j,hi(2)

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng+1)

      print *,'BC ',bc

      if (bc(1,1) .eq. WALL) then
         do j = lo(2),hi(2)
            vel(lo(1)-ng:lo(1)-1,j,1) = ZERO
            vel(lo(1)-ng:lo(1)-1,j,2) = &
              ( 15.0_dp_t * vel(lo(1)  ,j,2) &
               -10.0_dp_t * vel(lo(1)+1,j,2) &
               + 3.0_dp_t * vel(lo(1)+2,j,2) ) * EIGHTH
         end do
      end if
      if (bc(1,2) .eq. WALL) then
         do j = lo(2),hi(2)
            vel(hi(1)+1:hi(1)+ng,j,1) = ZERO
            vel(hi(1)+1:hi(1)+ng,j,2) = &
              ( 15.0_dp_t * vel(hi(1)  ,j,2) &
               -10.0_dp_t * vel(hi(1)-1,j,2) &
               + 3.0_dp_t * vel(hi(1)-2,j,2) ) * EIGHTH
         end do
      end if
      if (bc(2,1) .eq. WALL) then
         do i = lo(1)-ng,hi(1)+ng
            vel(i,lo(2)-ng:lo(2)-1,2) = ZERO
            vel(i,lo(2)-ng:lo(2)-1,1) = &
              ( 15.0_dp_t * vel(i,lo(2)  ,1) &
               -10.0_dp_t * vel(i,lo(2)+1,1) &
               + 3.0_dp_t * vel(i,lo(2)+2,1) ) * EIGHTH
         end do
      end if
      if (bc(2,2) .eq. WALL) then
         do i = lo(1)-ng,hi(1)+ng
            vel(i,hi(2)+1:hi(2)+ng,2) = ZERO
            vel(i,hi(2)+1:hi(2)+ng,1) = &
              ( 15.0_dp_t * vel(i,hi(2)  ,1) &
               -10.0_dp_t * vel(i,hi(2)-1,1) &
               + 3.0_dp_t * vel(i,hi(2)-2,1) ) * EIGHTH
         end do
      end if

      end subroutine setvelbc_2d

      subroutine setvelbc_3d(vel,lo,ng,bc)

      integer        , intent(in   ) :: lo(3),ng
      real(kind=dp_t), intent(inout) :: vel(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:, :)
      integer, intent(in) ::  bc(3,2)

!     Local variables
      integer :: i,j,k,hi(3)

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(vel,dim=3) - (2*ng+1)

      if (bc(1,1) .eq. WALL) then
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            vel(lo(1)-ng:lo(1)-1,j,k,1) = ZERO
            vel(lo(1)-ng:lo(1)-1,j,k,2) = &
              ( 15.0_dp_t * vel(lo(1)  ,j,k,2) &
               -10.0_dp_t * vel(lo(1)+1,j,k,2) &
               + 3.0_dp_t * vel(lo(1)+2,j,k,2) ) * EIGHTH
            vel(lo(1)-ng:lo(1)-1,j,k,3) = &
              ( 15.0_dp_t * vel(lo(1)  ,j,k,3) &
               -10.0_dp_t * vel(lo(1)+1,j,k,3) &
               + 3.0_dp_t * vel(lo(1)+2,j,k,3) ) * EIGHTH
         end do
         end do
      end if
      if (bc(1,2) .eq. WALL) then
         do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            vel(hi(1)+1:hi(1)+ng,j,k,1) = ZERO
            vel(hi(1)+1:hi(1)+ng,j,k,2) = &
              ( 15.0_dp_t * vel(hi(1)  ,j,k,2) &
               -10.0_dp_t * vel(hi(1)-1,j,k,2) &
               + 3.0_dp_t * vel(hi(1)-2,j,k,2) ) * EIGHTH
            vel(hi(1)+1:hi(1)+ng,j,k,3) = &
              ( 15.0_dp_t * vel(hi(1)  ,j,k,3) &
               -10.0_dp_t * vel(hi(1)-1,j,k,3) &
               + 3.0_dp_t * vel(hi(1)-2,j,k,3) ) * EIGHTH
         end do
         end do
      end if
      if (bc(2,1) .eq. WALL) then
         do k = lo(3),hi(3)
         do i = lo(1)-ng,hi(1)+ng
            vel(i,lo(2)-ng:lo(2)-1,k,2) = ZERO
            vel(i,lo(2)-ng:lo(2)-1,k,1) = &
              ( 15.0_dp_t * vel(i,lo(2)  ,k,1) &
               -10.0_dp_t * vel(i,lo(2)+1,k,1) &
               + 3.0_dp_t * vel(i,lo(2)+2,k,1) ) * EIGHTH
            vel(i,lo(2)-ng:lo(2)-1,k,3) = &
              ( 15.0_dp_t * vel(i,lo(2)  ,k,3) &
               -10.0_dp_t * vel(i,lo(2)+1,k,3) &
               + 3.0_dp_t * vel(i,lo(2)+2,k,3) ) * EIGHTH
         end do
         end do
      end if
      if (bc(2,2) .eq. WALL) then
         do k = lo(3),hi(3)
         do i = lo(1)-ng,hi(1)+ng
            vel(i,hi(2)+1:hi(2)+ng,k,2) = ZERO
            vel(i,hi(2)+1:hi(2)+ng,k,1) = &
              ( 15.0_dp_t * vel(i,hi(2)  ,k,1) &
               -10.0_dp_t * vel(i,hi(2)-1,k,1) &
               + 3.0_dp_t * vel(i,hi(2)-2,k,1) ) * EIGHTH
            vel(i,hi(2)+1:hi(2)+ng,k,3) = &
              ( 15.0_dp_t * vel(i,hi(2)  ,k,3) &
               -10.0_dp_t * vel(i,hi(2)-1,k,3) &
               + 3.0_dp_t * vel(i,hi(2)-2,k,3) ) * EIGHTH
         end do
         end do
      end if
      if (bc(3,1) .eq. WALL) then
         do j = lo(2)-ng,hi(2)+ng
         do i = lo(1)-ng,hi(1)+ng
            vel(i,j,lo(3)-ng:lo(3)-1,3) = ZERO
            vel(i,j,lo(3)-ng:lo(3)-1,1) = &
              ( 15.0_dp_t * vel(i,j,lo(3)  ,1) &
               -10.0_dp_t * vel(i,j,lo(3)+1,1) &
               + 3.0_dp_t * vel(i,j,lo(3)+2,1) ) * EIGHTH
            vel(i,j,lo(3)-ng:lo(3)-1,2) = &
              ( 15.0_dp_t * vel(i,j,lo(3)  ,2) &
               -10.0_dp_t * vel(i,j,lo(3)+1,2) &
               + 3.0_dp_t * vel(i,j,lo(3)+2,2) ) * EIGHTH
         end do
         end do
      end if
      if (bc(3,2) .eq. WALL) then
         do j = lo(2)-ng,hi(2)+ng
         do i = lo(1)-ng,hi(1)+ng
            vel(i,j,hi(3)+1:hi(3)+ng,3) = ZERO
            vel(i,j,hi(3)+1:hi(3)+ng,1) = &
              ( 15.0_dp_t * vel(i,j,hi(3)  ,1) &
               -10.0_dp_t * vel(i,j,hi(3)-1,1) &
               + 3.0_dp_t * vel(i,j,hi(3)-2,1) ) * EIGHTH
            vel(i,j,hi(3)+1:hi(3)+ng,2) = &
              ( 15.0_dp_t * vel(i,j,hi(3)  ,2) &
               -10.0_dp_t * vel(i,j,hi(3)-1,2) &
               + 3.0_dp_t * vel(i,j,hi(3)-2,2) ) * EIGHTH
         end do
         end do
      end if

      end subroutine setvelbc_3d

end module setvelbc_module
