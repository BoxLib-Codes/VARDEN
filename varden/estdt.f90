module estdt_module

  use bl_types
  use multifab_module

  implicit none

contains

   subroutine estdt (u,s,dx,cflfac,dtold,dt)

      type(multifab) , intent( in) :: u,s
      real(kind=dp_t), intent( in) :: dx(:)
      real(kind=dp_t), intent( in) :: cflfac,dtold
      real(kind=dp_t), intent(out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_hold
      integer :: i

      ng = u%ng
      dm = u%dim

      dt_hold = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, &
                            ng, dx, cflfac, dtold, dt)
            case (3)
              call estdt_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, &
                            ng, dx, cflfac, dtold, dt)
         end select
         dt = min(dt_hold,dt)
      end do

   end subroutine estdt

   subroutine estdt_2d (u,s,lo,hi,ng,dx,cflfac,dtold,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) :: u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: cflfac
      real (kind = dp_t), intent(in ) :: dtold
      real (kind = dp_t), intent(out) :: dt
    
!     real (kind = dp_t)     gp(lo_1-1:hi_1+1,lo_2-1:hi_2+1,2)
!     real (kind = dp_t)  force(lo_1-1:hi_1+1,lo_2-1:hi_2+1,2)

!     Local variables
      real (kind = dp_t)  spdx,spdy
      real (kind = dp_t)  pforcex,pforcey
      real (kind = dp_t)  dtchange
      real (kind = dp_t)  eps
      integer :: i, j

      dtchange = 1.1d0
      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,1))/dx(1))
          spdy    = max(spdy ,abs(u(i,j,2))/dx(2))
!         pforcex = max(pforcex,abs(gp(i,j,1)/s(i,j)-force(i,j,1)))
!         pforcey = max(pforcey,abs(gp(i,j,2)/s(i,j)-force(i,j,2)))
        enddo
      enddo

      if (spdx.lt.eps .and. spdy.lt.eps) then

        dt = min(dx(1),dx(2))

      else

        dt = 1.0D0  / max(spdx,spdy)

      endif

      if (pforcex .gt. eps) then
        dt = min(dt,sqrt(2.0D0 *dx(1)/pforcex))
      endif

      if (pforcey .gt. eps) then
        dt = min(dt,sqrt(2.0D0 *dx(2)/pforcey))
      endif

      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

   end subroutine estdt_2d

   subroutine estdt_3d (u,s,lo,hi,ng,dx,cflfac,dtold,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: cflfac
      real (kind = dp_t), intent(in ) :: dtold
      real (kind = dp_t), intent(out) :: dt
    
!     real (kind = dp_t)     gp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1, &
!                               lo(3)-1:hi(3)+1,3)
!     real (kind = dp_t)  force(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1, &
!                               lo(3)-1:hi(3)+1,3)

!     Local variables
      real (kind = dp_t)  spdx,spdy
      real (kind = dp_t)  pforcex,pforcey,pforcez
      real (kind = dp_t)  dtchange
      real (kind = dp_t)  eps
      integer :: i, j, k

      dtchange = 1.1d0
      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 
      pforcez = 0.0D0 

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,k,1))/dx(1))
          spdy    = max(spdy ,abs(u(i,j,k,2))/dx(2))
!         pforcex = max(pforcex,abs(gp(i,j,k,1)/s(i,j,k)-force(i,j,k,1)))
!         pforcey = max(pforcey,abs(gp(i,j,k,2)/s(i,j,k)-force(i,j,k,2)))
        enddo
      enddo
      enddo

      if (spdx.lt.eps .and. spdy.lt.eps) then

        dt = min(dx(1),dx(2))

      else

        dt = 1.0D0  / max(spdx,spdy)

      endif

!     if (pforcex .gt. eps) then
!       dt = min(dt,sqrt(2.0D0 *dx(1)/pforcex))
!     endif

!     if (pforcey .gt. eps) then
!       dt = min(dt,sqrt(2.0D0 *dx(2)/pforcey))
!     endif

!     if (pforcez .gt. eps) then
!       dt = min(dt,sqrt(2.0D0 *dx(3)/pforcez))
!     endif

      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

   end subroutine estdt_3d

end module estdt_module
