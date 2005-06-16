module estdt_module

  use bl_types
  use multifab_module

  implicit none

contains

   subroutine estdt (u,s,gp,force,dx,cflfac,dtold,dt)

      type(multifab) , intent( in) :: u,s,gp,force
      real(kind=dp_t), intent( in) :: dx(:)
      real(kind=dp_t), intent( in) :: cflfac,dtold
      real(kind=dp_t), intent(out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:),  fp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_hold
      real(kind=dp_t) :: dtchange
      integer         :: i

      ng = u%ng
      dm = u%dim

      dtchange = 1.1d0
      dt_hold  = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         gpp => dataptr(gp, i)
          fp => dataptr(force, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,1), gpp(:,:,1,:), fp(:,:,1,:),&
                            lo, hi, ng, dx, dt)
            case (3)
              call estdt_3d(uop(:,:,:,:), sop(:,:,:,1), gpp(:,:,:,:), fp(:,:,:,:),&
                            lo, hi, ng, dx, dt)
         end select
         dt_hold = min(dt_hold,dt)
      end do

      dt = dt_hold

      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

   end subroutine estdt

   subroutine estdt_2d (u,s,gp,force,lo,hi,ng,dx,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in ) ::    gp(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  spdx,spdy
      real (kind = dp_t)  pforcex,pforcey
      real (kind = dp_t)  eps
      integer :: i, j

      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,1))/dx(1))
          spdy    = max(spdy ,abs(u(i,j,2))/dx(2))
          pforcex = max(pforcex,abs(gp(i,j,1)/s(i,j)-force(i,j,1)))
          pforcey = max(pforcey,abs(gp(i,j,2)/s(i,j)-force(i,j,2)))
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

   end subroutine estdt_2d

   subroutine estdt_3d (u,s,gp,force,lo,hi,ng,dx,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(in ) ::    gp(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  spdx,spdy
      real (kind = dp_t)  pforcex,pforcey,pforcez
      real (kind = dp_t)  eps
      integer :: i, j, k

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
          pforcex = max(pforcex,abs(gp(i,j,k,1)/s(i,j,k)-force(i,j,k,1)))
          pforcey = max(pforcey,abs(gp(i,j,k,2)/s(i,j,k)-force(i,j,k,2)))
        enddo
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

      if (pforcez .gt. eps) then
        dt = min(dt,sqrt(2.0D0 *dx(3)/pforcez))
      endif

   end subroutine estdt_3d

end module estdt_module
