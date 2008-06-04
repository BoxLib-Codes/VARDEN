module estdt_module

  use bl_types
  use multifab_module
  use probin_module, only : verbose

  implicit none

  private

  public :: estdt

contains

   subroutine estdt (lev,u, s, gp, ext_vel_force, dx, dtold, dt)

      use probin_module, only: max_dt_growth, cflfac

      type(multifab) , intent( in) :: u,s,gp,ext_vel_force
      real(kind=dp_t), intent( in) :: dx(:)
      real(kind=dp_t), intent( in) :: dtold
      real(kind=dp_t), intent(out) :: dt
      integer        , intent( in) :: lev

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:),  fp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_proc, dt_grid, dt_start
      integer         :: i

      ng = u%ng
      dm = u%dim

      dt_proc  = 1.d20
      dt_start = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         gpp => dataptr(gp, i)
          fp => dataptr(ext_vel_force, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))

         dt_grid = 1.d20
         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,1), gpp(:,:,1,:), fp(:,:,1,:),&
                            lo, hi, ng, dx, dt_grid)
            case (3)
              call estdt_3d(uop(:,:,:,:), sop(:,:,:,1), gpp(:,:,:,:), fp(:,:,:,:),&
                            lo, hi, ng, dx, dt_grid)
         end select

         dt_proc = min(dt_grid, dt_proc)
      end do

      ! This sets dt to be the min of dt_proc over all processors.
      call parallel_reduce(dt ,dt_proc ,MPI_MIN)

      if (dt .eq. dt_start) then
          dt = min(dx(1),dx(2))
          if (dm .eq. 3) dt = min(dt,dx(3))
      end if

      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,max_dt_growth*dtold)

      if (parallel_IOProcessor() .and. verbose .ge. 1) then
        write(6,1000) lev,dt
      end if
1000  format("Computing dt at level ",i2," to be ... ",e15.8)

   end subroutine estdt

   subroutine estdt_2d (vel,s,gp,ext_vel_force,lo,hi,ng,dx,dt)

      integer, intent(in) :: lo(:), hi(:), ng

      real (kind = dp_t), intent(in   ) :: vel(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) :: gp(lo(1)-1:,lo(2)-1:,:)  
      real (kind = dp_t), intent(in   ) :: ext_vel_force(lo(1)-1:,lo(2)-1:,:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(inout) :: dt

!     Local variables
      real (kind = dp_t)  u,v,fx,fy
      real (kind = dp_t)  eps
      integer :: i, j

      eps = 1.0e-8

      u  = 0.0D0 
      v  = 0.0D0 
      fx = 0.0D0 
      fy = 0.0D0 

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          u  = max(u ,abs(vel(i,j,1)))
          v  = max(v ,abs(vel(i,j,2)))
          fx = max(fx,abs(gp(i,j,1)/s(i,j)-ext_vel_force(i,j,1)))
          fy = max(fy,abs(gp(i,j,2)/s(i,j)-ext_vel_force(i,j,2)))
        enddo
      enddo

      if (u .gt. eps) dt = min(dt,dx(1)/u)
      if (v .gt. eps) dt = min(dt,dx(2)/v)

      if (fx > eps) &
        dt = min(dt,sqrt(2.0D0 *dx(1)/fx))

      if (fy > eps) &
        dt = min(dt,sqrt(2.0D0 *dx(2)/fy))

   end subroutine estdt_2d

   subroutine estdt_3d (vel,s,gp,ext_vel_force,lo,hi,ng,dx,dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) :: vel(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(in   ) :: gp(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) :: ext_vel_force(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(inout) :: dt

!     Local variables
      real (kind = dp_t)  u,v,w,fx,fy,fz
      real (kind = dp_t)  eps
      integer :: i, j, k

      eps = 1.0e-8

      u  = 0.0D0 
      v  = 0.0D0 
      w  = 0.0D0
      fx = 0.0D0 
      fy = 0.0D0 
      fz = 0.0D0 

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          u  = max(u ,abs(vel(i,j,k,1)))
          v  = max(v ,abs(vel(i,j,k,2)))
          w  = max(w ,abs(vel(i,j,k,3)))
          fx = max(fx,abs(gp(i,j,k,1)/s(i,j,k)-ext_vel_force(i,j,k,1)))
          fy = max(fy,abs(gp(i,j,k,2)/s(i,j,k)-ext_vel_force(i,j,k,2)))
          fz = max(fz,abs(gp(i,j,k,3)/s(i,j,k)-ext_vel_force(i,j,k,3)))
        enddo
      enddo
      enddo

      if (u .gt. eps) dt = min(dt,dx(1)/u)
      if (v .gt. eps) dt = min(dt,dx(2)/v)
      if (w .gt. eps) dt = min(dt,dx(3)/w)

      if (fx > eps) &
        dt = min(dt,sqrt(2.0D0 *dx(1)/fx))

      if (fy > eps) &
        dt = min(dt,sqrt(2.0D0 *dx(2)/fy))

      if (fz > eps) &
        dt = min(dt,sqrt(2.0D0 *dx(3)/fz))

   end subroutine estdt_3d

end module estdt_module
