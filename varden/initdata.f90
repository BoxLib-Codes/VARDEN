module init_module

  use bl_types
  use bc_module
  use mg_bc_module
  use multifab_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t

contains

   subroutine initdata (u,s,dx,prob_hi)

      type(multifab) , intent( in) :: u,s
      real(kind=dp_t), intent( in) :: dx(:)
      real(kind=dp_t), intent( in) :: prob_hi(:)

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      integer :: i

      ng = u%ng
      dm = u%dim

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx, prob_hi)
            case (3)
              call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx, prob_hi)
         end select
      end do

   end subroutine initdata

   subroutine initdata_2d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(2), hi(2), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(2)
      real (kind = dp_t), intent(in ) :: prob_hi(2)

!     Local variables
      integer :: i, j, n
      real (kind = dp_t) :: x,y,r,cpx,cpy,spx,spy,Pi
      real (kind = dp_t) :: velfact

      Pi = 4.0_dp_t*atan(1.0) 
      velfact = 1.0_dp_t

      u = ZERO
      s = ZERO

      do j = lo(2), hi(2)
        y = (float(j)+0.5) * dx(2) / prob_hi(2)
        do i = lo(1), hi(1)
           x = (float(i)+0.5) * dx(1) / prob_hi(1)
           spx = sin(Pi*x)
           spy = sin(Pi*y)
           cpx = cos(Pi*x)
           cpy = cos(Pi*y)
           u(i,j,1) =  2.0_dp_t*velfact*spy*cpy*spx*spx
           u(i,j,2) = -2.0_dp_t*velfact*spx*cpx*spy*spy
           s(i,j,1) = 1.0_dp_t
!          r = sqrt((x-0.5)**2 + (y-0.5)**2)
!          s(i,j,1) = merge(1.2_dp_t,1.0_dp_t,r .lt. 0.15)
        enddo
      enddo

      if (size(s,dim=3).gt.1) then
        do n = 2, size(s,dim=3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
!          s(i,j,n) = 1.0_dp_t
           y = (float(j)+0.5) * dx(2) / prob_hi(2)
           x = (float(i)+0.5) * dx(1) / prob_hi(1)
           r = sqrt((x-0.5)**2 + (y-0.5)**2)
           s(i,j,n) = r
        end do
        end do
        end do
      end if

   end subroutine initdata_2d

   subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(3), hi(3), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(3)
      real (kind = dp_t), intent(in ) :: prob_hi(3)
    
!     Local variables
      integer :: i, j, k, n
      real (kind = dp_t) :: x,y,cpx,cpy,spx,spy,Pi
      real (kind = dp_t) :: velfact

      Pi = 4.0_dp_t*atan(1.0)
      velfact = 1.0_dp_t

      u = ZERO
      s = ZERO

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        y = (float(j)+0.5) * dx(2) / prob_hi(2)
        do i = lo(1)-ng, hi(1)+ng
         x = (float(i)+0.5) * dx(1) / prob_hi(1)
         spx = sin(Pi*x)
         spy = sin(Pi*y)
         cpx = cos(Pi*x)
         cpy = cos(Pi*y)
         u(i,j,k,1) =  2.0_dp_t*velfact*spy*cpy*spx*spx
         u(i,j,k,2) = -2.0_dp_t*velfact*spx*cpx*spy*spy
         u(i,j,k,3) = 0.0_dp_t
         s(i,j,k,1) = 1.0_dp_t
      enddo
      enddo
      enddo

      if (size(s,dim=4).gt.1) then
        do n = 1, size(s,dim=4)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           s(i,j,k,n) = 1.0_dp_t
        end do
        end do
        end do
        end do
      end if

   end subroutine initdata_3d

   subroutine define_bcs(phys_bc,norm_vel_bc,tang_vel_bc,scal_bc,press_bc,visc_coef)

   integer        , intent(in   ) :: phys_bc(:,:)
   integer        , intent(  out) :: norm_vel_bc(:,:)
   integer        , intent(  out) :: tang_vel_bc(:,:)
   integer        , intent(  out) :: scal_bc(:,:)
   integer        , intent(  out) :: press_bc(:,:)
   real(kind=dp_t), intent(in   ) :: visc_coef

   integer :: dm,i,n

   dm = size(phys_bc,dim=1)
 
   do n = 1, dm
   do i = 1, 2
      if (phys_bc(n,i) == WALL) then
          norm_vel_bc(n,i) = BC_DIR
          if (visc_coef > 0.0_dp_t) then
             tang_vel_bc(n,i) = BC_DIR
          else
             tang_vel_bc(n,i) = BC_NEU
          end if
              scal_bc(n,i) = BC_NEU
             press_bc(n,i) = BC_NEU
      else if (phys_bc(n,i) == INLET) then
          norm_vel_bc(n,i) = BC_DIR
          tang_vel_bc(n,i) = BC_DIR
              scal_bc(n,i) = BC_DIR
             press_bc(n,i) = BC_NEU
      else if (phys_bc(n,i) == OUTLET) then
          norm_vel_bc(n,i) = BC_NEU
          tang_vel_bc(n,i) = BC_NEU
              scal_bc(n,i) = BC_NEU
             press_bc(n,i) = BC_DIR
      end if
   end do
   end do

   end subroutine define_bcs

end module init_module
