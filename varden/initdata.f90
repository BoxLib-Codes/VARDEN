module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
  use setbc_module
  use define_bc_module
  use multifab_module

  implicit none

contains

   subroutine initdata (u,s,dx,prob_hi,bc,nscal)

      type(multifab) , intent(inout) :: u,s
      real(kind=dp_t), intent(in   ) :: dx(:)
      real(kind=dp_t), intent(in   ) :: prob_hi(:)
      type(bc_level) , intent(in   ) :: bc
      integer        , intent(in   ) :: nscal

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      integer :: i,n
      logical :: is_vel

      ng = u%ng
      dm = u%dim
 
      is_vel = .true.
      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx, prob_hi)
              do n = 1,dm
                call setbc_2d(uop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
              end do
              do n = 1,nscal
                call setbc_2d(sop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
            case (3)
              call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx, prob_hi)
              do n = 1,dm
                call setbc_3d(uop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
              end do
              do n = 1,nscal
                call setbc_3d(sop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
         end select
      end do
      call multifab_fill_boundary(u)
      call multifab_fill_boundary(s)

   end subroutine initdata

   subroutine initdata_2d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: prob_hi(:)

!     Local variables
      integer :: i, j, n
      real (kind = dp_t) :: x,y,r,cpx,cpy,spx,spy,Pi
      real (kind = dp_t) :: velfact
      real (kind = dp_t) :: ro

      Pi = 4.0_dp_t*atan(1.0) 
      velfact = 1.0_dp_t

!     ro is the density of air
      ro = 1.2e-3

      u = ZERO
      s = ZERO

      do j = lo(2), hi(2)
        y = (float(j)+HALF) * dx(2) / prob_hi(2)
        do i = lo(1), hi(1)
           x = (float(i)+HALF) * dx(1) / prob_hi(1)

!          Initial data for Poiseuille flow.
!          u(i,j,2) = ONE * (x) * (ONE - x)

!          Initial data for vortex-in-a-box
!          spx = sin(Pi*x)
!          spy = sin(Pi*y)
!          cpx = cos(Pi*x)
!          cpy = cos(Pi*y)
!          u(i,j,1) =  TWO*velfact*spy*cpy*spx*spx
!          u(i,j,2) = -TWO*velfact*spx*cpx*spy*spy

           s(i,j,1) = ONE
           r = sqrt((x-HALF)**2 + (y-HALF)**2)
           s(i,j,2) = merge(1.2_dp_t,ONE,r .lt. 0.15)
        enddo
      enddo

!     Impose inflow conditions if grid touches inflow boundary.
!     if (lo(2) .eq. 0) then
!        u(lo(1)       :hi(1)        ,lo(2)-1,1) = INLET_VX
!        u(lo(1)       :hi(1)        ,lo(2)-1,2) = INLET_VY
!        s(lo(1)       :hi(1)        ,lo(2)-1,1) = INLET_DEN
!        s(lo(1)       :hi(1)        ,lo(2)-1,2) = INLET_TRA
!        s(lo(1)-1:hi(1)+1,lo(2)-1,2) = ONE
!     end if

      if (size(s,dim=3).gt.2) then
        do n = 3, size(s,dim=3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
!          s(i,j,n) = ONE
           y = (float(j)+HALF) * dx(2) / prob_hi(2)
           x = (float(i)+HALF) * dx(1) / prob_hi(1)
           r = sqrt((x-HALF)**2 + (y-HALF)**2)
           s(i,j,n) = r
        end do
        end do
        end do
      end if

   end subroutine initdata_2d

   subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: prob_hi(:)
    
!     Local variables
      integer :: i, j, k, n
      real (kind = dp_t) :: x,y,cpx,cpy,spx,spy,Pi
      real (kind = dp_t) :: velfact

      Pi = 4.0_dp_t*atan(1.0)
      velfact = ONE

      u = ZERO
      s = ZERO

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        y = (float(j)+HALF) * dx(2) / prob_hi(2)
        do i = lo(1)-ng, hi(1)+ng
         x = (float(i)+HALF) * dx(1) / prob_hi(1)
         spx = sin(Pi*x)
         spy = sin(Pi*y)
         cpx = cos(Pi*x)
         cpy = cos(Pi*y)
         u(i,j,k,1) =  TWO*velfact*spy*cpy*spx*spx
         u(i,j,k,2) = -TWO*velfact*spx*cpx*spy*spy
         u(i,j,k,3) = ZERO
         s(i,j,k,1) = ONE
      enddo
      enddo
      enddo

      if (size(s,dim=4).gt.1) then
        do n = 2, size(s,dim=4)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           s(i,j,k,n) = ONE
        end do
        end do
        end do
        end do
      end if

   end subroutine initdata_3d

   subroutine impose_pressure_bcs(p,la_tower,mult)

     type(multifab ), intent(inout) :: p(:)
     type(layout   ), intent(in   ) :: la_tower(:)
     real(kind=dp_t), intent(in   ) :: mult
 
     type(box) :: bx,pd
     integer :: i,n,nlevs
     
     nlevs = size(p,dim=1)

     do n = 1,nlevs
        pd = layout_get_pd(la_tower(n))
        do i = 1, p(n)%nboxes; if ( remote(p(n),i) ) cycle
           bx = get_ibox(p(n),i)
           if (bx%lo(2) == pd%lo(2)) then
             bx%hi(2) = bx%lo(2)
             call setval(p(n),mult,bx)
           end if
        end do
     end do

   end subroutine impose_pressure_bcs

end module init_module
