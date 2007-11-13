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
            case (3)
              call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx, prob_hi)
         end select
      end do

      call multifab_fill_boundary(u)
      call multifab_fill_boundary(s)

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         lo =  lwb(get_box(u, i))
         select case (dm)
            case (2)
              do n = 1,dm
                call setbc_2d(uop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
              end do
              do n = 1,nscal
                call setbc_2d(sop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
            case (3)
              do n = 1,dm
                call setbc_3d(uop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
              end do
              do n = 1,nscal
                call setbc_3d(sop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
              end do
         end select
      end do

   end subroutine initdata

   subroutine initdata_2d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: prob_hi(:)

!     Local variables
      integer :: i, j
      real(kind=dp_t) :: xloc,yloc,dist

      ! zero initial velocity
      ! density = 1
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            u(i,j,1) = ZERO
            u(i,j,2) = ZERO
            s(i,j,1) = ONE
            s(i,j,2) = ZERO

         enddo
      enddo

      ! add two "bubbles" of higher density
      ! one centered over fine grid
      ! one centered over coarse grid
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)

            xloc = (i+HALF)*dx(1)
            yloc = (j+HALF)*dx(2)

            dist = sqrt((xloc-0.75d0)**2 + (yloc-0.5d0)**2)
            s(i,j,1) = s(i,j,1) + (1.0d0 - tanh(dist/0.05d0))

            dist = sqrt((xloc-0.25d0)**2 + (yloc-0.5d0)**2)
            s(i,j,1) = s(i,j,1) + (1.0d0 - tanh(dist/0.05d0))

         enddo
      enddo

   end subroutine initdata_2d

   subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(in ) :: prob_hi(:)
    
!     Local variables
      integer :: i, j, k

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               u(i,j,k,1) = ZERO
               u(i,j,k,2) = ZERO
               u(i,j,k,3) = ZERO
               s(i,j,k,1) = ONE
               s(i,j,k,2) = ZERO

            enddo
         enddo
      enddo

   end subroutine initdata_3d

   subroutine impose_pressure_bcs(p,mla,mult)

     type(multifab ), intent(inout) :: p(:)
     type(ml_layout), intent(in   ) :: mla
     real(kind=dp_t), intent(in   ) :: mult
 
     type(box)           :: bx,pd
     integer             :: i,n,nlevs
     
     nlevs = size(p,dim=1)

     do n = 1,nlevs
        pd = layout_get_pd(mla%la(n))
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
