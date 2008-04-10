module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
  use setbc_module
  use define_bc_module
  use multifab_module
  use multifab_fill_ghost_module
  use ml_restriction_module
  use ml_layout_module

  implicit none

  private
  public :: initdata, impose_pressure_bcs

contains

  subroutine initdata(nlevs,u,s,dx,prob_hi,bc,nscal,mla)

    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc(:)
    integer        , intent(in   ) :: nscal
    type(ml_layout), intent(inout) :: mla


    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    integer :: ng,dm,i,n

    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx(n,:), prob_hi)
          case (3)
             call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx(n,:), prob_hi)
          end select
       end do

       call multifab_fill_boundary(u(n))
       call multifab_fill_boundary(s(n))

       call multifab_physbc(u(n),1,1,   dm,   bc(n))
       call multifab_physbc(s(n),1,dm+1,nscal,bc(n))

    enddo

    do n=nlevs,2,-1
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,1,dm)
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,dm+1,nscal)
    enddo

  end subroutine initdata

  subroutine initdata_2d(u,s,lo,hi,ng,dx,prob_hi)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t)  :: x,y,dist
    real (kind = dp_t)  :: xblob = 0.5d0, yblob = 0.5d0, densfact = 2.0d0
    real (kind = dp_t)  :: blobrad = 0.1d0
   
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          u(i,j,1) = ZERO
          u(i,j,2) = ZERO
          s(i,j,1) = ONE
          s(i,j,2) = ZERO
       enddo
    enddo

    ! add a density perturbation
!    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,1) = TWO
!    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,1) = TWO
!    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,1) = TWO
!    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,1) = TWO
    do j=lo(2),hi(2)
!       y = lo(2) + dx(2)*((j-lo(2)) + HALF)
       y = dx(2)*(j + HALF)
       do i=lo(1),hi(1)
!          x = lo(1) + dx(1)*((i-lo(1)) + HALF)
          x = dx(1)*((i) + HALF)
          dist = SQRT((x-xblob)**2 + (y-yblob)**2)
          s(i,j,1) = ONE + HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
          s(i,j,2) = s(i,j,1)
!          write(*,*) i,j,s(i,j,1), s(i,j,2)
       enddo
    enddo

    ! add a tracer perturbation
!    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,2) = ONE
!    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,2) = ONE
!    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,2) = ONE
!    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,2) = ONE

  end subroutine initdata_2d

  subroutine initdata_3d(u,s,lo,hi,ng,dx,prob_hi)

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
             u(i,j,k,2) = ONE
             u(i,j,k,3) = ZERO
             s(i,j,k,1) = ONE
             s(i,j,k,2) = ZERO
          enddo
       enddo
    enddo

    ! add a density perturbation
    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,(hi(3)+1)/2-1,1) = TWO
    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,(hi(3)+1)/2-1,1) = TWO
    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,(hi(3)+1)/2-1,1) = TWO
    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,(hi(3)+1)/2-1,1) = TWO
    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,(hi(3)+1)/2  ,1) = TWO
    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,(hi(3)+1)/2  ,1) = TWO
    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,(hi(3)+1)/2  ,1) = TWO
    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,(hi(3)+1)/2  ,1) = TWO

    ! add a tracer perturbation
    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,(hi(3)+1)/2-1,2) = ONE
    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,(hi(3)+1)/2-1,2) = ONE
    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,(hi(3)+1)/2-1,2) = ONE
    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,(hi(3)+1)/2-1,2) = ONE
    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,(hi(3)+1)/2  ,2) = ONE
    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,(hi(3)+1)/2  ,2) = ONE
    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,(hi(3)+1)/2  ,2) = ONE
    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,(hi(3)+1)/2  ,2) = ONE
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
