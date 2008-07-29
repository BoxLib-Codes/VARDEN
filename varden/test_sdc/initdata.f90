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
  use probin_module,  only : nscal, prob_hi_y

  implicit none

  private
  public :: initdata, initdata_on_level

contains

  subroutine initdata(nlevs,u,s,dx,bc,mla)


    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    integer :: i,n,ng,dm


    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n), i) ) cycle
          uop => dataptr(u(n), i)
          sop => dataptr(s(n), i)
          lo =  lwb(get_box(u(n), i))
          hi =  upb(get_box(u(n), i))
          select case (dm)
          case (2)
             call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx(n,:))
          case (3)
             call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx(n,:))
          end select
       end do

       call multifab_fill_boundary(u(n))
       call multifab_fill_boundary(s(n))

       call multifab_physbc(u(n),1,1,   dm,   bc(n),dx(n,:),t=zero)
       call multifab_physbc(s(n),1,dm+1,nscal,bc(n),dx(n,:),t=zero)

    enddo

    do n=2,nlevs
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
            bc(n-1),bc(n),1,1,dm,dx(n-1:n,:),t=zero)
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
            bc(n-1),bc(n),1,dm+1,nscal,dx(n-1:n,:),t=zero)
    enddo

  end subroutine initdata

  subroutine initdata_on_level(u,s,dx,bc,la)
    
    use multifab_physbc_module
    use probin_module, only : nscal
    
    type(multifab) , intent(inout) :: u,s
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc
    type(layout)   , intent(inout) :: la
    
    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim)
    integer :: i,ng,dm,n
    
    ng = u%ng
    dm = u%dim
    
    do i = 1, u%nboxes
       if ( multifab_remote(u,i) ) cycle
       uop => dataptr(u,i)
       sop => dataptr(s,i)
       lo =  lwb(get_box(u,i))
       hi =  upb(get_box(u,i))
       select case (dm)
       case (2)
          call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx)
       case (3)
          call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx)
       end select
    end do

    call multifab_fill_boundary(u)
    call multifab_fill_boundary(s)
    
    call multifab_physbc(u,1,1,   dm,   bc,dx,t=zero)
    call multifab_physbc(s,1,dm+1,nscal,bc,dx,t=zero)
    
  end subroutine initdata_on_level

  subroutine initdata_2d (u,s,lo,hi,ng,dx)

    use probin_module, only: boussinesq

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t), parameter  :: lambda = fourth
    real (kind = dp_t), parameter  :: Um = one
    real (kind = dp_t), parameter  :: delta = 2.d-2
!    real (kind = dp_t)  :: x,y,dist
!    real (kind = dp_t)  :: xblob1 = 0.5d0, yblob1 = 0.5d0, densfact = 2.0d0
!    real (kind = dp_t)  :: xblob2 = 0.2d0, yblob2 = 0.2d0
!    real (kind = dp_t)  :: xblob3 = 0.2d0, yblob3 = 0.8d0
!    real (kind = dp_t)  :: blobrad = 0.1d0
   

    ! shear layer w/ roll-up
!    do i=lo(1),hi(1)   
!       do j=lo(2),hi(2)
!          u(i,j,1) = Um*(one + lambda * tanh(two*(dx(2)*j-half)/delta))
!          u(i,j,2) = zero
!          s(i,j,1) = ONE     ! density          
!          s(i,j,4) = zero    ! species C
!          if (j*dx(2) .le. half*prob_hi_y) then
!             s(i,j,2) = one    ! species A
!             s(i,j,3) = zero    ! species B
!          else
!             s(i,j,2) = zero    ! species A
!             s(i,j,3) = ten !one    ! species B
!       end if
!       enddo
!    enddo

    ! simple shear layer/plume
    do i=lo(1),hi(1)   
       do j=lo(2),hi(2)
          if (j*dx(2) < half*prob_hi_y) then
             u(i,j,1) = zero !one
             u(i,j,2) = zero
             s(i,j,1) = ONE     ! density
             s(i,j,2) = half*(one + sin(two*M_PI*i*dx(1)))
!one    ! species A
             s(i,j,3) = half*(one + sin(two*M_PI*i*dx(1))) !zero    ! species B
             s(i,j,4) = zero    ! species C
          else
             u(i,j,1) = zero !two !ten
             u(i,j,2) = zero
             s(i,j,1) = ONE     ! density
             s(i,j,2) = half*(one + sin(two*M_PI*i*dx(1))) !zero    ! species A
             s(i,j,3) = half*(one + sin(two*M_PI*i*dx(1)))
!one    ! species B
          s(i,j,4) = zero    ! species C
       end if
       enddo
    enddo

    


    ! add a density perturbation
!    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,1) = TWO
!    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,1) = TWO
!    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,1) = TWO
!    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,1) = TWO

    ! add a tracer perturbation
!    s((hi(1)+1)/2-1,(hi(2)+1)/2-1,2) = ONE
!    s((hi(1)+1)/2  ,(hi(2)+1)/2-1,2) = ONE
!    s((hi(1)+1)/2-1,(hi(2)+1)/2  ,2) = ONE
!    s((hi(1)+1)/2  ,(hi(2)+1)/2  ,2) = ONE

    ! bubble
!    do j=lo(2),hi(2)
       !y = lo(2) + dx(2)*((j-lo(2)) + HALF)
!       y = dx(2)*(j + HALF)
!       do i=lo(1),hi(1)
         !x = lo(1) + dx(1)*((i-lo(1)) + HALF)
!          x = dx(1)*((i) + HALF)
!          dist = SQRT((x-xblob1)**2 + (y-yblob1)**2)
          !s(i,j,1) = ONE + HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
!          s(i,j,2) = HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
!          dist = SQRT((x-xblob2)**2 + (y-yblob2)**2)
!          s(i,j,3) = HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
!          dist = SQRT((x-xblob3)**2 + (y-yblob3)**2)
!          s(i,j,4) = HALF*(densfact-ONE)*(ONE-TANH(30.*(dist-blobrad))) 
!       enddo
!    enddo
  
  end subroutine initdata_2d

  subroutine initdata_3d (u,s,lo,hi,ng,dx)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j, k
    real(kind=dp_t) :: xloc,yloc,zloc,dist

    ! zero initial velocity
    ! density = 1
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             u(i,j,k,1) = ZERO
             u(i,j,k,2) = ZERO
             s(i,j,k,1) = ONE
             s(i,j,k,2) = ZERO

          enddo
       enddo
    enddo

    ! add two "bubbles" of higher density
    ! one centered over fine grid
    ! one centered over coarse grid
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             xloc = (i+HALF)*dx(1)
             yloc = (j+HALF)*dx(2)
             zloc = (k+HALF)*dx(3)

             ! use this for one bubble problem
             dist = sqrt((xloc-0.5d0)**2 + (yloc-0.5d0)**2 + (zloc-0.5d0)**2)
             s(i,j,k,1) = s(i,j,k,1) + (1.0d0 - tanh(dist/0.05d0))

             ! use this for two bubble problem            
             !             dist = sqrt((xloc-0.75d0)**2 + (yloc-0.5d0)**2 + (zloc-0.5d0)**2)
             !             s(i,j,k,1) = s(i,j,k,1) + (1.0d0 - tanh(dist/0.05d0))
             !             dist = sqrt((xloc-0.25d0)**2 + (yloc-0.5d0)**2 + (zloc-0.5d0)**2)
             !             s(i,j,k,1) = s(i,j,k,1) + (1.0d0 - tanh(dist/0.05d0))

          enddo
       enddo
    enddo

  end subroutine initdata_3d

end module init_module
