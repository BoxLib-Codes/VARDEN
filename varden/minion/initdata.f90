module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
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

    do n=2,nlevs
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
            bc(n-1),bc(n),1,1,dm)
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
            bc(n-1),bc(n),1,dm+1,nscal)
    enddo

  end subroutine initdata

  subroutine initdata_2d (u,s,lo,hi,ng,dx,prob_hi)

    use probin_module, only: boussinesq

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)

    !     Local variables
    integer :: i, j
    real(kind=dp_t) :: xloc,yloc,dist
    real(kind=dp_t) :: sig
    real(kind=dp_t) :: rho0,rho1,xBubble,yBubble,rBubble,rhoBubble,tracerConc,yJump

    ! initial velocity = 0
    ! initial  density = 1
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          u(i,j,1) = ZERO
          u(i,j,2) = ZERO
          s(i,j,1) = ONE
          s(i,j,2) = ZERO

       enddo
    enddo

    ! minion initial data
    rho0    = 1.0100d0
    rho1    = 1.0000d0

    xBubble = 0.5d0
    yBubble = 0.25d0
    rBubble = 0.05d0
    yJump   = 0.5d0

    rhoBubble = 1.005d0
    tracerConc = 1.00d0

    sig = 0.02d0

    if (boussinesq .eq. 1) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          xloc = (i+HALF)*dx(1)
          yloc = (j+HALF)*dx(2)

          ! use this for one bubble problem
          dist = sqrt((xloc-xBubble)**2 + (yloc-yBubble)**2)
          if (dist .lt. rBubble) then
             s(i,j,2) = tracerConc
          else
             s(i,j,2) = s(i,j,2) + (tracerConc)*exp(-(dist-rBubble)/sig)
             s(i,j,2) = max(s(i,j,2),0.d0)
          endif
       enddo
       enddo
    else
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          xloc = (i+HALF)*dx(1)
          yloc = (j+HALF)*dx(2)

          ! use this for one bubble problem
          if(yloc.le. 0.5)then
             s(i,j,1) = rho0
          else
             s(i,j,1) = rho0+(rho1-rho0)*(yloc-yJump)/(1.0d0-yJump)
          endif
          dist = sqrt((xloc-xBubble)**2 + (yloc-yBubble)**2)
          if (dist .lt. rBubble) then
             s(i,j,1) = rhoBubble
             s(i,j,2) = tracerConc
          else
             s(i,j,1) = s(i,j,1) + (rhoBubble-rho0)*exp(-(dist-rBubble)/sig)
             s(i,j,2) = s(i,j,2) +     (tracerConc)*exp(-(dist-rBubble)/sig)
             s(i,j,2) = max(s(i,j,2),0.d0)
          endif
       enddo
       enddo
    end if

  end subroutine initdata_2d

  subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_hi)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)

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
