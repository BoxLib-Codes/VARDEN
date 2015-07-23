module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
  use define_bc_module
  use multifab_module
  use multifab_fill_ghost_module
  use ml_cc_restriction_module
  use ml_layout_module

  implicit none

  private
  public :: initdata, impose_pressure_bcs, initdata_on_level

contains

  subroutine initdata_on_level(u,s,dx,bc)

    use multifab_physbc_module
    use probin_module, only : nscal

    type(multifab) , intent(inout) :: u,s
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim)
    integer :: i,ng,dm,n

    ng = u%ng
    dm = u%dim

    do i = 1, nfabs(u)
       uop => dataptr(u,i)
       sop => dataptr(s,i)
       lo =  lwb(get_box(u,i))
       hi =  upb(get_box(u,i))
       select case (dm)
       case (2)
          call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx)
       case (3)
          call bl_error('Three-D not coded yet')
       end select
    end do

    call multifab_fill_boundary(u)
    call multifab_fill_boundary(s)

    call multifab_physbc(u,1,1,   dm,   bc)
    call multifab_physbc(s,1,dm+1,nscal,bc)

  end subroutine initdata_on_level

  subroutine initdata(nlevs,u,s,dx,bc,mla)

    use multifab_physbc_module
    use probin_module, only: nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla


    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    integer :: ng,dm,i,n

    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx(n,:))
          case (3)
             call bl_error('Three-D not coded yet')
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

  subroutine initdata_2d(u,s,lo,hi,ng,dx)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: x, y, r, rsq, cos_theta, sin_theta
    real (kind = dp_t) :: x_c, y_c, rho_c, u_c, v_c, bigrsq

    rho_c = HALF
      x_c = HALF
      y_c = HALF

    ! u_c = ONE
    ! v_c = ONE

      u_c = ZERO
      v_c = ZERO

    bigrsq = 0.4**2

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          x   = (dble(i)+HALF)*dx(1)
          y   = (dble(j)+HALF)*dx(2)
          rsq = ( (x-x_c)**2 + (y-y_c)**2 ) / bigrsq
          r   = dsqrt(rsq)

          sin_theta = (y-y_c) / r
          cos_theta = (x-x_c) / r

          if (rsq < 1) then
              u(i,j,1) =   u_c - 1024.d0 * sin_theta * (ONE - r)**6 * rsq**3
              u(i,j,2) =   v_c + 1024.d0 * cos_theta * (ONE - r)**6 * rsq**3
              s(i,j,1) = rho_c + 0.5d0 * (ONE - rsq)**6
          else
              u(i,j,1) = u_c 
              u(i,j,2) = v_c
              s(i,j,1) = rho_c 
          end if

          s(i,j,2) = ZERO

       enddo
    enddo


  end subroutine initdata_2d

  subroutine impose_pressure_bcs(p,mla,mult)

    type(multifab ), intent(inout) :: p(:)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: mult

    type(box)           :: bx,pd
    integer             :: i,n,nlevs

    nlevs = size(p,dim=1)

    do n = 1,nlevs
       pd = layout_get_pd(mla%la(n))
       do i = 1, nfabs(p(n))
          bx = get_ibox(p(n),i)
          if (bx%lo(2) == pd%lo(2)) then
             bx%hi(2) = bx%lo(2)
             call setval(p(n),mult,bx)
          end if
       end do
    end do

  end subroutine impose_pressure_bcs

end module init_module
