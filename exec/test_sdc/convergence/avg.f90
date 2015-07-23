module average

  use bl_types
  use multifab_module

  implicit none

  real(dp_t), private, parameter :: ZERO = 0.0_dp_t
  real(dp_t), private, parameter :: ONE  = 1.0_dp_t


contains

  subroutine ml_cc_restriction_c(crse, cc, fine, cf, ir, nc)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer, intent(in)           :: cc, cf, ir(:)
    integer, intent(in), optional :: nc

    integer             :: i, n, lnc, lo(fine%dim), hi(fine%dim), lof(fine%dim)
    real(dp_t), pointer :: fp(:,:,:,:), cp(:,:,:,:)
    type(layout)        :: lacfine
    type(multifab)      :: cfine

    lnc = 1; if ( present(nc) ) lnc = nc

    call layout_build_coarse(lacfine, fine%la, ir)

    call build(cfine, lacfine, nc = lnc, ng = 0)

    do i = 1, fine%nboxes
       if ( remote(fine, i) ) cycle
       lof = lwb(get_pbox(fine, i))
       lo  = lwb(get_ibox(cfine,i))
       hi  = upb(get_ibox(cfine,i))
       do n = 1, lnc
          fp => dataptr(fine,  i, n+cf-1, 1)
          cp => dataptr(cfine, i, n,      1)
          select case (cfine%dim)
          case (1)
             call cc_restriction_1d(cp(:,1,1,1), lo, fp(:,1,1,1), lof, lo, hi, ir)
          case (2)
             call cc_restriction_2d(cp(:,:,1,1), lo, fp(:,:,1,1), lof, lo, hi, ir)
          case (3)
             call cc_restriction_3d(cp(:,:,:,1), lo, fp(:,:,:,1), lof, lo, hi, ir)
          end select
       end do
    end do

    call copy(crse, cc, cfine, 1, lnc)

    call destroy(cfine)

    call multifab_fill_boundary_c(crse,cc,nc)

  end subroutine ml_cc_restriction_c

  subroutine ml_cc_restriction(crse, fine, ir)
    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    integer,        intent(in)    :: ir(:)
    if ( crse%nc .ne. fine%nc ) then
       call bl_error('ml_cc_restriction: crse & fine must have same # of components')
    end if
    call ml_cc_restriction_c(crse, 1, fine, 1, ir, crse%nc)
  end subroutine ml_cc_restriction

  subroutine cc_restriction_1d(cc, loc, ff, lof, lo, hi, ir)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:), hi(:)
    real (dp_t), intent(inout) :: cc(loc(1):)
    real (dp_t), intent(in)    :: ff(lof(1):)
    integer,     intent(in)    :: ir(:)

    real (dp_t) :: fac
    integer     :: i, l

    fac = one/real(product(ir),kind=dp_t)

    do i = lo(1), hi(1)
       cc(i) = zero
       do l = 0, ir(1)-1
          cc(i) = cc(i) + ff(ir(1)*i+l)
       end do
       cc(i) = cc(i)*fac
    end do

  end subroutine cc_restriction_1d

  subroutine cc_restriction_2d(cc, loc, ff, lof, lo, hi, ir)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:), hi(:)
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):)
    real (dp_t), intent(in)    :: ff(lof(1):,lof(2):)
    integer,     intent(in)    :: ir(:)

    real (dp_t) :: fac
    integer     :: i, j, l, m

    fac = one/real(product(ir),kind=dp_t)

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          cc(i,j) = zero
          do m = 0, ir(2)-1
             do l = 0, ir(1)-1
                cc(i,j) = cc(i,j) + ff(ir(1)*i+l,ir(2)*j+m)
             end do
          end do
          cc(i,j) = cc(i,j)*fac
       end do
    end do

  end subroutine cc_restriction_2d

  subroutine cc_restriction_3d(cc, loc, ff, lof, lo, hi, ir)
    integer,     intent(in)    :: loc(:)
    integer,     intent(in)    :: lof(:)
    integer,     intent(in)    :: lo(:),hi(:)
    real (dp_t), intent(inout) :: cc(loc(1):,loc(2):,loc(3):)
    real (dp_t), intent(in)    :: ff(lof(1):,lof(2):,lof(3):)
    integer,     intent(in)    :: ir(:)

    real (dp_t) :: fac
    integer     :: i, j, k, l, m, n

    fac = one/real(product(ir),kind=dp_t)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             cc(i,j,k) = zero
             do n = 0, ir(3)-1
                do m = 0, ir(2)-1
                   do l = 0, ir(1)-1
                      cc(i,j,k) = cc(i,j,k) + ff(ir(1)*i+l,ir(2)*j+m,ir(3)*k+n)
                   end do
                end do
             end do
             cc(i,j,k) = cc(i,j,k)*fac
          end do
       end do
    end do

  end subroutine cc_restriction_3d

end module average
