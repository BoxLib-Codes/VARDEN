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
  public :: initdata, initdata_on_level

contains

  subroutine initdata_on_level(u,s,dx,bc)

    use multifab_physbc_module
    use probin_module, only : nscal

    type(multifab) , intent(inout) :: u,s
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim)
    integer :: i,ng,dm

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

    call multifab_physbc(u,1,1,   dm,   bc)
    call multifab_physbc(s,1,dm+1,nscal,bc)

  end subroutine initdata_on_level

  subroutine initdata(nlevs,u,s,dx,bc,mla)

    use multifab_physbc_module
    use probin_module, only : nscal

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:),s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim)
    integer :: i,ng,dm,n

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
             call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx(n,:))
          case (3)
             call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx(n,:))
          end select
       end do

       call multifab_fill_boundary(u(n))
       call multifab_fill_boundary(s(n))

       call multifab_physbc(u(n),1,1,   dm,   bc(n))
       call multifab_physbc(s(n),1,dm+1,nscal,bc(n))

    enddo

    do n = nlevs,2,-1
       call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
    enddo

    do n = 2,nlevs
       call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,1,dm)
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      bc(n-1),bc(n),1,dm+1,nscal)
    enddo

  end subroutine initdata

  subroutine initdata_2d (u,s,lo,hi,ng,dx)

    use probin_module, only : prob_hi

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)

    !     Local variables
    integer :: i, j, comp, jhalf
    real (kind = dp_t) :: x,y,r,cpx,cpy,spx,spy,Pi
    real (kind = dp_t) :: velfact
    real (kind = dp_t) :: ro,r_pert
    real (kind = dp_t) :: r0,denfact

    Pi = 4.0_dp_t*atan(1.0) 
    velfact = 1.0_dp_t

    !     ro is the density of air
    ro = 1.2e-3

    r_pert = .025

    u = ZERO
    s = ZERO

    jhalf = (lo(2)+hi(2))/2

    if (.false.) then

       do j = lo(2), hi(2)
          !       y = (float(j)+HALF) * dx(2) / prob_hi(2)
          y = (float(j)+HALF) * dx(2) 
          do i = lo(1), hi(1)
             !          x = (float(i)+HALF) * dx(1) / prob_hi(1)
             x = (float(i)+HALF) * dx(1)

             !          Initial data for Poiseuille flow.
             !          u(i,j,2) = ONE * (x) * (ONE - x)

             !          Initial data for vortex-in-a-box
             !          if (x .le. 0.5) then
             !            spx = sin(Pi*x)
             !            cpx = cos(Pi*x)
             !          else 
             !            spx =  sin(Pi*(1.0-x))
             !            cpx = -cos(Pi*(1.0-x))
             !          end if
             !          if (y .le. 0.5) then
             !            spy = sin(Pi*y)
             !            cpy = cos(Pi*y)
             !          else 
             !            spy =  sin(Pi*(1.0-y))
             !            cpy = -cos(Pi*(1.0-y))
             !          end if

             !          spx = sin(Pi*x)
             !          spy = sin(Pi*y)
             !          cpx = cos(Pi*x)
             !          cpy = cos(Pi*y)

             !          u(i,j,1) =  TWO*velfact*spy*cpy*spx*spx
             !          u(i,j,2) = -TWO*velfact*spx*cpx*spy*spy

             u(i,j,1) = tanh(30.0_dp_T*(0.25_dp_t - abs(y-0.5_dp_t)))
             u(i,j,2) = 0.05d0 * sin(2.0_dp_t*Pi*x)

             !          u(i,j,1) = sin(y)
             !          u(i,j,2) = cos(x)

             s(i,j,1) = ONE
             r = sqrt((x-HALF)**2 + (y-HALF)**2)
             s(i,j,2) = merge(1.2_dp_t,ONE,r .lt. 0.15)

          enddo
       enddo

    else if (.false.) then

       u = ZERO
       s = ONE
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = (float(i)+HALF) * dx(1)
             y = (float(j)+HALF) * dx(2)
             if (y.lt.0.50) then
                u(i,j,1) = ONE
             else
                u(i,j,1) = -ONE
             end if
          enddo
       enddo

    else

       u = ZERO
       s = ONE
       r0 = 0.15d0
       denfact = 20.d0
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             y = (float(j)+HALF) * dx(2) / prob_hi(2)
             x = (float(i)+HALF) * dx(1) / prob_hi(1)
             r = sqrt((x-HALF)**2 + (y-HALF)**2)
             s(i,j,1) = ONE + HALF*(denfact-ONE)*(ONE-tanh(30.*(r-r0)))
             s(i,j,1) = ONE / s(i,j,1)

          end do
       end do

    end if

    !     Impose inflow conditions if grid touches inflow boundary.
    !     do i = lo(1), hi(1)
    !       x = (float(i)+HALF) * dx(1) / prob_hi(1)
    !       u(lo(1)       :hi(1)        ,lo(2)-1,2) = INLET_VY * FOUR*x*(ONE-x)
    !     end do

    !     Impose inflow conditions if grid touches inflow boundary.
    if (lo(2) .eq. 0) then
       u(lo(1)       :hi(1)        ,lo(2)-1,1) = INLET_VX
       u(lo(1)       :hi(1)        ,lo(2)-1,2) = INLET_VY
       s(lo(1)       :hi(1)        ,lo(2)-1,1) = INLET_DEN
       s(lo(1)       :hi(1)        ,lo(2)-1,2) = INLET_TRA
       s(lo(1)-1:hi(1)+1,lo(2)-1,2) = ONE
    end if

    if (size(s,dim=3).gt.2) then
       do comp = 3, size(s,dim=3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                !          s(i,j,comp) = ONE
                y = (float(j)+HALF) * dx(2) / prob_hi(2)
                x = (float(i)+HALF) * dx(1) / prob_hi(1)
                r = sqrt((x-HALF)**2 + (y-HALF)**2)
                s(i,j,comp) = r
             end do
          end do
       end do
    end if

  end subroutine initdata_2d

  subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_hi)

    use probin_module, only : prob_hi

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)

    !     Local variables
    integer :: i, j, k, comp
    integer :: imid,jmid,kmid
    real (kind = dp_t) :: Pi
    real (kind = dp_t) :: x,y,z
    real (kind = dp_t) :: cpx,cpy,cpz,spx,spy,spz
    real (kind = dp_t) :: velfact

    Pi = 4.0_dp_t*atan(1.0)
    velfact = ONE

    u = ZERO
    s = ZERO

    do k = lo(3), hi(3)
       z = (float(k)+HALF) * dx(3) / prob_hi(3)
       do j = lo(2), hi(2)
          y = (float(j)+HALF) * dx(2) / prob_hi(2)
          do i = lo(1)-ng, hi(1)+ng
             x = (float(i)+HALF) * dx(1) / prob_hi(1)

             !        spy = sin(Pi*y)
             !        spz = sin(Pi*z)
             !        cpy = cos(Pi*y)
             !        cpz = cos(Pi*z)
             !        u(i,j,k,2) =  TWO*velfact*spz*cpz*spy*spy
             !        u(i,j,k,3) = -TWO*velfact*spy*cpy*spz*spz
             !        u(i,j,k,1) = ZERO

             !        spz = sin(Pi*z)
             !        spx = sin(Pi*x)
             !        cpz = cos(Pi*z)
             !        cpx = cos(Pi*x)
             !        u(i,j,k,3) =  TWO*velfact*spx*cpx*spz*spz
             !        u(i,j,k,1) = -TWO*velfact*spz*cpz*spx*spx
             !        u(i,j,k,2) = ZERO

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
       do comp = 2, size(s,dim=4)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   s(i,j,k,comp) = ONE
                end do
             end do
          end do
       end do
    end if

    u = ZERO
    s = ONE
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             x = (float(i)+HALF) * dx(1)
             if (x.lt.0.50) u(i,j,k,1) = ONE
          enddo
       enddo
    enddo

    !     imid = (lo(1)+hi(1))/2+1
    !     jmid = (lo(2)+hi(2))/2+1
    !     kmid = (lo(3)+hi(3))/2+1
    !     u(imid  ,jmid  ,kmid  ,:) =  ONE
    !     u(imid-1,jmid-1,kmid-1,:) = -ONE

  end subroutine initdata_3d

end module init_module
