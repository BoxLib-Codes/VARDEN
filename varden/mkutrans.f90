module mkutrans_module

  use bl_types
  use multifab_module
  use slope_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

contains

      subroutine mkutrans_2d(vel,utrans,force,lo,dx,dt,ng,bc,irz)

      integer, intent(in) :: lo(2),ng

      real(kind=dp_t), intent(in) ::     vel(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in) ::   force(lo(1)- 1:,lo(2)- 1:,:)

      real(kind=dp_t), intent(inout) ::  utrans(lo(1)-2:,lo(2)-2:,:)

      real(kind=dp_t),intent(in) :: dt,dx(2)
      integer        ,intent(in) :: bc(2,2)
      integer        ,intent(in) :: irz

      real(kind=dp_t), allocatable::  velx(:,:,:)
      real(kind=dp_t), allocatable::  vely(:,:,:)

!     Local variables
      real(kind=dp_t) hx, hy, dth
      real(kind=dp_t) ulft,urgt,vbot,vtop

      real(kind=dp_t) :: eps

      integer :: hi(2)
      integer :: i,j,is,js,ie,je
      integer :: slope_order = 4

      logical :: test, do_refl

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng+1)

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      eps = 1.0e-8              ! FIXME what should EPS really be?

      dth = HALF * dt

      hx = dx(1)
      hy = dx(2)

      allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

      do_refl = (irz.eq.1)
      call slopex_2d(vel,velx,lo,ng,2,bc,do_refl,slope_order)
      call slopey_2d(vel,vely,lo,ng,2,bc,        slope_order)

!     Create the x-velocity to be used for transverse derivatives.
      do j = js-1,je+1 
        do i = is,ie+1 

          urgt = vel(i  ,j,1) - (HALF + dth*vel(i  ,j,1)/hx) * velx(i  ,j,1)
!    $           + dth * force(i  ,j,1)
          ulft = vel(i-1,j,1) + (HALF - dth*vel(i-1,j,1)/hx) * velx(i-1,j,1)
!    $           + dth * force(i-1,j,1)

          urgt = merge(vel(is-1,j,1),urgt,i.eq.is   .and. bc(1,1) .eq. INLET)
          urgt = merge(vel(ie+1,j,1),urgt,i.eq.ie+1 .and. bc(1,2) .eq. INLET)
          urgt = merge(ZERO     ,urgt,i.eq.is   .and. bc(1,1) .eq. WALL)
          urgt = merge(ZERO     ,urgt,i.eq.ie+1 .and. bc(1,2) .eq. WALL)

          ulft = merge(vel(is-1,j,1),ulft,i.eq.is   .and. bc(1,1) .eq. INLET)
          ulft = merge(vel(ie+1,j,1),ulft,i.eq.ie+1 .and. bc(1,2) .eq. INLET)
          ulft = merge(ZERO     ,ulft,i.eq.is   .and. bc(1,1) .eq. WALL)
          ulft = merge(ZERO     ,ulft,i.eq.ie+1 .and. bc(1,2) .eq. WALL)

          utrans(i,j,1) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
          test = ( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or.  &
                   (abs(ulft+urgt) .lt. eps) )
          utrans(i,j,1) = merge(ZERO,utrans(i,j,1),test)

        enddo
      enddo

!     Create the y-velocity to be used for transverse derivatives.
      do j = js,je+1 
        do i = is-1,ie+1 

          vtop = vel(i,j  ,2) - (HALF + dth*vel(i,j  ,2)/hy) * vely(i,j  ,2)
!    $           + dth * force(i,j  ,2)
          vbot = vel(i,j-1,2) + (HALF - dth*vel(i,j-1,2)/hy) * vely(i,j-1,2)
!    $           + dth * force(i,j-1,2)

          vtop = merge(vel(i,js-1,2),vtop,j.eq.js   .and. bc(2,1) .eq. INLET)
          vtop = merge(vel(i,je+1,2),vtop,j.eq.je+1 .and. bc(2,2) .eq. INLET)
          vtop = merge(ZERO     ,vtop,j.eq.js   .and. bc(2,1) .eq. WALL)
          vtop = merge(ZERO     ,vtop,j.eq.je+1 .and. bc(2,2) .eq. WALL)

          vbot = merge(vel(i,js-1,2),vbot,j.eq.js   .and. bc(2,1) .eq. INLET)
          vbot = merge(vel(i,je+1,2),vbot,j.eq.je+1 .and. bc(2,2) .eq. INLET)
          vbot = merge(ZERO     ,vbot,j.eq.js   .and. bc(2,1) .eq. WALL)
          vbot = merge(ZERO     ,vbot,j.eq.je+1 .and. bc(2,2) .eq. WALL)

          utrans(i,j,2)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
          test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or.  &
                   (abs(vbot+vtop) .lt. eps))
          utrans(i,j,2) = merge(ZERO,utrans(i,j,2),test)
        enddo
      enddo

      end subroutine mkutrans_2d

      subroutine mkutrans_3d(vel,utrans,force,lo,dx,dt,ng,bc)

      integer, intent(in) :: lo(3),ng

      real(kind=dp_t), intent(in   ) ::    vel(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real(kind=dp_t), intent(inout) :: utrans(lo(1)- 2:,lo(2)- 2:,lo(3)- 2:,:)

      real(kind=dp_t),intent(in) :: dt,dx(3)
      integer        ,intent(in) :: bc(3,2)

      real(kind=dp_t), allocatable::  velx(:,:,:,:),vely(:,:,:,:),velz(:,:,:,:)

!     Local variables
      real(kind=dp_t) ulft,urgt,vbot,vtop,wbot,wtop
      real(kind=dp_t) hx, hy, hz, dth

      real(kind=dp_t) :: eps

      logical :: test
      logical :: do_refl

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke

      integer :: slope_order = 4

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(vel,dim=2) - (2*ng+1)
 
      is = lo(1)
      js = lo(2)
      ks = lo(3)
      ie = hi(1)
      je = hi(2)
      ke = hi(3)

      eps = 1.0e-8              ! FIXME what should EPS really be?

      dth = HALF * dt

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(velz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

      do k = lo(3),hi(3)
         call slopex_2d(vel(:,:,k,:),velx(:,:,k,:),lo,ng,3,bc,do_refl,slope_order)
         call slopey_2d(vel(:,:,k,:),vely(:,:,k,:),lo,ng,3,bc,        slope_order)
      end do
      call slopez_3d(vel,velz,lo,ng,  3,bc,    slope_order)

!     Create the x-velocity to be used for transverse derivatives.
      do k = ks-1,ke+1
      do j = js-1,je+1
        do i = is,ie+1

          urgt = vel(i,j,k  ,1) - (HALF + dth*vel(i  ,j,k,1)/hx) * velx(i  ,j,k,1)
!    $           + dth * force(i  ,j,k,1)
          ulft = vel(i-1,j,k,1) + (HALF - dth*vel(i-1,j,k,1)/hx) * velx(i-1,j,k,1)
!    $           + dth * force(i-1,j,k,1)

          urgt = merge(vel(is-1,j,k,1),urgt,i.eq.is   .and. bc(1,1) .eq. INLET)
          urgt = merge(vel(ie+1,j,k,1),urgt,i.eq.ie+1 .and. bc(1,2) .eq. INLET)
          urgt = merge(ZERO           ,urgt,i.eq.is   .and. bc(1,1) .eq. WALL)
          urgt = merge(ZERO           ,urgt,i.eq.ie+1 .and. bc(1,2) .eq. WALL)

          ulft = merge(vel(is-1,j,k,1),ulft,i.eq.is   .and. bc(1,1) .eq. INLET)
          ulft = merge(vel(ie+1,j,k,1),ulft,i.eq.ie+1 .and. bc(1,2) .eq. INLET)
          ulft = merge(ZERO           ,ulft,i.eq.is   .and. bc(1,1) .eq. WALL)
          ulft = merge(ZERO           ,ulft,i.eq.ie+1 .and. bc(1,2) .eq. WALL)

          utrans(i,j,k,1) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
          test=( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or. &
                (abs(ulft+urgt) .lt. eps) )
          utrans(i,j,k,1) = merge(ZERO,utrans(i,j,k,1),test)

        enddo
        enddo
      enddo

!     Create the y-velocity to be used for transverse derivatives.
      do j = js,je+1
        do k = ks-1,ke+1
        do i = is-1,ie+1

          vtop = vel(i,j  ,k,2) - (HALF + dth*vel(i,j  ,k,2)/hy) * vely(i,j  ,k,2)
!    $           + dth * force(i,j  ,k,2)
          vbot = vel(i,j-1,k,2) + (HALF - dth*vel(i,j-1,k,2)/hy) * vely(i,j-1,k,2)
!    $           + dth * force(i,j-1,k,2)

          vtop = merge(vel(i,js-1,k,2),vtop,j.eq.js   .and. bc(2,1) .eq. INLET)
          vtop = merge(vel(i,je+1,k,2),vtop,j.eq.je+1 .and. bc(2,2) .eq. INLET)
          vtop = merge(ZERO           ,vtop,j.eq.js   .and. bc(2,1) .eq. WALL)
          vtop = merge(ZERO           ,vtop,j.eq.je+1 .and. bc(2,2) .eq. WALL)

          vbot = merge(vel(i,js-1,k,2),vbot,j.eq.js   .and. bc(2,1) .eq. INLET)
          vbot = merge(vel(i,je+1,k,2),vbot,j.eq.je+1 .and. bc(2,2) .eq. INLET)
          vbot = merge(ZERO           ,vbot,j.eq.js   .and. bc(2,1) .eq. WALL)
          vbot = merge(ZERO           ,vbot,j.eq.je+1 .and. bc(2,2) .eq. WALL)

          utrans(i,j,k,2)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
          test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or. &
                   (abs(vbot+vtop) .lt. eps))
          utrans(i,j,k,2) = merge(ZERO,utrans(i,j,k,2),test)

        enddo
        enddo
      enddo

!     Create the z-velocity to be used for transverse derivatives.
      do k = ks,ke+1
        do j = js-1,je+1
        do i = is-1,ie+1

          wtop = vel(i,j,k  ,3) - (HALF + dth*vel(i,j,k  ,3)/hz) * velz(i,j,k  ,3)
!    $           + dth * force(i,j,k  ,3)
          wbot = vel(i,j,k-1,3) + (HALF - dth*vel(i,j,k-1,3)/hz) * velz(i,j,k-1,3)
!    $           + dth * force(i,j,k-1,3)

          wtop = merge(vel(i,j,ks-1,3),wtop,k.eq.ks   .and. bc(3,1) .eq. INLET)
          wtop = merge(vel(i,j,ke+1,3),wtop,k.eq.ke+1 .and. bc(3,2) .eq. INLET)
          wtop = merge(ZERO           ,wtop,k.eq.ks   .and. bc(3,1) .eq. WALL)
          wtop = merge(ZERO           ,wtop,k.eq.ke+1 .and. bc(3,2) .eq. WALL)

          wbot = merge(vel(i,j,ks-1,3),wbot,k.eq.ks   .and. bc(3,1) .eq. INLET)
          wbot = merge(vel(i,j,ke+1,3),wbot,k.eq.ke+1 .and. bc(3,2) .eq. INLET)
          wbot = merge(ZERO           ,wbot,k.eq.ks   .and. bc(3,1) .eq. WALL)
          wbot = merge(ZERO           ,wbot,k.eq.ke+1 .and. bc(3,2) .eq. WALL)

          utrans(i,j,k,3)=merge(wbot,wtop,(wbot+wtop).gt.ZERO)
          test = ( (wbot .le. ZERO  .and.  wtop .ge. ZERO)  .or. &
                   (abs(wbot+wtop) .lt. eps))
          utrans(i,j,k,3) = merge(ZERO,utrans(i,j,k,3),test)

        enddo
        enddo
      enddo

      end subroutine mkutrans_3d

end module mkutrans_module
