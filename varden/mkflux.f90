module mkflux_module

  use bl_types
  use multifab_module
  use slope_module
  use cvmg_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

contains

      subroutine mkflux_2d(s,u,sedgex,sedgey,&
                           uadv,utrans,force,lo,dx,dt,is_vel, &
                           visc_coef,irz,bc,velpred,ng)

      integer, intent(in) :: lo(2),ng

      real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), intent(inout) ::   uadv(lo(1)-1:,lo(2)-1:,:)
      real(kind=dp_t), intent(in   ) :: utrans(lo(1)-2:,lo(2)-2:,:)

      real(kind=dp_t),intent(in) :: dt,dx(2),visc_coef
      integer        ,intent(in) :: irz,velpred,bc(2,2)
      logical        ,intent(in) :: is_vel

      real(kind=dp_t), allocatable::  slopex(:,:,:),slopey(:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:)

!     Local variables
      real(kind=dp_t) ubardth, vbardth
      real(kind=dp_t) hx, hy, dth
      real(kind=dp_t) splus,sminus
      real(kind=dp_t) savg,st
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt

      integer :: hi(2)
      integer :: n,ncomp
      integer :: slope_order = 4
      logical :: test

      real(kind=dp_t) :: eps

      integer :: i,j,is,js,ie,je
      logical :: do_refl

      ncomp = size(s,dim=3)
      do_refl = irz .eq. 1

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      allocate(s_l(lo(1)-1:hi(1)+2))
      allocate(s_r(lo(1)-1:hi(1)+2))
      allocate(s_b(lo(2)-1:hi(2)+2))
      allocate(s_t(lo(2)-1:hi(2)+2))

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))

      call slopex_2d(s,slopex,lo,ng,ncomp,bc,do_refl,slope_order)
      call slopey_2d(s,slopey,lo,ng,ncomp,bc,        slope_order)

      eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)

      dth = HALF*dt

      hx = dx(1)
      hy = dx(2)

!
!     Loop for fluxes on x-edges.
!
      do n = 1,ncomp
       do j = js,je 
        if (velpred .eq. 0 .or. n .eq. 1) then
        do i = is-1,ie+1 

          spbot = s(i,j  ,n) + (HALF - dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
!    $            + dth * force(i,  j,n)
          sptop = s(i,j+1,n) - (HALF + dth*u(i,j+1,2)/hy) * slopey(i,j+1,n)
!    $            + dth * force(i,j+1,n)

          sptop = cvmgt(s(i,je+1,n),sptop,j.eq.je .and. bc(2,2) .eq. INLET)
          spbot = cvmgt(s(i,je+1,n),spbot,j.eq.je .and. bc(2,2) .eq. INLET)

          if (j .eq. je .and. bc(2,2) .eq. WALL) then
            if (is_vel .and. n .eq. 2) then
              sptop = ZERO
              spbot = ZERO
            elseif (is_vel .and. n .eq. 1) then
              sptop = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
              spbot = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
            else
              sptop = spbot
            endif
          endif

          splus = cvmgp(spbot,sptop,utrans(i,j+1,2))
          savg  = HALF * (spbot + sptop)
          splus = cvmgt(splus, savg, abs(utrans(i,j+1,2)) .gt. eps)

          smtop = s(i,j  ,n) - (HALF + dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
!    $            + dth * force(i,j  ,n)
          smbot = s(i,j-1,n) + (HALF - dth*u(i,j-1,2)/hy) * slopey(i,j-1,n)
!    $            + dth * force(i,j-1,n)

          smtop = cvmgt(s(i,js-1,n),smtop,j.eq.js .and. bc(2,1) .eq. INLET)
          smbot = cvmgt(s(i,js-1,n),smbot,j.eq.js .and. bc(2,1) .eq. INLET)

          if (j .eq. js .and. bc(2,1) .eq. WALL) then
            if (is_vel .and. (n .eq. 2)) then
              smtop = ZERO
              smbot = ZERO
            elseif (is_vel .and. (n .ne. 2)) then
              smbot = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
              smtop = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
            else
              smbot = smtop
            endif
          endif

          sminus = cvmgp(smbot,smtop,utrans(i,j,2))
          savg   = HALF * (smbot + smtop)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,2)) .gt. eps)

          st = force(i,j,n) - &
                HALF * (utrans(i,j,2)+utrans(i,j+1,2))*(splus - sminus) / hy

          ubardth = dth*u(i,j,1)/hx

          s_l(i+1)= s(i,j,n) + (HALF-ubardth)*slopex(i,j,n) + dth*st
          s_r(i  )= s(i,j,n) - (HALF+ubardth)*slopex(i,j,n) + dth*st

         enddo

         if (velpred .eq. 1) then
           do i = is, ie+1 
             savg = HALF*(s_r(i) + s_l(i))
             test = ( (s_l(i) .le. ZERO  .and. &
                       s_r(i) .ge. ZERO)  .or. &
                     (abs(s_l(i) + s_r(i)) .lt. eps) )
             sedgex(i,j,n)=cvmgp(s_l(i),s_r(i),savg)
             sedgex(i,j,n)=cvmgt(savg,sedgex(i,j,n),test)
           enddo
         else
           do i = is, ie+1 
             sedgex(i,j,n)=cvmgp(s_l(i),s_r(i),uadv(i,j,1))
             savg = HALF*(s_r(i) + s_l(i))
             sedgex(i,j,n)=cvmgt(savg,sedgex(i,j,n),abs(uadv(i,j,1)) .lt. eps)
           enddo
         endif

         if (bc(1,1) .eq. WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(is,j,n) = ZERO
           elseif (is_vel .and. n .ne. 1) then
             sedgex(is,j,n) = cvmgt(ZERO,s_r(is),visc_coef.gt.0.0.and.irz.eq.0)
           else 
             sedgex(is,j,n) = s_r(is)
           endif
         elseif (bc(1,1) .eq. INLET) then
           sedgex(is,j,n) = s(is-1,j,n)
         elseif (bc(1,1) .eq. OUTLET) then
           sedgex(is,j,n) = s_r(is)
         endif
         if (bc(1,2) .eq. WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(ie+1,j,n) = ZERO
           else if (is_vel .and. n .ne. 1) then
             sedgex(ie+1,j,n) = cvmgt(ZERO,s_l(ie+1),visc_coef .gt. 0.0)
           else 
             sedgex(ie+1,j,n) = s_l(ie+1)
           endif
         elseif (bc(1,2) .eq. INLET) then
           sedgex(ie+1,j,n) = s(ie+1,j,n)
         elseif (bc(1,2) .eq. OUTLET) then
           sedgex(ie+1,j,n) = s_l(ie+1)
         endif

         if (velpred .eq. 1) then
           do i = is, ie+1 
             uadv(i,j,1) = sedgex(i,j,1)
           enddo
         endif
         endif
       enddo

      enddo

!
!     Loop for fluxes on y-edges.
!
      do n = 1,ncomp
       do i = is, ie 
        if (velpred .eq. 0 .or. n .eq. 2) then
        do j = js-1, je+1 

          splft = s(i,j  ,n) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
!    $            + dth * force(i  ,j,n)
          sprgt = s(i+1,j,n) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,n)
!    $            + dth * force(i+1,j,n)

          sprgt = cvmgt(s(ie+1,j,n),sprgt,i.eq.ie .and. bc(1,2) .eq. INLET)
          splft = cvmgt(s(ie+1,j,n),splft,i.eq.ie .and. bc(1,2) .eq. INLET)

          if (i .eq. ie .and. bc(1,2) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              splft = ZERO
              sprgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              sprgt = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
              splft = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
            else
              sprgt = splft
            endif
          endif

          splus = cvmgp(splft,sprgt,utrans(i+1,j,1))
          savg  = HALF * (splft + sprgt)
          splus = cvmgt(splus, savg, abs(utrans(i+1,j,1)) .gt. eps)

          smrgt = s(i  ,j,n) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
!    $            + dth * force(i  ,j,n)
          smlft = s(i-1,j,n) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,n)
!    $            + dth * force(i-1,j,n)

          smrgt = cvmgt(s(is-1,j,n),smrgt,i.eq.is .and. bc(1,1) .eq. INLET)
          smlft = cvmgt(s(is-1,j,n),smlft,i.eq.is .and. bc(1,1) .eq. INLET)

          if (i .eq. is .and. bc(1,1) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              smlft = ZERO
              smrgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              smlft = cvmgt(ZERO,smrgt,visc_coef.gt.ZERO.and.irz.eq.0)
              smrgt = cvmgt(ZERO,smrgt,visc_coef.gt.ZERO.and.irz.eq.0)
            else
              smlft = smrgt
            endif
          endif

          sminus = cvmgp(smlft,smrgt,utrans(i,j,1))
          savg   = HALF * (smlft + smrgt)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,1)) .gt. eps)

          st = force(i,j,n) - &
               HALF * (utrans(i,j,1)+utrans(i+1,j,1))*(splus - sminus) / hx

          vbardth = dth*u(i,j,2)/hy

          s_b(j+1)= s(i,j,n) + (HALF-vbardth)*slopey(i,j,n) + dth*st
          s_t(j  )= s(i,j,n) - (HALF+vbardth)*slopey(i,j,n) + dth*st

        enddo

        if (velpred .eq. 1) then
          do j = js, je+1 
            savg = HALF*(s_b(j) + s_t(j))
            test = ( (s_b(j) .le. ZERO  .and. &
                      s_t(j) .ge. ZERO)  .or. &
                   (abs(s_b(j) + s_t(j)) .lt. eps) )
            sedgey(i,j,n)=cvmgp(s_b(j),s_t(j),savg)
            sedgey(i,j,n)=cvmgt(savg,sedgey(i,j,n),test)
          enddo
        else

          do j = js, je+1 
            sedgey(i,j,n)=cvmgp(s_b(j),s_t(j),uadv(i,j,2))
            savg = HALF*(s_b(j) + s_t(j))
            sedgey(i,j,n)=cvmgt(savg,sedgey(i,j,n),abs(uadv(i,j,2)) .lt. eps)
          enddo
        endif

        if (bc(2,1) .eq. WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,js,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,js,n) = cvmgt(ZERO,s_t(js),visc_coef .gt. 0.0)
          else 
            sedgey(i,js,n) = s_t(js)
          endif
        elseif (bc(2,1) .eq. INLET) then
          sedgey(i,js,n) = s(i,js-1,n)
        elseif (bc(2,1) .eq. OUTLET) then
          sedgey(i,js,n) = s_t(js)
        endif

        if (bc(2,2) .eq. WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,je+1,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,je+1,n) = cvmgt(ZERO,s_b(je+1),visc_coef .gt. 0.0)
          else 
            sedgey(i,je+1,n) = s_b(je+1)
          endif
        elseif (bc(2,2) .eq. INLET) then
          sedgey(i,je+1,n) = s(i,je+1,n)
        elseif (bc(2,2) .eq. OUTLET) then
          sedgey(i,je+1,n) = s_b(je+1)
        endif

        if (velpred .eq. 1) then
          do j = js, je+1 
            uadv(i,j,2) = sedgey(i,j,2)
          enddo
        endif
        endif
       enddo
      enddo

      end subroutine mkflux_2d

      subroutine mkflux_3d(s,u,sedgex,sedgey,sedgez,&
                           uadv,utrans,force,lo,dx,dt,is_vel, &
                           visc_coef,bc,velpred,ng)


      implicit none

      integer, intent(in) :: lo(3),ng

      real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::  force(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
      real(kind=dp_t),intent(inout) :: sedgex(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
      real(kind=dp_t),intent(inout) :: sedgey(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
      real(kind=dp_t),intent(inout) :: sedgez(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
      real(kind=dp_t),intent(inout) ::   uadv(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
      real(kind=dp_t),intent(in   ) :: utrans(lo(1)-2:,lo(2)-2:,lo(3)-2:,:)

      real(kind=dp_t),intent(in) :: dt,dx(3),visc_coef
      integer        ,intent(in) :: velpred,bc(3,2)
      logical        ,intent(in) :: is_vel

      real(kind=dp_t), allocatable::  slopex(:,:,:,:)
      real(kind=dp_t), allocatable::  slopey(:,:,:,:)
      real(kind=dp_t), allocatable::  slopez(:,:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)

!     Local variables
      real(kind=dp_t) ubardth, vbardth, wbardth
      real(kind=dp_t) hx, hy, hz, dth
      real(kind=dp_t) splus,sminus,st,str,savg
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
      logical test

      real(kind=dp_t) :: eps

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke,n
      integer :: slope_order = 4
      integer :: ncomp
      logical :: do_refl = .false.

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

      allocate(s_l(lo(1)-1:hi(1)+2))
      allocate(s_r(lo(1)-1:hi(1)+2))
      allocate(s_b(lo(2)-1:hi(2)+2))
      allocate(s_t(lo(2)-1:hi(2)+2))
      allocate(s_u(lo(3)-1:hi(3)+2))
      allocate(s_d(lo(3)-1:hi(3)+2))

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

      ncomp = size(s,dim=4)
      do k = lo(3),hi(3)
         call slopex_2d(s(:,:,k,:),slopex(:,:,k,:),lo,ng,ncomp,bc,do_refl,slope_order)
         call slopey_2d(s(:,:,k,:),slopey(:,:,k,:),lo,ng,ncomp,bc,        slope_order)
      end do
      call slopez_3d(s,slopez,lo,ng,ncomp,bc,slope_order)

      eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

      dth = HALF*dt

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

!     loop for x fluxes
      do n = 1,ncomp
       if (velpred .eq. 0 .or. n .eq. 1) then
        do k = ks,ke 
        do j = js,je 
         do i = is-1,ie+1 

!        ******************************************************************
!         MAKE TRANSVERSE DERIVATIVES IN Y-DIRECTION
!        ******************************************************************

          spbot = s(i,j  ,k,n) + (HALF - dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
!    $            + dth * force(i,  j,k,n)
          sptop = s(i,j+1,k,n) - (HALF + dth*u(i,j+1,k,2)/hy)*slopey(i,j+1,k,n)
!    $            + dth * force(i,j+1,k,n)

          sptop = cvmgt(s(i,je+1,k,n),sptop,j.eq.je .and. bc(2,2) .eq. INLET)
          spbot = cvmgt(s(i,je+1,k,n),spbot,j.eq.je .and. bc(2,2) .eq. INLET)

          if (j .eq. je .and. bc(2,2) .eq. WALL) then
            if (is_vel .and. n .eq. 2) then
              sptop = ZERO
              spbot = ZERO
            else if (is_vel .and. n .ne. 2) then
              sptop = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
              spbot = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
            else
              sptop = spbot
            endif
          endif

          splus = cvmgp(spbot,sptop,utrans(i,j+1,k,2))
          savg  = HALF * (spbot + sptop)
          splus = cvmgt(splus, savg, abs(utrans(i,j+1,k,2)) .gt. eps)

          smtop = s(i,j  ,k,n) - (HALF + dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
!    $            + dth * force(i,j  ,k,n)
          smbot = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n)
!    $            + dth * force(i,j-1,k,n)

          smtop = cvmgt(s(i,js-1,k,n),smtop,j.eq.js .and. bc(2,1) .eq. INLET)
          smbot = cvmgt(s(i,js-1,k,n),smbot,j.eq.js .and. bc(2,1) .eq. INLET)

          if (j .eq. js .and. bc(2,1) .eq. WALL) then
            if (is_vel .and. n .eq. 2) then
              smtop = ZERO
              smbot = ZERO
            else if (is_vel .and. n .ne. 2) then
              smtop = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
              smbot = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
            else
              smbot = smtop
            endif
          endif

          sminus = cvmgp(smbot,smtop,utrans(i,j,k,2))
          savg   = HALF * (smbot + smtop)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,2)) .gt. eps)

          str =  HALF * (utrans(i,j,k,2)+utrans(i,j+1,k,2))*(splus - sminus) / hy

!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Z-DIRECTION
!        ******************************************************************

          spbot = s(i,j,k  ,n) + (HALF - dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
!    $            + dth * force(i,j,k  ,n)
          sptop = s(i,j,k+1,n) - (HALF + dth*u(i,j,k+1,3)/hz)*slopez(i,j,k+1,n)
!    $            + dth * force(i,j,k+1,n)

          sptop = cvmgt(s(i,j,ke+1,n),sptop,k.eq.ke .and. bc(3,2) .eq. INLET)
          spbot = cvmgt(s(i,j,ke+1,n),spbot,k.eq.ke .and. bc(3,2) .eq. INLET)

          if (k .eq. ke .and. bc(3,2) .eq. WALL) then
            if (is_vel .and. n .eq. 3) then
              sptop = ZERO
              spbot = ZERO
            else if (is_vel .and. n .ne. 3) then
              sptop = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
              spbot = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
            else
              sptop = spbot
            endif
          endif

          splus = cvmgp(spbot,sptop,utrans(i,j,k+1,3))
          savg  = HALF * (spbot + sptop)
          splus = cvmgt(splus, savg, abs(utrans(i,j,k+1,3)) .gt. eps)

          smtop = s(i,j,k  ,n) - (HALF + dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
!    $            + dth * force(i,j,k  ,n)
          smbot = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n)
!    $            + dth * force(i,j,k-1,n)

          smtop = cvmgt(s(i,j,ks-1,n),smtop,k.eq.ks .and. bc(3,1) .eq. INLET)
          smbot = cvmgt(s(i,j,ks-1,n),smbot,k.eq.ks .and. bc(3,1) .eq. INLET)

          if (k .eq. ks .and. bc(3,1) .eq. WALL) then
            if (is_vel .and. n .eq. 3) then
              smtop = ZERO
              smbot = ZERO
            else if (is_vel .and. n .ne. 3) then
              smtop = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
              smbot = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
            else 
              smbot = smtop
            endif
          endif

          sminus = cvmgp(smbot,smtop,utrans(i,j,k,3))
          savg   = HALF * (smbot + smtop)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,3)) .gt. eps)

          str = str + HALF * (utrans(i,j,k,3)+utrans(i,j,k+1,3))* &
                             (splus - sminus) / hz

!        ******************************************************************
!        MAKE LEFT AND RIGHT STATES
!        ******************************************************************

          st = force(i,j,k,n) - str
          ubardth = dth*u(i,j,k,1)/hx

          s_l(i+1)= s(i,j,k,n) + (HALF-ubardth)*slopex(i,j,k,n) + dth*st
          s_r(i  )= s(i,j,k,n) - (HALF+ubardth)*slopex(i,j,k,n) + dth*st

        enddo

        if (velpred .eq. 1) then
          do i = is, ie+1 
            savg = HALF*(s_r(i) + s_l(i))
            test = ( (s_l(i) .le. ZERO  .and. &
                      s_r(i) .ge. ZERO)  .or. &
                    (abs(s_l(i) + s_r(i)) .lt. eps) )
            sedgex(i,j,k,n)=cvmgp(s_l(i),s_r(i),savg)
            sedgex(i,j,k,n)=cvmgt(savg,sedgex(i,j,k,n),test)
          enddo
        else
          do i = is, ie+1 
            sedgex(i,j,k,n)=cvmgp(s_l(i),s_r(i),uadv(i,j,k,1))
            savg = HALF*(s_r(i) + s_l(i))
            sedgex(i,j,k,n)=cvmgt(savg,sedgex(i,j,k,n),abs(uadv(i,j,k,1)) .lt. eps)
          enddo
        endif

        if (bc(1,1) .eq. WALL) then
          if (is_vel .and. n .eq. 1) then
            sedgex(is,j,k,n) = ZERO
          else if (is_vel .and. n .ne. 1) then
            sedgex(is,j,k,n) = cvmgt(ZERO,s_r(is),visc_coef .gt. 0.0)
          else 
            sedgex(is,j,k,n) = s_r(is)
          endif
        elseif (bc(1,1) .eq. INLET) then
          sedgex(is,j,k,n) = s(is-1,j,k,n)
        elseif (bc(1,1) .eq. OUTLET) then
          sedgex(is,j,k,n) = s_r(is)
        endif

        if (bc(1,2) .eq. WALL) then
          if (is_vel .and. n .eq. 1) then
            sedgex(ie+1,j,k,n) = ZERO
          else if (is_vel .and. n .ne. 1) then
            sedgex(ie+1,j,k,n) = cvmgt(ZERO,s_l(ie+1),visc_coef .gt. 0.0)
          else 
            sedgex(ie+1,j,k,n) = s_l(ie+1)
          endif
        elseif (bc(1,2) .eq. INLET) then
          sedgex(ie+1,j,k,n) = s(ie+1,j,k,n)
        elseif (bc(1,2) .eq. OUTLET) then
          sedgex(ie+1,j,k,n) = s_l(ie+1)
        endif

        if (velpred .eq. 1) then
          do i = is, ie+1 
            uadv(i,j,k,1) = sedgex(i,j,k,1)
          enddo
        endif

        enddo
        enddo
       endif
      enddo

!        ******************************************************************
!        ******************************************************************
!        ******************************************************************

!     loop for y fluxes

      do n = 1,ncomp
       if (velpred .eq. 0 .or. n .eq. 2) then
       do k = ks, ke 
       do i = is, ie 
        do j = js-1, je+1 

!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN X-DIRECTION
!        ******************************************************************

          splft = s(i  ,j,k,n) + (HALF - dth*u(i  ,j,k,1)/hx)*slopex(i  ,j,k,n)
!    $            + dth * force(i  ,j,k,n)
          sprgt = s(i+1,j,k,n) - (HALF + dth*u(i+1,j,k,1)/hx)*slopex(i+1,j,k,n)
!    $            + dth * force(i+1,j,k,n)

          sprgt = cvmgt(s(ie+1,j,k,n),sprgt,i.eq.ie .and. bc(1,2) .eq. INLET)
          splft = cvmgt(s(ie+1,j,k,n),splft,i.eq.ie .and. bc(1,2) .eq. INLET)

          if (i .eq. ie .and. bc(1,2) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              sprgt = ZERO
              splft = ZERO
            else if (is_vel .and. n .ne. 1) then
              sprgt = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
              splft = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
            else
              sprgt = splft
            endif
          endif

          splus = cvmgp(splft,sprgt,utrans(i+1,j,k,1))
          savg  = HALF * (splft + sprgt)
          splus = cvmgt(splus, savg, abs(utrans(i+1,j,k,1)) .gt. eps)

          smrgt = s(i  ,j,k,n) - (HALF + dth*u(i  ,j,k,1)/hx)*slopex(i  ,j,k,n)
!    $            + dth * force(i  ,j,k,n)
          smlft = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n)
!    $            + dth * force(i-1,j,k,n)

          smrgt = cvmgt(s(is-1,j,k,n),smrgt,i.eq.is .and. bc(1,1) .eq. INLET)
          smlft = cvmgt(s(is-1,j,k,n),smlft,i.eq.is .and. bc(1,1) .eq. INLET)

          if (i .eq. is .and. bc(1,1) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              smrgt = ZERO
              smlft = ZERO
            else if (is_vel .and. n .ne. 1) then
              smrgt = cvmgt(ZERO,smrgt,visc_coef .gt. ZERO)
              smlft = cvmgt(ZERO,smrgt,visc_coef .gt. ZERO)
            else
              smlft = smrgt
            endif
          endif
 
          sminus = cvmgp(smlft,smrgt,utrans(i,j,k,1))
          savg   = HALF * (smlft + smrgt)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,1)) .gt. eps)

          str    = HALF * (utrans(i,j,k,1)+utrans(i+1,j,k,1))*(splus - sminus) / hx

!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Z-DIRECTION
!        ******************************************************************

          spbot = s(i,j,k  ,n) + (HALF - dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
!    $            + dth * force(i,j,k  ,n)
          sptop = s(i,j,k+1,n) - (HALF + dth*u(i,j,k+1,3)/hz)*slopez(i,j,k+1,n)
!    $            + dth * force(i,j,k+1,n)

          sptop = cvmgt(s(i,j,ke+1,n),sptop,k.eq.ke .and. bc(3,2) .eq. INLET)
          spbot = cvmgt(s(i,j,ke+1,n),spbot,k.eq.ke .and. bc(3,2) .eq. INLET)

          if (k .eq. ke .and. bc(3,2) .eq. WALL) then
            if (is_vel .and. n .eq. 3) then
              sptop = ZERO
              spbot = ZERO
            else if (is_vel .and. n .ne. 3) then
              sptop = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
              spbot = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
            else
              sptop = spbot
            endif
          endif

          splus = cvmgp(spbot,sptop,utrans(i,j,k+1,3))
          savg  = HALF * (spbot + sptop)
          splus = cvmgt(splus, savg, abs(utrans(i,j,k+1,3)) .gt. eps)

          smtop = s(i,j,k  ,n) - (HALF + dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
!    $            + dth * force(i,j,k  ,n)
          smbot = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n)
!    $            + dth * force(i,j,k-1,n)

          smtop = cvmgt(s(i,j,ks-1,n),smtop,k.eq.ke .and. bc(3,1) .eq. INLET)
          smbot = cvmgt(s(i,j,ks-1,n),smbot,k.eq.ke .and. bc(3,1) .eq. INLET)

          if (k .eq. ks  .and.  bc(3,1) .eq. WALL) then
            if (is_vel .and. n .eq. 3) then
              smtop = ZERO
              smbot = ZERO
            else if (is_vel .and. n .ne. 3) then
              smbot = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
              smtop = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
            else
              smbot = smtop
            endif
          endif

          sminus = cvmgp(smbot,smtop,utrans(i,j,k,3))
          savg   = HALF * (smbot + smtop)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,3)) .gt. eps)

          str = str + HALF * (utrans(i,j,k,3)+utrans(i,j,k+1,3))*(splus - sminus) / hz

!        ******************************************************************
!        MAKE TOP AND BOTTOM STATES
!        ******************************************************************

          st = force(i,j,k,n) - str

          vbardth = dth*u(i,j,k,2)/hy

          s_b(j+1)= s(i,j,k,n) + (HALF-vbardth)*slopey(i,j,k,n) + dth*st
          s_t(j  )= s(i,j,k,n) - (HALF+vbardth)*slopey(i,j,k,n) + dth*st
        enddo

        if (velpred .eq. 1) then
          do j = js, je+1 
            savg = HALF*(s_t(j) + s_b(j))
            test = ( (s_b(j) .le. ZERO  .and. &
                      s_t(j) .ge. ZERO)  .or. &
                   (abs(s_b(j) + s_t(j)) .lt. eps) )
            sedgey(i,j,k,n)=cvmgp(s_b(j),s_t(j),savg)
            sedgey(i,j,k,n)=cvmgt(savg,sedgey(i,j,k,n),test)
          enddo
        else
          do j = js, je+1 
            sedgey(i,j,k,n)=cvmgp(s_b(j),s_t(j),uadv(i,j,k,2))
            savg = HALF*(s_t(j) + s_b(j))
            sedgey(i,j,k,n)=cvmgt(savg,sedgey(i,j,k,n),abs(uadv(i,j,k,2)) .lt. eps)
          enddo
        endif

        if (bc(2,1) .eq. WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,js,k,n) = ZERO
          else if (is_vel .and. n .ne. 2) then
            sedgey(i,js,k,n) = cvmgt(ZERO,s_t(js),visc_coef .gt. 0.0)
          else
            sedgey(i,js,k,n) = s_t(js)
          endif
        elseif (bc(2,1) .eq. INLET) then
          sedgey(i,js,k,n) = s(i,js-1,k,n)
        elseif (bc(2,1) .eq. OUTLET) then
          sedgey(i,js,k,n) = s_t(js)
        endif

        if (bc(2,2) .eq. WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,je+1,k,n) = ZERO
          else if (is_vel .and. n .ne. 2) then
            sedgey(i,je+1,k,n) = cvmgt(ZERO,s_b(je+1),visc_coef .gt. 0.0)
          else
            sedgey(i,je+1,k,n) = s_b(je+1)
          endif
        elseif (bc(2,2) .eq. INLET) then
          sedgey(i,je+1,k,n) = s(i,je+1,k,n)
        elseif (bc(2,2) .eq. OUTLET) then
          sedgey(i,je+1,k,n) = s_b(je+1)
        endif

        if (velpred .eq. 1) then
          do j = js, je+1 
            uadv(i,j,k,2) = sedgey(i,j,k,2)
          enddo
        endif

        enddo
        enddo
       endif
      enddo

!     ******************************************************************
!     ******************************************************************
!     ******************************************************************

!     loop for z fluxes

      do n = 1,ncomp
        if (velpred .eq. 0 .or. n .eq. 3) then
        do j = js, je 
        do i = is, ie 
          do k = ks-1, ke+1 

!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN X-DIRECTION
!        ******************************************************************

          splft = s(i  ,j,k,n) + (HALF - dth*u(i  ,j,k,1)/hx) * slopex(i  ,j,k,n)
!    $            + dth * force(i  ,j,k,n)
          sprgt = s(i+1,j,k,n) - (HALF + dth*u(i+1,j,k,1)/hx) * slopex(i+1,j,k,n)
!    $            + dth * force(i+1,j,k,n)

          sprgt = cvmgt(s(ie+1,j,k,n),sprgt,i.eq.ie .and. bc(1,2) .eq. INLET)
          splft = cvmgt(s(ie+1,j,k,n),splft,i.eq.ie .and. bc(1,2) .eq. INLET)

          if (i .eq. ie  .and.  bc(1,2) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              sprgt = ZERO
              splft = ZERO
            else if (is_vel .and. n .ne. 1) then
              sprgt = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
              splft = cvmgt(ZERO,splft,visc_coef .gt. ZERO)
            else
              sprgt = splft
            endif
          endif

          splus = cvmgp(splft,sprgt,utrans(i+1,j,k,1))
          savg  = HALF * (splft + sprgt)
          splus = cvmgt(splus, savg, abs(utrans(i+1,j,k,1)) .gt. eps)

          smrgt = s(i  ,j,k,n) - (HALF + dth*u(i  ,j,k,1)/hx)*slopex(i  ,j,k,n)
!    $            + dth * force(i  ,j,k,n)
          smlft = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n)
!    $            + dth * force(i-1,j,k,n)

          smrgt = cvmgt(s(is-1,j,k,n),smrgt,i.eq.is .and. bc(1,1) .eq. INLET)
          smlft = cvmgt(s(is-1,j,k,n),smlft,i.eq.is .and. bc(1,1) .eq. INLET)

          if (i .eq. is  .and.  bc(1,1) .eq. WALL) then
            if (is_vel .and. n .eq. 1) then
              smrgt = ZERO
              smlft = ZERO
            else if (is_vel .and. n .ne. 1) then
              smrgt = cvmgt(ZERO,smrgt,visc_coef .gt. ZERO)
              smlft = cvmgt(ZERO,smrgt,visc_coef .gt. ZERO)
            else
              smlft = smrgt
            endif
          endif
 
          sminus = cvmgp(smlft,smrgt,utrans(i,j,k,1))
          savg   = HALF * (smlft + smrgt)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,1)) .gt. eps)

          str = HALF * (utrans(i,j,k,1)+utrans(i+1,j,k,1))*(splus - sminus) / hx

!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Y-DIRECTION
!        ******************************************************************

          spbot = s(i,j  ,k,n) + (HALF - dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
!    $            + dth * force(i,  j,k,n)
          sptop = s(i,j+1,k,n) - (HALF + dth*u(i,j+1,k,2)/hy)*slopey(i,j+1,k,n)
!    $            + dth * force(i,j+1,k,n)

          sptop = cvmgt(s(i,je+1,k,n),sptop,j.eq.je .and. bc(2,2) .eq. INLET)
          spbot = cvmgt(s(i,je+1,k,n),spbot,j.eq.je .and. bc(2,2) .eq. INLET)

          if (j .eq. je .and. bc(2,2) .eq. WALL) then
            if (is_vel .and. n .eq. 2) then
              sptop = ZERO
              spbot = ZERO
            else if (is_vel .and. n .ne. 2) then
              sptop = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
              spbot = cvmgt(ZERO,spbot,visc_coef .gt. ZERO)
            else
              sptop = spbot
            endif
          endif

          splus = cvmgp(spbot,sptop,utrans(i,j+1,k,2))
          savg  = HALF * (spbot + sptop)
          splus = cvmgt(splus, savg, abs(utrans(i,j+1,k,2)) .gt. eps)

          smtop = s(i,j  ,k,n) - (HALF + dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
!    $            + dth * force(i,j  ,k,n)
          smbot = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n)
!    $            + dth * force(i,j-1,k,n)

          smtop = cvmgt(s(i,js-1,k,n),smtop,j.eq.js .and. bc(2,1) .eq. INLET)
          smbot = cvmgt(s(i,js-1,k,n),smbot,j.eq.js .and. bc(2,1) .eq. INLET)

          if (j .eq. js  .and.  bc(2,1) .eq. WALL) then
            if (is_vel .and. n .eq. 2) then
              smtop = ZERO
              smbot = ZERO
            else if (is_vel .and. n .ne. 2) then
              smtop = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
              smbot = cvmgt(ZERO,smtop,visc_coef .gt. ZERO)
            else
              smbot = smtop
            endif
          endif

          sminus = cvmgp(smbot,smtop,utrans(i,j,k,2))
          savg   = HALF * (smbot + smtop)
          sminus = cvmgt(sminus, savg, abs(utrans(i,j,k,2)) .gt. eps)

          str =  str + HALF * (utrans(i,j,k,2)+utrans(i,j+1,k,2))*(splus - sminus) / hy

!        ******************************************************************
!        MAKE DOWN AND UP STATES
!        ******************************************************************

          st = force(i,j,k,n) - str

          wbardth = dth*u(i,j,k,3)/hz

          s_d(k+1)= s(i,j,k,n) + (HALF-wbardth)*slopez(i,j,k,n) + dth*st
          s_u(k  )= s(i,j,k,n) - (HALF+wbardth)*slopez(i,j,k,n) + dth*st

        enddo

        if (velpred .eq. 1) then
          do k = ks, ke+1 
            savg = HALF*(s_d(k) + s_u(k))
            test = ( (s_d(k) .le. ZERO  .and. &
                      s_u(k) .ge. ZERO)  .or. &
                   (abs(s_d(k) + s_u(k)) .lt. eps) )
            sedgez(i,j,k,n)=cvmgp(s_d(k),s_u(k),savg)
            sedgez(i,j,k,n)=cvmgt(savg,sedgez(i,j,k,n),test)
          enddo
        else
          do k = ks, ke+1 
            sedgez(i,j,k,n)=cvmgp(s_d(k),s_u(k),uadv(i,j,k,3))
            savg = HALF*(s_d(k) + s_u(k))
            sedgez(i,j,k,n)=cvmgt(savg,sedgez(i,j,k,n),abs(uadv(i,j,k,3)) .lt. eps)
          enddo
        endif

        if (bc(3,1) .eq. WALL) then
          if (is_vel .and. n .eq. 3) then
            sedgez(i,j,ks,n) = ZERO
          else if (is_vel .and. n .ne. 3) then
            sedgez(i,j,ks,n) = cvmgt(ZERO,s_u(ks),visc_coef .gt. 0.0)
          else
            sedgez(i,j,ks,n) = s_u(ks)
          endif
        elseif (bc(3,1) .eq. INLET) then
          sedgez(i,j,ks,n) = s(i,j,ks-1,n)
        elseif (bc(3,1) .eq. OUTLET) then
          sedgez(i,j,ks,n) = s_u(ks)
        endif

        if (bc(3,2) .eq. WALL) then
          if (is_vel .and. n .eq. 3) then
            sedgez(i,j,ke+1,n) = ZERO
          else if (is_vel .and. n .ne. 3) then
            sedgez(i,j,ke+1,n) = cvmgt(ZERO,s_d(ke+1),visc_coef .gt. 0.0)
          else
            sedgez(i,j,ke+1,n) = s_d(ke+1)
          endif
        elseif (bc(3,2) .eq. INLET) then
          sedgez(i,j,ke+1,n) = s(i,j,ke+1,n)
        elseif (bc(3,2) .eq. OUTLET) then
          sedgez(i,j,ke+1,n) = s_d(ke+1)
        endif

        if (velpred .eq. 1) then
          do k = ks, ke+1 
            uadv(i,j,k,3) = sedgez(i,j,k,3)
          enddo
        endif

        enddo
        enddo
       endif
      enddo

      end subroutine mkflux_3d

end module mkflux_module
