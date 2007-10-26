module velpred_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module

  implicit none

contains

      subroutine velpred_2d(s,u,umac,vmac,&
                            force,lo,dx,dt,&
                            phys_bc,adv_bc,ng)

      integer, intent(in) :: lo(2),ng

      real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:,:,:),slopey(:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:)

      real(kind=dp_t) hx, hy, dth
      real(kind=dp_t) splus,sminus
      real(kind=dp_t) savg,st
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt

      integer :: hi(2)
      integer :: n,ncomp
      integer :: slope_order = 4
      logical :: test

      real(kind=dp_t) :: abs_eps, eps, umax

      integer :: i,j,is,js,ie,je

      ncomp = size(s,dim=3)

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      allocate(s_l(lo(1)-1:hi(1)+2))
      allocate(s_r(lo(1)-1:hi(1)+2))
      allocate(s_b(lo(2)-1:hi(2)+2))
      allocate(s_t(lo(2)-1:hi(2)+2))

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))

      call slopex_2d(s,slopex,lo,ng,ncomp,adv_bc,slope_order)
      call slopey_2d(s,slopey,lo,ng,ncomp,adv_bc,slope_order)

      abs_eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)

      dth = HALF*dt

      hx = dx(1)
      hy = dx(2)
         
      umax = abs(u(is,js,1))
      do j = js,je
         do i = is,ie
            umax = max(umax,abs(u(i,j,1)))
            umax = max(umax,abs(u(i,j,2)))
         end do
      end do
      
      eps = abs_eps * umax
      
      !
      !     Loop to obtain velocities on x-edges.
      !
      do n = 1,ncomp

         if(n .eq. 1) then

            do j = js,je 
               do i = is-1,ie+1 
                  spbot = s(i,j  ,n) + (HALF - dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
                  sptop = s(i,j+1,n) - (HALF + dth*u(i,j+1,2)/hy) * slopey(i,j+1,n)
                  
                  sptop = merge(s(i,je+1,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  spbot = merge(s(i,je+1,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  
                  if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                     if (n .eq. 2) then
                        sptop = ZERO
                        spbot = ZERO
                     elseif (n .eq. 1) then
                        sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                        spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                     endif
                  endif
                  
                  savg = HALF * (spbot + sptop)
                  test = ((spbot .le. ZERO .and. sptop .ge. ZERO) .or. (abs(spbot+sptop) .lt. eps))
                  splus = merge(spbot,sptop,spbot+sptop.gt.ZERO)
                  splus = merge(savg,splus,test)
                  
                  smtop = s(i,j  ,n) - (HALF + dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
                  smbot = s(i,j-1,n) + (HALF - dth*u(i,j-1,2)/hy) * slopey(i,j-1,n)
                  
                  smtop = merge(s(i,js-1,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  smbot = merge(s(i,js-1,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  
                  if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                     if (n .eq. 2) then
                        smtop = ZERO
                        smbot = ZERO
                     elseif (n .ne. 2) then
                        smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                        smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                     endif
                  endif
                  
                  savg  = HALF * (smbot + smtop)
                  test = ((smbot .le. ZERO .and. smtop .ge. ZERO) .or. (abs(smbot+smtop) .lt. eps))
                  sminus = merge(smbot,smtop,smbot+smtop.gt.ZERO)
                  sminus = merge(savg,sminus,test)
                  
                  st = force(i,j,n) - &
                       HALF * (splus + sminus)*(splus - sminus) / hy
                  
                  s_l(i+1)= s(i,j,n) + (HALF-dth*u(i,j,1)/hx)*slopex(i,j,1) + dth*st
                  s_r(i  )= s(i,j,n) - (HALF+dth*u(i,j,1)/hx)*slopex(i,j,1) + dth*st
               enddo
               
               do i = is, ie+1 
                  savg = HALF*(s_r(i) + s_l(i))
                  test = ((s_l(i) .le. ZERO  .and. s_r(i) .ge. ZERO) .or. &
                       (abs(s_l(i) + s_r(i)) .lt. eps))
                  umac(i,j)=merge(s_l(i),s_r(i),savg.gt.ZERO)
                  umac(i,j)=merge(savg,umac(i,j),test)
               enddo
            
               if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                  if (n .eq. 1) then
                     umac(is,j) = ZERO
                  elseif (n .ne. 1) then
                     umac(is,j) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
                  endif
               elseif (phys_bc(1,1) .eq. INLET) then
                  umac(is,j) = s(is-1,j,n)
               elseif (phys_bc(1,1) .eq. OUTLET) then
                  umac(is,j) = s_r(is)
               endif
               if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                  if (n .eq. 1) then
                     umac(ie+1,j) = ZERO
                  else if (n .ne. 1) then
                     umac(ie+1,j) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
                  endif
               elseif (phys_bc(1,2) .eq. INLET) then
                  umac(ie+1,j) = s(ie+1,j,n)
               elseif (phys_bc(1,2) .eq. OUTLET) then
                  umac(ie+1,j) = s_l(ie+1)
               endif
            enddo
         endif
      enddo

      !
      !     Loop to obtain velocities on y-edges.
      !
      do n = 1,ncomp

         if(n .eq. 2) then

            do i = is, ie 
               do j = js-1, je+1 
                  splft = s(i,j  ,n) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
                  sprgt = s(i+1,j,n) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,n)
                  
                  sprgt = merge(s(ie+1,j,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  splft = merge(s(ie+1,j,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  
                  if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                     if (n .eq. 1) then
                        splft = ZERO
                        sprgt = ZERO
                     elseif (n .ne. 1) then
                        sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                        splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                     endif
                  endif
                  
                  savg = HALF * (splft + sprgt)
                  test = ((splft .le. ZERO .and. sprgt .ge. ZERO) .or. (abs(splft+sprgt) .lt. eps))
                  splus = merge(splft,sprgt,splft+sprgt.gt.ZERO)
                  splus = merge(savg,splus,test)
                  
                  smrgt = s(i  ,j,n) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
                  smlft = s(i-1,j,n) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,n)
                  
                  smrgt = merge(s(is-1,j,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  smlft = merge(s(is-1,j,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  
                  if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                     if (n .eq. 1) then
                        smlft = ZERO
                        smrgt = ZERO
                     elseif (n .ne. 1) then
                        smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                        smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                     endif
                  endif
                  
                  savg  = HALF * (smlft + smrgt)
                  test = ((smlft .le. ZERO .and. smrgt .ge. ZERO) .or. (abs(smlft+smrgt) .lt. eps))
                  sminus = merge(smlft,smrgt,smlft+smrgt.gt.ZERO)
                  sminus = merge(savg,sminus,test)
                  
                  st = force(i,j,n) - &
                       HALF * (splus + sminus)*(splus - sminus) / hx
                  
                  s_b(j+1)= s(i,j,n) + (HALF-dth*u(i,j,2)/hy)*slopey(i,j,1) + dth*st
                  s_t(j  )= s(i,j,n) - (HALF+dth*u(i,j,2)/hy)*slopey(i,j,1) + dth*st
               enddo
               
               do j = js, je+1 
                  savg = HALF*(s_b(j) + s_t(j))
                  test = ((s_b(j) .le. ZERO  .and. s_t(j) .ge. ZERO) .or. &
                       (abs(s_b(j) + s_t(j)) .lt. eps))
                  vmac(i,j)=merge(s_b(j),s_t(j),savg.gt.ZERO)
                  vmac(i,j)=merge(savg,vmac(i,j),test)
               enddo
               
               if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                  if (n .eq. 2) then
                     vmac(i,js) = ZERO
                  elseif (n .ne. 2) then
                     vmac(i,js) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
                  endif
               elseif (phys_bc(2,1) .eq. INLET) then
                  vmac(i,js) = s(i,js-1,n)
               elseif (phys_bc(2,1) .eq. OUTLET) then
                  vmac(i,js) = s_t(js)
               endif
               
               if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                  if (n .eq. 2) then
                     vmac(i,je+1) = ZERO
                  elseif (n .ne. 2) then
                     vmac(i,je+1) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
                  endif
               elseif (phys_bc(2,2) .eq. INLET) then
                  vmac(i,je+1) = s(i,je+1,n)
               elseif (phys_bc(2,2) .eq. OUTLET) then
                  vmac(i,je+1) = s_b(je+1)
               endif
            enddo
         endif
      enddo
      
      deallocate(s_l)
      deallocate(s_r)
      deallocate(s_b)
      deallocate(s_t)

      deallocate(slopex)
      deallocate(slopey)

      end subroutine velpred_2d

      subroutine velpred_3d(s,u,&
                            umac,vmac,wmac, &
                            force,lo,dx,dt, &
                            phys_bc,adv_bc,ng)

      integer, intent(in) :: lo(:),ng

      real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:,:,:,:)
      real(kind=dp_t), allocatable::  slopey(:,:,:,:)
      real(kind=dp_t), allocatable::  slopez(:,:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)

      real(kind=dp_t) hx, hy, hz, dth, dt4, dt6
      real(kind=dp_t) splus,sminus,st,str,savg
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
      logical test
      logical cornerCoupling

      real(kind=dp_t) :: abs_eps,eps,umax

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke,n
      integer :: slope_order = 4
      integer :: ncomp

      ! these correspond to s_L^x, etc.
      real(kind=dp_t), allocatable:: slx(:,:,:),srx(:,:,:)
      real(kind=dp_t), allocatable:: sly(:,:,:),sry(:,:,:)
      real(kind=dp_t), allocatable:: slz(:,:,:),srz(:,:,:)

      ! these correspond to s_{\i-\half\e_x}^x, etc.
      real(kind=dp_t), allocatable:: simhx(:,:,:),simhy(:,:,:),simhz(:,:,:)

      ! these correspond to s_L^{x|y}, etc.
      real(kind=dp_t), allocatable:: slxy(:,:,:),srxy(:,:,:),slxz(:,:,:),srxz(:,:,:)
      real(kind=dp_t), allocatable:: slyx(:,:,:),sryx(:,:,:),slyz(:,:,:),sryz(:,:,:)
      real(kind=dp_t), allocatable:: slzx(:,:,:),srzx(:,:,:),slzy(:,:,:),srzy(:,:,:)

      ! these correspond to s_{\i-\half\e_x}^{x|y}, etc.
      real(kind=dp_t), allocatable:: simhxy(:,:,:),simhxz(:,:,:)
      real(kind=dp_t), allocatable:: simhyx(:,:,:),simhyz(:,:,:)
      real(kind=dp_t), allocatable:: simhzx(:,:,:),simhzy(:,:,:)

      ! these correspond to \mathrm{sedge}_L^x, etc.
      real(kind=dp_t), allocatable:: sedgelx(:,:,:),sedgerx(:,:,:)
      real(kind=dp_t), allocatable:: sedgely(:,:,:),sedgery(:,:,:)
      real(kind=dp_t), allocatable:: sedgelz(:,:,:),sedgerz(:,:,:)

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

      ! these are allocated from lo:hi+1 in the normal direction
      ! and from lo-1:hi+1 in the transverse directions
      allocate(slx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(srx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(sly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(sry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(slz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(srz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

      ! these are allocated from lo:hi+1 in the normal direction
      ! and from lo-1:hi+1 in the transverse directions
      allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))


      ! The general pattern is
      ! lo:hi+1 in normal direction
      ! lo:hi in transverse direction
      ! lo-1:hi+1 in unused direction
      allocate(slxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
      allocate(srxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
      allocate(simhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

      allocate(slxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
      allocate(srxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
      allocate(simhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

      allocate(slyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(sryx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

      allocate(slyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(sryz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(simhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

      allocate(slzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(srzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(simhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

      allocate(slzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
      allocate(srzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
      allocate(simhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

      ! these are allocated from lo:hi+1 in the normal direction
      ! and from lo:hi in the transverse directions
      allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
      allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
      allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(sedgelz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
      allocate(sedgerz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

      ncomp = size(s,dim=4)
      do k = lo(3)-1,hi(3)+1
         call slopex_2d(s(:,:,k,:),slopex(:,:,k,:),lo,ng,ncomp,adv_bc,slope_order)
         call slopey_2d(s(:,:,k,:),slopey(:,:,k,:),lo,ng,ncomp,adv_bc,slope_order)
      end do
      call slopez_3d(s,slopez,lo,ng,ncomp,adv_bc,slope_order)

      abs_eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

      dth = HALF*dt
      dt4 = dt/4.0d0
      dt6 = dt/6.0d0

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      umax = abs(u(is,js,ks,1))
      do k = ks,ke
         do j = js,je
            do i = is,ie
               umax = max(umax,abs(u(i,j,k,1)))
               umax = max(umax,abs(u(i,j,k,2)))
               umax = max(umax,abs(u(i,j,k,3)))
            end do
         end do
      end do
           
      eps = abs_eps * umax
      
      ! loop over components
      do n = 1,ncomp

!******************************************************************
! Create s_{\i-\half\e_x}^x, etc.
!******************************************************************
         
         ! loop over appropriate x-faces
         do k=ks-1,ke+1
            do j=js-1,je+1
               do i=is,ie+1
                  ! make slx, srx with 1D extrapolation
                  slx(i,j,k) = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n)
                  srx(i,j,k) = s(i  ,j,k,n) - (HALF + dth*u(i,  j,k,1)/hx)*slopex(i,  j,k,n)
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slx(i,j,k) = merge(s(is-1,j,k,n),slx(i,j,k),phys_bc(1,1) .eq. INLET)
                     srx(i,j,k) = merge(s(is-1,j,k,n),srx(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 1) then
                           slx(i,j,k) = ZERO
                           srx(i,j,k) = ZERO
                        else if(n .ne. 1) then
                           slx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slx(i,j,k) = merge(s(ie+1,j,k,n),slx(i,j,k),phys_bc(1,2) .eq. INLET)
                     srx(i,j,k) = merge(s(ie+1,j,k,n),srx(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 1) then
                           slx(i,j,k) = ZERO
                           srx(i,j,k) = ZERO
                        else if (n .ne. 1) then
                           slx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhx by solving Riemann problem
                  simhx(i,j,k) = merge(slx(i,j,k),srx(i,j,k),slx(i,j,k)+srx(i,j,k) .gt. ZERO)
                  savg = HALF*(slx(i,j,k)+srx(i,j,k))
                  simhx(i,j,k) = merge(simhx(i,j,k),savg,abs(slx(i,j,k)+srx(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate y-faces
         do k=ks-1,ke+1
            do j=js,je+1
               do i=is-1,ie+1
                  ! make sly, sry with 1D extrapolation
                  sly(i,j,k) = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n)
                  sry(i,j,k) = s(i,j,  k,n) - (HALF + dth*u(i,j,  k,2)/hy)*slopey(i,j,  k,n)
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     sly(i,j,k) = merge(s(is,j-1,k,n),sly(i,j,k),phys_bc(2,1) .eq. INLET)
                     sry(i,j,k) = merge(s(is,j-1,k,n),sry(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 2) then
                           sly(i,j,k) = ZERO
                           sry(i,j,k) = ZERO
                        else if(n .ne. 2) then
                           sly(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sry(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     sly(i,j,k) = merge(s(i,je+1,k,n),sly(i,j,k),phys_bc(2,2) .eq. INLET)
                     sry(i,j,k) = merge(s(i,je+1,k,n),sry(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 2) then
                           sly(i,j,k) = ZERO
                           sry(i,j,k) = ZERO
                        else if (n .ne. 2) then
                           sly(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sry(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhy by solving Riemann problem
                  simhy(i,j,k) = merge(sly(i,j,k),sry(i,j,k),sly(i,j,k)+sry(i,j,k) .gt. ZERO)
                  savg = HALF*(sly(i,j,k)+sry(i,j,k))
                  simhy(i,j,k) = merge(simhy(i,j,k),savg,abs(sly(i,j,k)+sry(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate z-faces
         do k=ks,ke+1
            do j=js-1,je+1
               do i=is-1,ie+1
                  ! make slz, srz with 1D extrapolation
                  slz(i,j,k) = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n)
                  srz(i,j,k) = s(i,j,k,  n) - (HALF + dth*u(i,j,k,  3)/hz)*slopez(i,j,k,  n)
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slz(i,j,k) = merge(s(is,j,k-1,n),slz(i,j,k),phys_bc(3,1) .eq. INLET)
                     srz(i,j,k) = merge(s(is,j,k-1,n),srz(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 3) then
                           slz(i,j,k) = ZERO
                           srz(i,j,k) = ZERO
                        else if(n .ne. 3) then
                           slz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slz(i,j,k) = merge(s(i,j,ke+1,n),slz(i,j,k),phys_bc(3,2) .eq. INLET)
                     srz(i,j,k) = merge(s(i,j,ke+1,n),srz(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 3) then
                           slz(i,j,k) = ZERO
                           srz(i,j,k) = ZERO
                        else if (n .ne. 3) then
                           slz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhz by solving Riemann problem
                  simhz(i,j,k) = merge(slz(i,j,k),srz(i,j,k),slz(i,j,k)+srz(i,j,k) .gt. ZERO)
                  savg = HALF*(slz(i,j,k)+srz(i,j,k))
                  simhz(i,j,k) = merge(simhz(i,j,k),savg,abs(slz(i,j,k)+srz(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo

!******************************************************************
! Create s_{\i-\half\e_x}^{x|y}, etc.
!******************************************************************

         ! loop over appropriate xy faces
         do k=ks-1,ke+1
            do j=js,je
               do i=is,ie+1
                  ! make slxy, srxy by updating 1D extrapolation
                  slxy(i,j,k) = slx(i,j,k) - (dt6/hy)*(simhy(i-1,j+1,k)+simhy(i-1,j,k))*(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                  srxy(i,j,k) = srx(i,j,k) - (dt6/hy)*(simhy(i,  j+1,k)+simhy(i,  j,k))*(simhy(i,  j+1,k)-simhy(i,  j,k))
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slxy(i,j,k) = merge(s(is-1,j,k,n),slxy(i,j,k),phys_bc(1,1) .eq. INLET)
                     srxy(i,j,k) = merge(s(is-1,j,k,n),srxy(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 1) then
                           slxy(i,j,k) = ZERO
                           srxy(i,j,k) = ZERO
                        else if(n .ne. 1) then
                           slxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slxy(i,j,k) = merge(s(ie+1,j,k,n),slxy(i,j,k),phys_bc(1,2) .eq. INLET)
                     srxy(i,j,k) = merge(s(ie+1,j,k,n),srxy(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 1) then
                           slxy(i,j,k) = ZERO
                           srxy(i,j,k) = ZERO
                        else if (n .ne. 1) then
                           slxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhxy by solving Riemann problem
                  simhxy(i,j,k) = merge(slxy(i,j,k),srxy(i,j,k),slxy(i,j,k)+srxy(i,j,k) .gt. ZERO)
                  savg = HALF*(slxy(i,j,k)+srxy(i,j,k))
                  simhxy(i,j,k) = merge(simhxy(i,j,k),savg,abs(slxy(i,j,k)+srxy(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate xz faces
         do k=ks,ke
            do j=js-1,je+1
               do i=is,ie+1
                  ! make slxz, srxz by updating 1D extrapolation
                  slxz(i,j,k) = slx(i,j,k) - (dt6/hz)*(simhz(i-1,j,k+1)+simhz(i-1,j,k))*(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                  srxz(i,j,k) = srx(i,j,k) - (dt6/hz)*(simhz(i,  j,k+1)+simhz(i,  j,k))*(simhz(i,  j,k+1)-simhz(i,  j,k))
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slxz(i,j,k) = merge(s(is-1,j,k,n),slxz(i,j,k),phys_bc(1,1) .eq. INLET)
                     srxz(i,j,k) = merge(s(is-1,j,k,n),srxz(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 1) then
                           slxz(i,j,k) = ZERO
                           srxz(i,j,k) = ZERO
                        else if(n .ne. 1) then
                           slxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slxz(i,j,k) = merge(s(ie+1,j,k,n),slxz(i,j,k),phys_bc(1,2) .eq. INLET)
                     srxz(i,j,k) = merge(s(ie+1,j,k,n),srxz(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 1) then
                           slxz(i,j,k) = ZERO
                           srxz(i,j,k) = ZERO
                        else if (n .ne. 1) then
                           slxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhxz by solving Riemann problem
                  simhxz(i,j,k) = merge(slxz(i,j,k),srxz(i,j,k),slxz(i,j,k)+srxz(i,j,k) .gt. ZERO)
                  savg = HALF*(slxz(i,j,k)+srxz(i,j,k))
                  simhxz(i,j,k) = merge(simhxz(i,j,k),savg,abs(slxz(i,j,k)+srxz(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo

         ! loop over appropriate yx faces
         do k=ks-1,ke+1
            do j=js,je+1
               do i=is,ie
                  ! make slyx, sryx by updating 1D extrapolation
                  slyx(i,j,k) = sly(i,j,k) - (dt6/hx)*(simhx(i+1,j-1,k)+simhx(i,j-1,k))*(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                  sryx(i,j,k) = sry(i,j,k) - (dt6/hx)*(simhx(i+1,j,  k)+simhx(i,j,  k))*(simhx(i+1,j,  k)-simhx(i,j,  k))
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     slyx(i,j,k) = merge(s(is,j-1,k,n),slyx(i,j,k),phys_bc(2,1) .eq. INLET)
                     sryx(i,j,k) = merge(s(is,j-1,k,n),sryx(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 2) then
                           slyx(i,j,k) = ZERO
                           sryx(i,j,k) = ZERO
                        else if(n .ne. 2) then
                           slyx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sryx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     slyx(i,j,k) = merge(s(i,je+1,k,n),slyx(i,j,k),phys_bc(2,2) .eq. INLET)
                     sryx(i,j,k) = merge(s(i,je+1,k,n),sryx(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 2) then
                           slyx(i,j,k) = ZERO
                           sryx(i,j,k) = ZERO
                        else if (n .ne. 2) then
                           slyx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sryx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhyx by solving Riemann problem
                  simhyx(i,j,k) = merge(slyx(i,j,k),sryx(i,j,k),slyx(i,j,k)+sryx(i,j,k) .gt. ZERO)
                  savg = HALF*(slyx(i,j,k)+sryx(i,j,k))
                  simhyx(i,j,k) = merge(simhyx(i,j,k),savg,abs(slyx(i,j,k)+sryx(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate yz faces
         do k=ks,ke
            do j=js,je+1
               do i=is-1,ie+1
                  ! make slyz, sryz by updating 1D extrapolation
                  slyz(i,j,k) = sly(i,j,k) - (dt6/hz)*(simhz(i,j-1,k+1)+simhz(i,j-1,k))*(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                  sryz(i,j,k) = sry(i,j,k) - (dt6/hz)*(simhz(i,j,  k+1)+simhz(i,j,  k))*(simhz(i,j,  k+1)-simhz(i,j,  k))
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     slyz(i,j,k) = merge(s(is,j-1,k,n),slyz(i,j,k),phys_bc(2,1) .eq. INLET)
                     sryz(i,j,k) = merge(s(is,j-1,k,n),sryz(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 2) then
                           slyz(i,j,k) = ZERO
                           sryz(i,j,k) = ZERO
                        else if(n .ne. 2) then
                           slyz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sryz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     slyz(i,j,k) = merge(s(i,je+1,k,n),slyz(i,j,k),phys_bc(2,2) .eq. INLET)
                     sryz(i,j,k) = merge(s(i,je+1,k,n),sryz(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 2) then
                           slyz(i,j,k) = ZERO
                           sryz(i,j,k) = ZERO
                        else if (n .ne. 2) then
                           slyz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sryz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhyz by solving Riemann problem
                  simhyz(i,j,k) = merge(slyz(i,j,k),sryz(i,j,k),slyz(i,j,k)+sryz(i,j,k) .gt. ZERO)
                  savg = HALF*(slyz(i,j,k)+sryz(i,j,k))
                  simhyz(i,j,k) = merge(simhyz(i,j,k),savg,abs(slyz(i,j,k)+sryz(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate zx faces
         do k=ks,ke+1
            do j=js-1,je+1
               do i=is,ie
                  ! make slzx, srzx by updating 1D extrapolation
                  slzx(i,j,k) = slz(i,j,k) - (dt6/hx)*(simhx(i+1,j,k-1)+simhx(i,j,k-1))*(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                  srzx(i,j,k) = srz(i,j,k) - (dt6/hx)*(simhx(i+1,j,k  )+simhx(i,j,k  ))*(simhx(i+1,j,k  )-simhx(i,j,k  ))
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slzx(i,j,k) = merge(s(is,j,k-1,n),slzx(i,j,k),phys_bc(3,1) .eq. INLET)
                     srzx(i,j,k) = merge(s(is,j,k-1,n),srzx(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 3) then
                           slzx(i,j,k) = ZERO
                           srzx(i,j,k) = ZERO
                        else if(n .ne. 3) then
                           slzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slzx(i,j,k) = merge(s(i,j,ke+1,n),slzx(i,j,k),phys_bc(3,2) .eq. INLET)
                     srzx(i,j,k) = merge(s(i,j,ke+1,n),srzx(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 3) then
                           slzx(i,j,k) = ZERO
                           srzx(i,j,k) = ZERO
                        else if (n .ne. 3) then
                           slzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhzx by solving Riemann problem
                  simhzx(i,j,k) = merge(slzx(i,j,k),srzx(i,j,k),slzx(i,j,k)+srzx(i,j,k) .gt. ZERO)
                  savg = HALF*(slzx(i,j,k)+srzx(i,j,k))
                  simhzx(i,j,k) = merge(simhzx(i,j,k),savg,abs(slzx(i,j,k)+srzx(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate zy faces
         do k=ks,ke+1
            do j=js,je
               do i=is-1,ie+1
                  ! make slzy, srzy by updating 1D extrapolation
                  slzy(i,j,k) = slz(i,j,k) - (dt6/hy)*(simhy(i,j+1,k-1)+simhy(i,j,k-1))*(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                  srzy(i,j,k) = srz(i,j,k) - (dt6/hy)*(simhy(i,j+1,k  )+simhy(i,j,k  ))*(simhy(i,j+1,k  )-simhy(i,j,k  ))
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slzy(i,j,k) = merge(s(is,j,k-1,n),slzy(i,j,k),phys_bc(3,1) .eq. INLET)
                     srzy(i,j,k) = merge(s(is,j,k-1,n),srzy(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 3) then
                           slzy(i,j,k) = ZERO
                           srzy(i,j,k) = ZERO
                        else if(n .ne. 3) then
                           slzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slzy(i,j,k) = merge(s(i,j,ke+1,n),slzy(i,j,k),phys_bc(3,2) .eq. INLET)
                     srzy(i,j,k) = merge(s(i,j,ke+1,n),srzy(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 3) then
                           slzy(i,j,k) = ZERO
                           srzy(i,j,k) = ZERO
                        else if (n .ne. 3) then
                           slzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make simhzy by solving Riemann problem
                  simhzy(i,j,k) = merge(slzy(i,j,k),srzy(i,j,k),slzy(i,j,k)+srzy(i,j,k) .gt. ZERO)
                  savg = HALF*(slzy(i,j,k)+srzy(i,j,k))
                  simhzy(i,j,k) = merge(simhzy(i,j,k),savg,abs(slzy(i,j,k)+srzy(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
!******************************************************************
! Create sedgelx, etc.
!******************************************************************

         ! loop over appropriate x-faces
         do k=ks,ke
            do j=js,je
               do i=is,ie+1
                  ! make sedgelx, sedgerx
                  sedgelx(i,j,k) = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n) &
                       - (dt4/hy)*(simhyz(i-1,j+1,k)+simhyz(i-1,j,k))*(simhyz(i-1,j+1,k)-simhyz(i-1,j,k)) &
                       - (dt4/hz)*(simhzy(i-1,j,k+1)+simhzy(i-1,j,k))*(simhzy(i-1,j,k+1)-simhzy(i-1,j,k)) &
                       + dth*force(i,j,k,n)
                  sedgerx(i,j,k) = s(i,  j,k,n) - (HALF + dth*u(i,  j,k,1)/hx)*slopex(i,  j,k,n) &
                       - (dt4/hy)*(simhyz(i,  j+1,k)+simhyz(i,  j,k))*(simhyz(i,  j+1,k)-simhyz(i,  j,k)) &
                       - (dt4/hz)*(simhzy(i,  j,k+1)+simhzy(i,  j,k))*(simhzy(i,  j,k+1)-simhzy(i,  j,k)) &
                       + dth*force(i,j,k,n)
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     sedgelx(i,j,k) = merge(s(is-1,j,k,n),sedgelx(i,j,k),phys_bc(1,1) .eq. INLET)
                     sedgerx(i,j,k) = merge(s(is-1,j,k,n),sedgerx(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 1) then
                           sedgelx(i,j,k) = ZERO
                           sedgerx(i,j,k) = ZERO
                        else if(n .ne. 1) then
                           sedgelx(i,j,k) = merge(ZERO,sedgerx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           sedgerx(i,j,k) = merge(ZERO,sedgerx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     sedgelx(i,j,k) = merge(s(ie+1,j,k,n),sedgelx(i,j,k),phys_bc(1,2) .eq. INLET)
                     sedgerx(i,j,k) = merge(s(ie+1,j,k,n),sedgerx(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 1) then
                           sedgelx(i,j,k) = ZERO
                           sedgerx(i,j,k) = ZERO
                        else if (n .ne. 1) then
                           sedgelx(i,j,k) = merge(ZERO,sedgelx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           sedgerx(i,j,k) = merge(ZERO,sedgelx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make umac by solving Riemann problem
                  savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
                  test = ((sedgelx(i,j,k) .le. ZERO .and. sedgerx(i,j,k) .ge. ZERO) .or. &
                       (abs(sedgelx(i,j,k)+sedgerx(i,j,k)) .lt. eps))
                  umac(i,j,k) = merge(sedgelx(i,j,k),sedgerx(i,j,k),savg .gt. ZERO)
                  umac(i,j,k) = merge(0.d0,umac(i,j,k),test)
               enddo
            enddo
         enddo
         
         ! loop over appropriate y-faces
         do k=ks,ke
            do j=js,je+1
               do i=is,ie
                  ! make sedgely, sedgery
                  sedgely(i,j,k) = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n) &
                       - (dt4/hx)*(simhxz(i+1,j-1,k)+simhxz(i,j-1,k))*(simhxz(i+1,j-1,k)-simhxz(i,j-1,k)) &
                       - (dt4/hz)*(simhzx(i,j-1,k+1)+simhzx(i,j-1,k))*(simhzx(i,j-1,k+1)-simhzx(i,j-1,k)) &
                       + dth*force(i,j,k,n)
                  sedgery(i,j,k) = s(i,j,  k,n) - (HALF + dth*u(i,j,  k,2)/hy)*slopey(i,j,  k,n) &
                       - (dt4/hx)*(simhxz(i+1,j,  k)+simhxz(i,j,  k))*(simhxz(i+1,j,  k)-simhxz(i,j,  k)) &
                       - (dt4/hz)*(simhzx(i,j,  k+1)+simhzx(i,j,  k))*(simhzx(i,j,  k+1)-simhzx(i,j,  k)) &
                       + dth*force(i,j,k,n)
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     sedgely(i,j,k) = merge(s(is,j-1,k,n),sedgely(i,j,k),phys_bc(2,1) .eq. INLET)
                     sedgery(i,j,k) = merge(s(is,j-1,k,n),sedgery(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 2) then
                           sedgely(i,j,k) = ZERO
                           sedgery(i,j,k) = ZERO
                        else if(n .ne. 2) then
                           sedgely(i,j,k) = merge(ZERO,sedgery(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sedgery(i,j,k) = merge(ZERO,sedgery(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     sedgely(i,j,k) = merge(s(i,je+1,k,n),sedgely(i,j,k),phys_bc(2,2) .eq. INLET)
                     sedgery(i,j,k) = merge(s(i,je+1,k,n),sedgery(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 2) then
                           sedgely(i,j,k) = ZERO
                           sedgery(i,j,k) = ZERO
                        else if (n .ne. 2) then
                           sedgely(i,j,k) = merge(ZERO,sedgely(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sedgery(i,j,k) = merge(ZERO,sedgely(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make vmac by solving Riemann problem
                  savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
                  test = ((sedgely(i,j,k) .le. ZERO .and. sedgery(i,j,k) .ge. ZERO) .or. &
                       (abs(sedgely(i,j,k)+sedgery(i,j,k)) .lt. eps))
                  vmac(i,j,k) = merge(sedgely(i,j,k),sedgery(i,j,k),savg .gt. ZERO)
                  vmac(i,j,k) = merge(0.d0,vmac(i,j,k),test)
               enddo
            enddo
         enddo
         
         ! loop over appropriate z-faces
         do k=ks,ke+1
            do j=js,je
               do i=is,ie
                  ! make sedgelz, sedgerz
                  sedgelz(i,j,k) = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n) &
                       - (dt4/hx)*(simhxy(i+1,j,k-1)+simhxy(i,j,k-1))*(simhxy(i+1,j,k-1)-simhxy(i,j,k-1)) &
                       - (dt4/hy)*(simhyx(i,j+1,k-1)+simhyx(i,j,k-1))*(simhyx(i,j+1,k-1)-simhyx(i,j,k-1)) &
                       + dth*force(i,j,k,n)
                  sedgerz(i,j,k) = s(i,j,k,  n) - (HALF + dth*u(i,j,k,  3)/hz)*slopez(i,j,k,  n) &
                       - (dt4/hx)*(simhxy(i+1,j,k  )+simhxy(i,j,k  ))*(simhxy(i+1,j,k  )-simhxy(i,j,k  )) &
                       - (dt4/hy)*(simhyx(i,j+1,k  )+simhyx(i,j,k  ))*(simhyx(i,j+1,k  )-simhyx(i,j,k  )) &
                       + dth*force(i,j,k,n)
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     sedgelz(i,j,k) = merge(s(is,j,k-1,n),sedgelz(i,j,k),phys_bc(3,1) .eq. INLET)
                     sedgerz(i,j,k) = merge(s(is,j,k-1,n),sedgerz(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(n .eq. 3) then
                           sedgelz(i,j,k) = ZERO
                           sedgerz(i,j,k) = ZERO
                        else if(n .ne. 3) then
                           sedgelz(i,j,k) = merge(ZERO,sedgerz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           sedgerz(i,j,k) = merge(ZERO,sedgerz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     sedgelz(i,j,k) = merge(s(i,j,ke+1,n),sedgelz(i,j,k),phys_bc(3,2) .eq. INLET)
                     sedgerz(i,j,k) = merge(s(i,j,ke+1,n),sedgerz(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (n .eq. 3) then
                           sedgelz(i,j,k) = ZERO
                           sedgerz(i,j,k) = ZERO
                        else if (n .ne. 3) then
                           sedgelz(i,j,k) = merge(ZERO,sedgelz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           sedgerz(i,j,k) = merge(ZERO,sedgelz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        endif
                     endif
                  endif
                  
                  ! make wmac by solving Riemann problem
                  savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
                  test = ((sedgelz(i,j,k) .le. ZERO .and. sedgerz(i,j,k) .ge. ZERO) .or. &
                       (abs(sedgelz(i,j,k)+sedgerz(i,j,k)) .lt. eps))
                  wmac(i,j,k) = merge(sedgelz(i,j,k),sedgerz(i,j,k),savg .gt. ZERO)
                  wmac(i,j,k) = merge(0.d0,wmac(i,j,k),test)
               enddo
            enddo
         enddo
      enddo
      
      deallocate(s_l)
      deallocate(s_r)
      deallocate(s_b)
      deallocate(s_t)
      deallocate(s_d)
      deallocate(s_u)

      deallocate(slopex)
      deallocate(slopey)
      deallocate(slopez)

      deallocate(slx)
      deallocate(srx)
      deallocate(sly)
      deallocate(sry)
      deallocate(slz)
      deallocate(srz)

      deallocate(simhx)
      deallocate(simhy)
      deallocate(simhz)

      deallocate(slxy)
      deallocate(srxy)
      deallocate(slxz)
      deallocate(srxz)
      deallocate(slyx)
      deallocate(sryx)
      deallocate(slyz)
      deallocate(sryz)
      deallocate(slzx)
      deallocate(srzx)
      deallocate(slzy)
      deallocate(srzy)

      deallocate(simhxy)
      deallocate(simhxz)
      deallocate(simhyx)
      deallocate(simhyz)
      deallocate(simhzx)
      deallocate(simhzy)

      deallocate(sedgelx)
      deallocate(sedgerx)
      deallocate(sedgely)
      deallocate(sedgery)
      deallocate(sedgelz)
      deallocate(sedgerz)

      end subroutine velpred_3d

end module velpred_module
