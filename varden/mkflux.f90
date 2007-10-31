module mkflux_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module

  implicit none

contains

      subroutine mkflux_2d(s,u,sedgex,sedgey,umac,vmac, &
                           force,lo,dx,dt,is_vel, &
                           phys_bc,adv_bc,ng)

      integer, intent(in) :: lo(2),ng

      real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) :: adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:,:,:),slopey(:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:)

      real(kind=dp_t) hx, hy, dt2
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

      dt2 = HALF*dt

      hx = dx(1)
      hy = dx(2)

      ! Compute eps, which is relative to the max mac velocity
      umax = abs(umac(is,js))
      do j = js,je
         do i = is,ie+1
            umax = max(umax,abs(umac(i,j)))
         end do
      end do
      do j = js,je+1
         do i = is,ie
            umax = max(umax,abs(vmac(i,j)))
         end do
      end do
      eps = abs_eps * umax
      
      !
      !     Loop for fluxes on x-edges.
      !
      do n = 1,ncomp
         do j = js,je 
            do i = is-1,ie+1 
               ! extrapolate to transverse faces on hi side of cell
               spbot = s(i,j  ,n) + (HALF - dt2*u(i,j  ,2)/hy) * slopey(i,j  ,n)
               sptop = s(i,j+1,n) - (HALF + dt2*u(i,j+1,2)/hy) * slopey(i,j+1,n)

               ! apply boundary conditions on hi-y side
               sptop = merge(s(i,je+1,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
               spbot = merge(s(i,je+1,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
               if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                  if (is_vel .and. n .eq. 2) then
                     sptop = ZERO
                     spbot = ZERO
                  elseif (is_vel .and. n .eq. 1) then
                     sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                     spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                  else
                     sptop = spbot
                  endif
               endif
               
               ! upwind based on vmac
               splus = merge(spbot,sptop,vmac(i,j+1).gt.ZERO)
               savg  = HALF * (spbot + sptop)
               splus = merge(splus, savg, abs(vmac(i,j+1)) .gt. eps)
               
               ! extrapolate to transverse faces on lo side of cell
               smtop = s(i,j  ,n) - (HALF + dt2*u(i,j  ,2)/hy) * slopey(i,j  ,n)
               smbot = s(i,j-1,n) + (HALF - dt2*u(i,j-1,2)/hy) * slopey(i,j-1,n)
               
               ! apply boundary conditions on lo-y side
               smtop = merge(s(i,js-1,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
               smbot = merge(s(i,js-1,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
               if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                  if (is_vel .and. (n .eq. 2)) then
                     smtop = ZERO
                     smbot = ZERO
                  elseif (is_vel .and. (n .ne. 2)) then
                     smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                     smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                  else
                     smbot = smtop
                  endif
               endif
               
               ! upwind based on vmac
               sminus = merge(smbot,smtop,vmac(i,j).gt.ZERO)
               savg   = HALF * (smbot + smtop)
               sminus = merge(sminus, savg, abs(vmac(i,j)) .gt. eps)
               
               ! create transverse plus force term
               st = force(i,j,n) - &
                    HALF * (vmac(i,j)+vmac(i,j+1))*(splus - sminus) / hy
               
               ! extrapolate to normal faces with transverse and force contribution
               s_l(i+1)= s(i,j,n) + (HALF-dt2*umac(i+1,j)/hx)*slopex(i,j,n) + dt2*st
               s_r(i  )= s(i,j,n) - (HALF+dt2*umac(i,j)/hx)*slopex(i,j,n) + dt2*st
            enddo
            
            ! upwind based on umac
            do i = is, ie+1 
               sedgex(i,j,n)=merge(s_l(i),s_r(i),umac(i,j).gt.ZERO)
               savg = HALF*(s_r(i) + s_l(i))
               sedgex(i,j,n)=merge(savg,sedgex(i,j,n),abs(umac(i,j)) .lt. eps)
            enddo
            
            ! apply boundary conditions on lo-x side of domain
            if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
               if (is_vel .and. n .eq. 1) then
                  sedgex(is,j,n) = ZERO
               elseif (is_vel .and. n .ne. 1) then
                  sedgex(is,j,n) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
               else 
                  sedgex(is,j,n) = s_r(is)
               endif
            elseif (phys_bc(1,1) .eq. INLET) then
               sedgex(is,j,n) = s(is-1,j,n)
            elseif (phys_bc(1,1) .eq. OUTLET) then
               sedgex(is,j,n) = s_r(is)
            endif

            ! apply boundary conditions on hi-x side of domain
            if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
               if (is_vel .and. n .eq. 1) then
                  sedgex(ie+1,j,n) = ZERO
               else if (is_vel .and. n .ne. 1) then
                  sedgex(ie+1,j,n) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
               else 
                  sedgex(ie+1,j,n) = s_l(ie+1)
               endif
            elseif (phys_bc(1,2) .eq. INLET) then
               sedgex(ie+1,j,n) = s(ie+1,j,n)
            elseif (phys_bc(1,2) .eq. OUTLET) then
               sedgex(ie+1,j,n) = s_l(ie+1)
            endif
         enddo
      enddo

      !
      !     Loop for fluxes on y-edges.
      !
      do n = 1,ncomp
         do i = is, ie 
            do j = js-1, je+1 
               ! extrapolate to transverse faces on hi side of cell
               splft = s(i,j  ,n) + (HALF - dt2*u(i  ,j,1)/hx) * slopex(i  ,j,n)
               sprgt = s(i+1,j,n) - (HALF + dt2*u(i+1,j,1)/hx) * slopex(i+1,j,n)
               
               ! apply boundary conditions on hi-x side
               sprgt = merge(s(ie+1,j,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
               splft = merge(s(ie+1,j,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
               if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                  if (is_vel .and. n .eq. 1) then
                     splft = ZERO
                     sprgt = ZERO
                  elseif (is_vel .and. n .ne. 1) then
                     sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                     splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                  else
                     sprgt = splft
                  endif
               endif
               
               ! upwind based on umac
               splus = merge(splft,sprgt,umac(i+1,j).gt.ZERO)
               savg  = HALF * (splft + sprgt)
               splus = merge(splus, savg, abs(umac(i+1,j)) .gt. eps)
               
               ! extrapolate to transverse faces on lo side of cell
               smrgt = s(i  ,j,n) - (HALF + dt2*u(i  ,j,1)/hx) * slopex(i  ,j,n)
               smlft = s(i-1,j,n) + (HALF - dt2*u(i-1,j,1)/hx) * slopex(i-1,j,n)
               
               ! apply boundary conditions on lo-x side
               smrgt = merge(s(is-1,j,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
               smlft = merge(s(is-1,j,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
               if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                  if (is_vel .and. n .eq. 1) then
                     smlft = ZERO
                     smrgt = ZERO
                  elseif (is_vel .and. n .ne. 1) then
                     smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                     smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                  else
                     smlft = smrgt
                  endif
               endif
               
               ! upwind based on umac
               sminus = merge(smlft,smrgt,umac(i,j).gt.ZERO)
               savg   = HALF * (smlft + smrgt)
               sminus = merge(sminus, savg, abs(umac(i,j)) .gt. eps)
               
               ! create transverse plus force term
               st = force(i,j,n) - &
                    HALF * (umac(i,j)+umac(i+1,j))*(splus - sminus) / hx
               
               ! extrapolate to normal faces with transverse and force contribution
               s_b(j+1)= s(i,j,n) + (HALF-dt2*vmac(i,j+1)/hy)*slopey(i,j,n) + dt2*st
               s_t(j  )= s(i,j,n) - (HALF+dt2*vmac(i,j)/hy)*slopey(i,j,n) + dt2*st
            enddo
            
            ! upwind based on vmac
            do j = js, je+1 
               sedgey(i,j,n)=merge(s_b(j),s_t(j),vmac(i,j).gt.ZERO)
               savg = HALF*(s_b(j) + s_t(j))
               sedgey(i,j,n)=merge(savg,sedgey(i,j,n),abs(vmac(i,j)) .lt. eps)
            enddo
            
            ! apply boundary conditions on lo-y side of domain
            if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
               if (is_vel .and. n .eq. 2) then
                  sedgey(i,js,n) = ZERO
               elseif (is_vel .and. n .ne. 2) then
                  sedgey(i,js,n) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
               else 
                  sedgey(i,js,n) = s_t(js)
               endif
            elseif (phys_bc(2,1) .eq. INLET) then
               sedgey(i,js,n) = s(i,js-1,n)
            elseif (phys_bc(2,1) .eq. OUTLET) then
               sedgey(i,js,n) = s_t(js)
            endif
            
            ! apply boundary conditions on hi-y side of domain
            if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
               if (is_vel .and. n .eq. 2) then
                  sedgey(i,je+1,n) = ZERO
               elseif (is_vel .and. n .ne. 2) then
                  sedgey(i,je+1,n) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
               else 
                  sedgey(i,je+1,n) = s_b(je+1)
               endif
            elseif (phys_bc(2,2) .eq. INLET) then
               sedgey(i,je+1,n) = s(i,je+1,n)
            elseif (phys_bc(2,2) .eq. OUTLET) then
               sedgey(i,je+1,n) = s_b(je+1)
            endif
         enddo
      enddo
      
      deallocate(s_l)
      deallocate(s_r)
      deallocate(s_b)
      deallocate(s_t)

      deallocate(slopex)
      deallocate(slopey)

      end subroutine mkflux_2d

      subroutine mkflux_3d(s,u,sedgex,sedgey,sedgez,&
                           umac,vmac,wmac, &
                           force,lo,dx,dt,is_vel, &
                           phys_bc,adv_bc,ng)

      integer, intent(in) :: lo(:),ng

      real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) :: adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:,:,:,:)
      real(kind=dp_t), allocatable::  slopey(:,:,:,:)
      real(kind=dp_t), allocatable::  slopez(:,:,:,:)

      real(kind=dp_t) hx, hy, hz, dt2, dt4, dt6
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

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

      ! Normal predictor states.
      ! Allocated from lo:hi+1 in the normal direction
      ! lo-1:hi+1 in the transverse directions
      allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

      allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

      allocate(slz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(srz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(simhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

      ! These are transverse terms.  The size allocation is tricky.
      ! lo:hi+1 in normal direction
      ! lo:hi in transverse direction
      ! lo-1:hi+1 in unused direction
      allocate(slxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
      allocate(srxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
      allocate(simhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

      allocate(slxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
      allocate(srxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
      allocate(simhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

      allocate(slyx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(sryx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(simhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

      allocate(slyz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(sryz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
      allocate(simhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

      allocate(slzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(srzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
      allocate(simhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

      allocate(slzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
      allocate(srzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
      allocate(simhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

      ! Final edge states.
      ! lo:hi+1 in the normal direction
      ! lo:hi in the transverse directions
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

      dt2 = HALF*dt
      dt4 = dt/4.0d0
      dt6 = dt/6.0d0

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      ! Compute eps, which is relative to the max mac velocity
      umax = abs(umac(is,js,ks))
      do k = ks,ke
         do j = js,je
            do i = is,ie+1
               umax = max(umax,abs(umac(i,j,k)))
            end do
         end do
      end do
      do k = ks,ke
         do j = js,je+1
            do i = is,ie
               umax = max(umax,abs(vmac(i,j,k)))
            end do
         end do
      end do
      do k = ks,ke+1
         do j = js,je
            do i = is,ie
               umax = max(umax,abs(wmac(i,j,k)))
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
                  slx(i,j,k) = s(i-1,j,k,n) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,n)
                  srx(i,j,k) = s(i  ,j,k,n) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i,  j,k,n)
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slx(i,j,k) = merge(s(is-1,j,k,n),slx(i,j,k),phys_bc(1,1) .eq. INLET)
                     srx(i,j,k) = merge(s(is-1,j,k,n),srx(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 1) then
                           slx(i,j,k) = ZERO
                           srx(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 1) then
                           slx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        else
                           slx(i,j,k) = srx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slx(i,j,k) = merge(s(ie+1,j,k,n),slx(i,j,k),phys_bc(1,2) .eq. INLET)
                     srx(i,j,k) = merge(s(ie+1,j,k,n),srx(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 1) then
                           slx(i,j,k) = ZERO
                           srx(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 1) then
                           slx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        else
                           srx(i,j,k) = slx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhx by solving Riemann problem
                  simhx(i,j,k) = merge(slx(i,j,k),srx(i,j,k),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(slx(i,j,k)+srx(i,j,k))
                  simhx(i,j,k) = merge(simhx(i,j,k),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate y-faces
         do k=ks-1,ke+1
            do j=js,je+1
               do i=is-1,ie+1
                  ! make sly, sry with 1D extrapolation
                  sly(i,j,k) = s(i,j-1,k,n) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,n)
                  sry(i,j,k) = s(i,j,  k,n) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j,  k,n)
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     sly(i,j,k) = merge(s(is,j-1,k,n),sly(i,j,k),phys_bc(2,1) .eq. INLET)
                     sry(i,j,k) = merge(s(is,j-1,k,n),sry(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           sly(i,j,k) = ZERO
                           sry(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           sly(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sry(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           sly(i,j,k) = sry(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     sly(i,j,k) = merge(s(i,je+1,k,n),sly(i,j,k),phys_bc(2,2) .eq. INLET)
                     sry(i,j,k) = merge(s(i,je+1,k,n),sry(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           sly(i,j,k) = ZERO
                           sry(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           sly(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sry(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sry(i,j,k) = sly(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhy by solving Riemann problem
                  simhy(i,j,k) = merge(sly(i,j,k),sry(i,j,k),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(sly(i,j,k)+sry(i,j,k))
                  simhy(i,j,k) = merge(simhy(i,j,k),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate z-faces
         do k=ks,ke+1
            do j=js-1,je+1
               do i=is-1,ie+1
                  ! make slz, srz with 1D extrapolation
                  slz(i,j,k) = s(i,j,k-1,n) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,n)
                  srz(i,j,k) = s(i,j,k,  n) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k,  n)
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slz(i,j,k) = merge(s(is,j,k-1,n),slz(i,j,k),phys_bc(3,1) .eq. INLET)
                     srz(i,j,k) = merge(s(is,j,k-1,n),srz(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 3) then
                           slz(i,j,k) = ZERO
                           srz(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 3) then
                           slz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        else
                           slz(i,j,k) = srz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slz(i,j,k) = merge(s(i,j,ke+1,n),slz(i,j,k),phys_bc(3,2) .eq. INLET)
                     srz(i,j,k) = merge(s(i,j,ke+1,n),srz(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 3) then
                           slz(i,j,k) = ZERO
                           srz(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 3) then
                           slz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        else
                           srz(i,j,k) = slz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhz by solving Riemann problem
                  simhz(i,j,k) = merge(slz(i,j,k),srz(i,j,k),wmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slz(i,j,k)+srz(i,j,k))
                  simhz(i,j,k) = merge(simhz(i,j,k),savg,abs(wmac(i,j,k)) .gt. eps)
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
                  slxy(i,j,k) = slx(i,j,k) - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                  srxy(i,j,k) = srx(i,j,k) - (dt6/hy)*(vmac(i,  j+1,k)+vmac(i,  j,k))*(simhy(i,  j+1,k)-simhy(i,  j,k))
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slxy(i,j,k) = merge(s(is-1,j,k,n),slxy(i,j,k),phys_bc(1,1) .eq. INLET)
                     srxy(i,j,k) = merge(s(is-1,j,k,n),srxy(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 1) then
                           slxy(i,j,k) = ZERO
                           srxy(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 1) then
                           slxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        else
                           slxy(i,j,k) = srxy(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slxy(i,j,k) = merge(s(ie+1,j,k,n),slxy(i,j,k),phys_bc(1,2) .eq. INLET)
                     srxy(i,j,k) = merge(s(ie+1,j,k,n),srxy(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 1) then
                           slxy(i,j,k) = ZERO
                           srxy(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 1) then
                           slxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        else
                           srxy(i,j,k) = slxy(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhxy by solving Riemann problem
                  simhxy(i,j,k) = merge(slxy(i,j,k),srxy(i,j,k),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(slxy(i,j,k)+srxy(i,j,k))
                  simhxy(i,j,k) = merge(simhxy(i,j,k),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate xz faces
         do k=ks,ke
            do j=js-1,je+1
               do i=is,ie+1
                  ! make slxz, srxz by updating 1D extrapolation
                  slxz(i,j,k) = slx(i,j,k) - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k))*(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                  srxz(i,j,k) = srx(i,j,k) - (dt6/hz)*(wmac(i,  j,k+1)+wmac(i,  j,k))*(simhz(i,  j,k+1)-simhz(i,  j,k))
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slxz(i,j,k) = merge(s(is-1,j,k,n),slxz(i,j,k),phys_bc(1,1) .eq. INLET)
                     srxz(i,j,k) = merge(s(is-1,j,k,n),srxz(i,j,k),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 1) then
                           slxz(i,j,k) = ZERO
                           srxz(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 1) then
                           slxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        else
                           slxz(i,j,k) = srxz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slxz(i,j,k) = merge(s(ie+1,j,k,n),slxz(i,j,k),phys_bc(1,2) .eq. INLET)
                     srxz(i,j,k) = merge(s(ie+1,j,k,n),srxz(i,j,k),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 1) then
                           slxz(i,j,k) = ZERO
                           srxz(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 1) then
                           slxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                        else
                           srxz(i,j,k) = slxz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhxz by solving Riemann problem
                  simhxz(i,j,k) = merge(slxz(i,j,k),srxz(i,j,k),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(slxz(i,j,k)+srxz(i,j,k))
                  simhxz(i,j,k) = merge(simhxz(i,j,k),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo

         ! loop over appropriate yx faces
         do k=ks-1,ke+1
            do j=js,je+1
               do i=is,ie
                  ! make slyx, sryx by updating 1D extrapolation
                  slyx(i,j,k) = sly(i,j,k) - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                  sryx(i,j,k) = sry(i,j,k) - (dt6/hx)*(umac(i+1,j,  k)+umac(i,j,  k))*(simhx(i+1,j,  k)-simhx(i,j,  k))
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     slyx(i,j,k) = merge(s(is,j-1,k,n),slyx(i,j,k),phys_bc(2,1) .eq. INLET)
                     sryx(i,j,k) = merge(s(is,j-1,k,n),sryx(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           slyx(i,j,k) = ZERO
                           sryx(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           slyx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sryx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           slyx(i,j,k) = sryx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     slyx(i,j,k) = merge(s(i,je+1,k,n),slyx(i,j,k),phys_bc(2,2) .eq. INLET)
                     sryx(i,j,k) = merge(s(i,je+1,k,n),sryx(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           slyx(i,j,k) = ZERO
                           sryx(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           slyx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sryx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sryx(i,j,k) = slyx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhyx by solving Riemann problem
                  simhyx(i,j,k) = merge(slyx(i,j,k),sryx(i,j,k),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slyx(i,j,k)+sryx(i,j,k))
                  simhyx(i,j,k) = merge(simhyx(i,j,k),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate yz faces
         do k=ks,ke
            do j=js,je+1
               do i=is-1,ie+1
                  ! make slyz, sryz by updating 1D extrapolation
                  slyz(i,j,k) = sly(i,j,k) - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k))*(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                  sryz(i,j,k) = sry(i,j,k) - (dt6/hz)*(wmac(i,j,  k+1)+wmac(i,j,  k))*(simhz(i,j,  k+1)-simhz(i,j,  k))
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     slyz(i,j,k) = merge(s(is,j-1,k,n),slyz(i,j,k),phys_bc(2,1) .eq. INLET)
                     sryz(i,j,k) = merge(s(is,j-1,k,n),sryz(i,j,k),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           slyz(i,j,k) = ZERO
                           sryz(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           slyz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sryz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           slyz(i,j,k) = sryz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     slyz(i,j,k) = merge(s(i,je+1,k,n),slyz(i,j,k),phys_bc(2,2) .eq. INLET)
                     sryz(i,j,k) = merge(s(i,je+1,k,n),sryz(i,j,k),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           slyz(i,j,k) = ZERO
                           sryz(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           slyz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sryz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sryz(i,j,k) = slyz(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhyz by solving Riemann problem
                  simhyz(i,j,k) = merge(slyz(i,j,k),sryz(i,j,k),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slyz(i,j,k)+sryz(i,j,k))
                  simhyz(i,j,k) = merge(simhyz(i,j,k),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate zx faces
         do k=ks,ke+1
            do j=js-1,je+1
               do i=is,ie
                  ! make slzx, srzx by updating 1D extrapolation
                  slzx(i,j,k) = slz(i,j,k) - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                  srzx(i,j,k) = srz(i,j,k) - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhx(i+1,j,k  )-simhx(i,j,k  ))
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slzx(i,j,k) = merge(s(is,j,k-1,n),slzx(i,j,k),phys_bc(3,1) .eq. INLET)
                     srzx(i,j,k) = merge(s(is,j,k-1,n),srzx(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 3) then
                           slzx(i,j,k) = ZERO
                           srzx(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 3) then
                           slzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        else
                           slzx(i,j,k) = srzx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slzx(i,j,k) = merge(s(i,j,ke+1,n),slzx(i,j,k),phys_bc(3,2) .eq. INLET)
                     srzx(i,j,k) = merge(s(i,j,ke+1,n),srzx(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 3) then
                           slzx(i,j,k) = ZERO
                           srzx(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 3) then
                           slzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        else
                           srzx(i,j,k) = slzx(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhzx by solving Riemann problem
                  simhzx(i,j,k) = merge(slzx(i,j,k),srzx(i,j,k),wmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slzx(i,j,k)+srzx(i,j,k))
                  simhzx(i,j,k) = merge(simhzx(i,j,k),savg,abs(wmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! loop over appropriate zy faces
         do k=ks,ke+1
            do j=js,je
               do i=is-1,ie+1
                  ! make slzy, srzy by updating 1D extrapolation
                  slzy(i,j,k) = slz(i,j,k) - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                  srzy(i,j,k) = srz(i,j,k) - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhy(i,j+1,k  )-simhy(i,j,k  ))
                  
                  ! impose lo side bc's
                  if(k .eq. ks) then
                     slzy(i,j,k) = merge(s(is,j,k-1,n),slzy(i,j,k),phys_bc(3,1) .eq. INLET)
                     srzy(i,j,k) = merge(s(is,j,k-1,n),srzy(i,j,k),phys_bc(3,1) .eq. INLET)
                     if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 3) then
                           slzy(i,j,k) = ZERO
                           srzy(i,j,k) = ZERO
                        else if(is_vel .and. n .ne. 3) then
                           slzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           srzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                        else
                           slzy(i,j,k) = srzy(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(k .eq. ke+1) then
                     slzy(i,j,k) = merge(s(i,j,ke+1,n),slzy(i,j,k),phys_bc(3,2) .eq. INLET)
                     srzy(i,j,k) = merge(s(i,j,ke+1,n),srzy(i,j,k),phys_bc(3,2) .eq. INLET)
                     if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 3) then
                           slzy(i,j,k) = ZERO
                           srzy(i,j,k) = ZERO
                        else if (is_vel .and. n .ne. 3) then
                           slzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                           srzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                        else
                           srzy(i,j,k) = slzy(i,j,k)
                        endif
                     endif
                  endif
                  
                  ! make simhzy by solving Riemann problem
                  simhzy(i,j,k) = merge(slzy(i,j,k),srzy(i,j,k),wmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slzy(i,j,k)+srzy(i,j,k))
                  simhzy(i,j,k) = merge(simhzy(i,j,k),savg,abs(wmac(i,j,k)) .gt. eps)
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
                  sedgelx(i,j,k) = s(i-1,j,k,n) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,n) &
                       - (dt4/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhyz(i-1,j+1,k)-simhyz(i-1,j,k)) &
                       - (dt4/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k))*(simhzy(i-1,j,k+1)-simhzy(i-1,j,k)) &
                       + dt2*force(i-1,j,k,n)
                  sedgerx(i,j,k) = s(i,  j,k,n) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i,  j,k,n) &
                       - (dt4/hy)*(vmac(i,  j+1,k)+vmac(i,  j,k))*(simhyz(i,  j+1,k)-simhyz(i,  j,k)) &
                       - (dt4/hz)*(wmac(i,  j,k+1)+wmac(i,  j,k))*(simhzy(i,  j,k+1)-simhzy(i,  j,k)) &
                       + dt2*force(i,j,k,n)
                  
                  ! make sedgex by solving Riemann problem
                  ! boundary conditions enforced outside of i,j,k loop
                  sedgex(i,j,k,n) = merge(sedgelx(i,j,k),sedgerx(i,j,k),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
                  sedgex(i,j,k,n) = merge(sedgex(i,j,k,n),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo

         ! sedgex boundary conditions
         do k=ks,ke
            do j=js,je
               ! lo side
               if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(is,j,k,n) = ZERO
                  elseif (is_vel .and. n .ne. 1) then
                     sedgex(is,j,k,n) = merge(ZERO,sedgerx(is,j,k),phys_bc(1,1).eq.NO_SLIP_WALL)
                  else
                     sedgex(is,j,k,n) = sedgerx(is,j,k)
                  endif
               elseif (phys_bc(1,1) .eq. INLET) then
                  sedgex(is,j,k,n) = s(is-1,j,k,n)
               elseif (phys_bc(1,1) .eq. OUTLET) then
                  if (is_vel .and. n.eq.1) then
                     sedgex(is,j,k,n) = MIN(sedgerx(is,j,k),ZERO)
                  else
                     sedgex(is,j,k,n) = sedgerx(is,j,k)
                  end if
               endif

               ! hi side
               if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(ie+1,j,k,n) = ZERO
                  else if (is_vel .and. n .ne. 1) then
                     sedgex(ie+1,j,k,n) = merge(ZERO,sedgelx(ie+1,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                  else 
                     sedgex(ie+1,j,k,n) = sedgelx(ie+1,j,k)
                  endif
               elseif (phys_bc(1,2) .eq. INLET) then
                  sedgex(ie+1,j,k,n) = s(ie+1,j,k,n)
               elseif (phys_bc(1,2) .eq. OUTLET) then
                  if (is_vel .and. n.eq.1) then
                     sedgex(ie+1,j,k,n) = MAX(sedgelx(ie+1,j,k),ZERO)
                  else
                     sedgex(ie+1,j,k,n) = sedgelx(ie+1,j,k)
                  end if
               endif
            enddo
         enddo
         
         ! loop over appropriate y-faces
         do k=ks,ke
            do j=js,je+1
               do i=is,ie
                  ! make sedgely, sedgery
                  sedgely(i,j,k) = s(i,j-1,k,n) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,n) &
                       - (dt4/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhxz(i+1,j-1,k)-simhxz(i,j-1,k)) &
                       - (dt4/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k))*(simhzx(i,j-1,k+1)-simhzx(i,j-1,k)) &
                       + dt2*force(i,j-1,k,n)
                  sedgery(i,j,k) = s(i,j,  k,n) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j,  k,n) &
                       - (dt4/hx)*(umac(i+1,j,  k)+umac(i,j,  k))*(simhxz(i+1,j,  k)-simhxz(i,j,  k)) &
                       - (dt4/hz)*(wmac(i,j,  k+1)+wmac(i,j,  k))*(simhzx(i,j,  k+1)-simhzx(i,j,  k)) &
                       + dt2*force(i,j,k,n)
                  
                  ! make sedgey by solving Riemann problem
                  ! boundary conditions enforced outside of i,j,k loop
                  sedgey(i,j,k,n) = merge(sedgely(i,j,k),sedgery(i,j,k),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
                  sedgey(i,j,k,n) = merge(sedgey(i,j,k,n),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo
         
         ! sedgey boundary conditions
         do k=ks,ke
            do i=is,ie
               ! lo side
               if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgey(i,js,k,n) = ZERO
                  elseif (is_vel .and. n .ne. 2) then
                     sedgey(i,js,k,n) = merge(ZERO,sedgery(i,js,k),phys_bc(2,1).eq.NO_SLIP_WALL)
                  else 
                     sedgey(i,js,k,n) = sedgery(i,js,k)
                  endif
               elseif (phys_bc(2,1) .eq. INLET) then
                  sedgey(i,js,k,n) = s(i,js-1,k,n)
               elseif (phys_bc(2,1) .eq. OUTLET) then
                  if (is_vel .and. n.eq.2) then
                     sedgey(i,js,k,n) = MIN(sedgery(i,js,k),ZERO)
                  else
                     sedgey(i,js,k,n) = sedgery(i,js,k)
                  end if
               endif
               
               ! hi side
               if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgey(i,je+1,k,n) = ZERO
                  elseif (is_vel .and. n .ne. 2) then
                     sedgey(i,je+1,k,n) = merge(ZERO,sedgely(i,je+1,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                  else 
                     sedgey(i,je+1,k,n) = sedgely(i,je+1,k)
                  endif
               elseif (phys_bc(2,2) .eq. INLET) then
                  sedgey(i,je+1,k,n) = s(i,je+1,k,n)
               elseif (phys_bc(2,2) .eq. OUTLET) then
                  if (is_vel .and. n.eq.2) then
                     sedgey(i,je+1,k,n) = MAX(sedgely(i,je+1,k),ZERO)
                  else
                     sedgey(i,je+1,k,n) = sedgely(i,je+1,k)
                  end if
               endif
            enddo
         enddo

         ! loop over appropriate z-faces
         do k=ks,ke+1
            do j=js,je
               do i=is,ie
                  ! make sedgelz, sedgerz
                  sedgelz(i,j,k) = s(i,j,k-1,n) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,n) &
                       - (dt4/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhxy(i+1,j,k-1)-simhxy(i,j,k-1)) &
                       - (dt4/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhyx(i,j+1,k-1)-simhyx(i,j,k-1)) &
                       + dt2*force(i,j,k-1,n)
                  sedgerz(i,j,k) = s(i,j,k,  n) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k,  n) &
                       - (dt4/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhxy(i+1,j,k  )-simhxy(i,j,k  )) &
                       - (dt4/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhyx(i,j+1,k  )-simhyx(i,j,k  )) &
                       + dt2*force(i,j,k,n)
                  
                  ! make sedgez by solving Riemann problem
                  ! boundary conditions enforced outside of i,j,k loop
                  sedgez(i,j,k,n) = merge(sedgelz(i,j,k),sedgerz(i,j,k),wmac(i,j,k) .gt. ZERO)
                  savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
                  sedgez(i,j,k,n) = merge(sedgez(i,j,k,n),savg,abs(wmac(i,j,k)) .gt. eps)
               enddo
            enddo
         enddo

         ! sedgez boundary conditions
         do j=js,je
            do i=is,ie
               ! lo side
               if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgez(i,j,ks,n) = ZERO
                  elseif (is_vel .and. n .ne. 2) then
                     sedgez(i,j,ks,n) = merge(ZERO,sedgerz(i,j,ks),phys_bc(3,1).eq.NO_SLIP_WALL)
                  else 
                     sedgez(i,j,ks,n) = sedgerz(i,j,ks)
                  endif
               elseif (phys_bc(3,1) .eq. INLET) then
                  sedgez(i,j,ks,n) = s(i,j,ks-1,n)
               elseif (phys_bc(3,1) .eq. OUTLET) then
                  if (is_vel .and. n.eq.3) then
                     sedgez(i,j,ks,n) = MIN(sedgerz(i,j,ks),ZERO)
                  else
                     sedgez(i,j,ks,n) = sedgerz(i,j,ks)
                  end if
               endif

               ! hi side
               if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgez(i,j,ke+1,n) = ZERO
                  elseif (is_vel .and. n .ne. 2) then
                     sedgez(i,j,ke+1,n) = merge(ZERO,sedgelz(i,j,ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
                  else 
                     sedgez(i,j,ke+1,n) = sedgelz(i,j,ke+1)
                  endif
               elseif (phys_bc(3,2) .eq. INLET) then
                  sedgez(i,j,ke+1,n) = s(i,j,ke+1,n)
               elseif (phys_bc(3,2) .eq. OUTLET) then
                  if (is_vel .and. n.eq.3) then
                     sedgez(i,j,ke+1,n) = MAX(sedgelz(i,j,ke+1),ZERO)
                  else
                     sedgez(i,j,ke+1,n) = sedgelz(i,j,ke+1)
                  end if
               endif
            enddo
         enddo

      enddo ! end loop over components
      
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

      end subroutine mkflux_3d

end module mkflux_module
