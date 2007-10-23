module mkflux_macvel_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module

  implicit none

contains

      subroutine mkflux_macvel_2d(s,u,sedgex,sedgey,uadv,vadv,utrans,vtrans,&
                                  force,lo,dx,dt,is_vel,is_cons,&
                                  phys_bc,adv_bc,ng)

      integer, intent(in) :: lo(2),ng

      real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) ::   uadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::   vadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel
      logical        ,intent(in) :: is_cons(:)

      real(kind=dp_t), allocatable::  slopex(:,:,:),slopey(:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:)

!     Local variables
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
         
      umax = abs(uadv(is,js))
      do j = js,je
         do i = is,ie+1
            umax = max(umax,abs(uadv(i,j)))
         end do
      end do
      do j = js,je+1
         do i = is,ie
            umax = max(umax,abs(vadv(i,j)))
         end do
      end do
      
      eps = abs_eps * umax
      
      !
      !     Loop for fluxes on x-edges.
      !
      do n = 1,ncomp
         do j = js,je 
            do i = is-1,ie+1 
               spbot = s(i,j  ,n) + (HALF - dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
               !    $            + dth * force(i,  j,n)
               sptop = s(i,j+1,n) - (HALF + dth*u(i,j+1,2)/hy) * slopey(i,j+1,n)
               !    $            + dth * force(i,j+1,n)
               
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
               
               splus = merge(spbot,sptop,vtrans(i,j+1).gt.ZERO)
               savg  = HALF * (spbot + sptop)
               splus = merge(splus, savg, abs(vtrans(i,j+1)) .gt. eps)
               
               smtop = s(i,j  ,n) - (HALF + dth*u(i,j  ,2)/hy) * slopey(i,j  ,n)
               !    $            + dth * force(i,j  ,n)
               smbot = s(i,j-1,n) + (HALF - dth*u(i,j-1,2)/hy) * slopey(i,j-1,n)
               !    $            + dth * force(i,j-1,n)
               
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
               
               sminus = merge(smbot,smtop,vtrans(i,j).gt.ZERO)
               savg   = HALF * (smbot + smtop)
               sminus = merge(sminus, savg, abs(vtrans(i,j)) .gt. eps)
               
               st = force(i,j,n) - &
                    HALF * (vtrans(i,j)+vtrans(i,j+1))*(splus - sminus) / hy
               
               s_l(i+1)= s(i,j,n) + (HALF-dth*uadv(i+1,j)/hx)*slopex(i,j,1) + dth*st
               s_r(i  )= s(i,j,n) - (HALF+dth*uadv(i,j)/hx)*slopex(i,j,1) + dth*st

            enddo
            
            do i = is, ie+1 
               sedgex(i,j,n)=merge(s_l(i),s_r(i),uadv(i,j).gt.ZERO)
               savg = HALF*(s_r(i) + s_l(i))
               sedgex(i,j,n)=merge(savg,sedgex(i,j,n),abs(uadv(i,j)) .lt. eps)
            enddo
            
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
               splft = s(i,j  ,n) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
               !    $            + dth * force(i  ,j,n)
               sprgt = s(i+1,j,n) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,n)
               !    $            + dth * force(i+1,j,n)
               
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
               
               splus = merge(splft,sprgt,utrans(i+1,j).gt.ZERO)
               savg  = HALF * (splft + sprgt)
               splus = merge(splus, savg, abs(utrans(i+1,j)) .gt. eps)
               
               smrgt = s(i  ,j,n) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,n)
               !    $            + dth * force(i  ,j,n)
               smlft = s(i-1,j,n) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,n)
               !    $            + dth * force(i-1,j,n)
               
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
               
               sminus = merge(smlft,smrgt,utrans(i,j).gt.ZERO)
               savg   = HALF * (smlft + smrgt)
               sminus = merge(sminus, savg, abs(utrans(i,j)) .gt. eps)
               
               st = force(i,j,n) - &
                    HALF * (utrans(i,j)+utrans(i+1,j))*(splus - sminus) / hx
               
               s_b(j+1)= s(i,j,n) + (HALF-dth*vadv(i,j+1)/hy)*slopey(i,j,1) + dth*st
               s_t(j  )= s(i,j,n) - (HALF+dth*vadv(i,j)/hy)*slopey(i,j,1) + dth*st
               
            enddo
            
            do j = js, je+1 
               sedgey(i,j,n)=merge(s_b(j),s_t(j),vadv(i,j).gt.ZERO)
               savg = HALF*(s_b(j) + s_t(j))
               sedgey(i,j,n)=merge(savg,sedgey(i,j,n),abs(vadv(i,j)) .lt. eps)
            enddo
            
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

      end subroutine mkflux_macvel_2d

      subroutine mkflux_macvel_3d(s,u,sedgex,sedgey,sedgez,&
                           uadv,vadv,wadv,utrans,vtrans,wtrans, &
                           force,lo,dx,dt,is_vel,is_cons, &
                           phys_bc,adv_bc,ng)


      implicit none

      integer, intent(in) :: lo(:),ng

      real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) ::   uadv(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   vadv(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::   wadv(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) :: wtrans(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel
      logical        ,intent(in) :: is_cons(:)

      real(kind=dp_t), allocatable::  slopex(:,:,:,:)
      real(kind=dp_t), allocatable::  slopey(:,:,:,:)
      real(kind=dp_t), allocatable::  slopez(:,:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)

!     Local variables
      real(kind=dp_t) hx, hy, hz, dth
      real(kind=dp_t) splus,sminus,st,str,savg
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
      logical test

      real(kind=dp_t) :: abs_eps,eps,umax

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke,n
      integer :: slope_order = 4
      integer :: ncomp

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

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      umax = abs(uadv(is,js,ks))
      do k = ks,ke
         do j = js,je
            do i = is,ie+1
               umax = max(umax,abs(uadv(i,j,k)))
            end do
         end do
      end do
      do k = ks,ke
         do j = js,je+1
            do i = is,ie
               umax = max(umax,abs(vadv(i,j,k)))
            end do
         end do
      end do
      do k = ks,ke+1
         do j = js,je
            do i = is,ie
               umax = max(umax,abs(wadv(i,j,k)))
            end do
         end do
      end do
      
      eps = abs_eps * umax
      
      !     loop for x fluxes
      do n = 1,ncomp
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
                  
                  sptop = merge(s(i,je+1,k,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  spbot = merge(s(i,je+1,k,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  
                  if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 2) then
                        sptop = ZERO
                        spbot = ZERO
                     else if (is_vel .and. n .ne. 2) then
                        sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                        spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                     else
                        sptop = spbot
                     endif
                  endif
                  
                  splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                  savg  = HALF * (spbot + sptop)
                  splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                  
                  smtop = s(i,j  ,k,n) - (HALF + dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
                  !    $            + dth * force(i,j  ,k,n)
                  smbot = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n)
                  !    $            + dth * force(i,j-1,k,n)
                  
                  smtop = merge(s(i,js-1,k,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  smbot = merge(s(i,js-1,k,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  
                  if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 2) then
                        smtop = ZERO
                        smbot = ZERO
                     else if (is_vel .and. n .ne. 2) then
                        smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                        smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                     else
                        smbot = smtop
                     endif
                  endif
                  
                  sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smbot + smtop)
                  sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                  
                  str =  HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
                  
!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Z-DIRECTION
!        ******************************************************************

                  spbot = s(i,j,k  ,n) + (HALF - dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
                  !    $            + dth * force(i,j,k  ,n)
                  sptop = s(i,j,k+1,n) - (HALF + dth*u(i,j,k+1,3)/hz)*slopez(i,j,k+1,n)
                  !    $            + dth * force(i,j,k+1,n)
                  
                  sptop = merge(s(i,j,ke+1,n),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                  spbot = merge(s(i,j,ke+1,n),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                  
                  if (k .eq. ke .and. (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 3) then
                        sptop = ZERO
                        spbot = ZERO
                     else if (is_vel .and. n .ne. 3) then
                        sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                        spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                     else
                        sptop = spbot
                     endif
                  endif
                  
                  splus = merge(spbot,sptop,wtrans(i,j,k+1).gt.ZERO)
                  savg  = HALF * (spbot + sptop)
                  splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                  
                  smtop = s(i,j,k  ,n) - (HALF + dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
                  !    $            + dth * force(i,j,k  ,n)
                  smbot = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n)
                  !    $            + dth * force(i,j,k-1,n)
                  
                  smtop = merge(s(i,j,ks-1,n),smtop,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                  smbot = merge(s(i,j,ks-1,n),smbot,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                  
                  if (k .eq. ks .and. (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 3) then
                        smtop = ZERO
                        smbot = ZERO
                     else if (is_vel .and. n .ne. 3) then
                        smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                        smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                     else 
                        smbot = smtop
                     endif
                  endif
                  
                  sminus = merge(smbot,smtop,wtrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smbot + smtop)
                  sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                  
                  str = str + HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))* &
                       (splus - sminus) / hz
                  
!        ******************************************************************
!        MAKE LEFT AND RIGHT STATES
!        ******************************************************************
                  
                  st = force(i,j,k,n) - str
                  
                  s_l(i+1)= s(i,j,k,n) + (HALF-dth*uadv(i+1,j,k))*slopex(i,j,k,1) + dth*st
                  s_r(i  )= s(i,j,k,n) - (HALF+dth*uadv(i  ,j,k))*slopex(i,j,k,1) + dth*st

               enddo
               
               do i = is, ie+1 
                  sedgex(i,j,k,n)=merge(s_l(i),s_r(i),uadv(i,j,k).gt.ZERO)
                  savg = HALF*(s_r(i) + s_l(i))
                  sedgex(i,j,k,n)=merge(savg,sedgex(i,j,k,n),abs(uadv(i,j,k)) .lt. eps)
               enddo
               
               if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(is,j,k,n) = ZERO
                  else if (is_vel .and. n .ne. 1) then
                     sedgex(is,j,k,n) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
                  else 
                     sedgex(is,j,k,n) = s_r(is)
                  endif
               elseif (phys_bc(1,1) .eq. INLET) then
                  sedgex(is,j,k,n) = s(is-1,j,k,n)
               elseif (phys_bc(1,1) .eq. OUTLET) then
                  sedgex(is,j,k,n) = s_r(is)
               endif
               
               if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(ie+1,j,k,n) = ZERO
                  else if (is_vel .and. n .ne. 1) then
                     sedgex(ie+1,j,k,n) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
                  else 
                     sedgex(ie+1,j,k,n) = s_l(ie+1)
                  endif
               elseif (phys_bc(1,2) .eq. INLET) then
                  sedgex(ie+1,j,k,n) = s(ie+1,j,k,n)
               elseif (phys_bc(1,2) .eq. OUTLET) then
                  sedgex(ie+1,j,k,n) = s_l(ie+1)
               endif
               
            enddo
         enddo
      enddo

!        ******************************************************************
!        ******************************************************************
!        ******************************************************************

!     loop for y fluxes
      do n = 1,ncomp
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
                  
                  sprgt = merge(s(ie+1,j,k,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  splft = merge(s(ie+1,j,k,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  
                  if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 1) then
                        sprgt = ZERO
                        splft = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                        splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                     else
                        sprgt = splft
                     endif
                  endif
                  
                  splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                  savg  = HALF * (splft + sprgt)
                  splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                  
                  smrgt = s(i  ,j,k,n) - (HALF + dth*u(i  ,j,k,1)/hx)*slopex(i  ,j,k,n)
                  !    $            + dth * force(i  ,j,k,n)
                  smlft = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n)
                  !    $            + dth * force(i-1,j,k,n)
                  
                  smrgt = merge(s(is-1,j,k,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  smlft = merge(s(is-1,j,k,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  
                  if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 1) then
                        smrgt = ZERO
                        smlft = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                        smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                     else
                        smlft = smrgt
                     endif
                  endif
                  
                  sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smlft + smrgt)
                  sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                  
                  str    = HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                  
!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Z-DIRECTION
!        ******************************************************************

                  spbot = s(i,j,k  ,n) + (HALF - dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
                  !    $            + dth * force(i,j,k  ,n)
                  sptop = s(i,j,k+1,n) - (HALF + dth*u(i,j,k+1,3)/hz)*slopez(i,j,k+1,n)
                  !    $            + dth * force(i,j,k+1,n)
                  
                  sptop = merge(s(i,j,ke+1,n),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                  spbot = merge(s(i,j,ke+1,n),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                  
                  if (k .eq. ke .and. (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 3) then
                        sptop = ZERO
                        spbot = ZERO
                     else if (is_vel .and. n .ne. 3) then
                        sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                        spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                     else
                        sptop = spbot
                     endif
                  endif
                  
                  splus = merge(spbot,sptop,wtrans(i,j,k+1).gt.ZERO)
                  savg  = HALF * (spbot + sptop)
                  splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                  
                  smtop = s(i,j,k  ,n) - (HALF + dth*u(i,j,k  ,3)/hz)*slopez(i,j,k  ,n)
                  !    $            + dth * force(i,j,k  ,n)
                  smbot = s(i,j,k-1,n) + (HALF - dth*u(i,j,k-1,3)/hz)*slopez(i,j,k-1,n)
                  !    $            + dth * force(i,j,k-1,n)
                  
                  smtop = merge(s(i,j,ks-1,n),smtop,k.eq.ke .and. phys_bc(3,1) .eq. INLET)
                  smbot = merge(s(i,j,ks-1,n),smbot,k.eq.ke .and. phys_bc(3,1) .eq. INLET)
                  
                  if (k .eq. ks .and. (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 3) then
                        smtop = ZERO
                        smbot = ZERO
                     else if (is_vel .and. n .ne. 3) then
                        smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                        smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                     else
                        smbot = smtop
                     endif
                  endif
                  
                  sminus = merge(smbot,smtop,wtrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smbot + smtop)
                  sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                  
                  str = str + HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz

!        ******************************************************************
!        MAKE TOP AND BOTTOM STATES
!        ******************************************************************
                  
                  st = force(i,j,k,n) - str
                  
                  s_b(j+1)= s(i,j,k,n) + (HALF-dth*vadv(i,j+1,k))*slopey(i,j,k,1) + dth*st
                  s_t(j  )= s(i,j,k,n) - (HALF+dth*vadv(i,j,  k))*slopey(i,j,k,1) + dth*st
                  
               enddo
               
               
               do j = js, je+1 
                  sedgey(i,j,k,n)=merge(s_b(j),s_t(j),vadv(i,j,k).gt.ZERO)
                  savg = HALF*(s_t(j) + s_b(j))
                  sedgey(i,j,k,n)=merge(savg,sedgey(i,j,k,n),abs(vadv(i,j,k)) .lt. eps)
               enddo
               
               
               if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgey(i,js,k,n) = ZERO
                  else if (is_vel .and. n .ne. 2) then
                     sedgey(i,js,k,n) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
                  else
                     sedgey(i,js,k,n) = s_t(js)
                  endif
               elseif (phys_bc(2,1) .eq. INLET) then
                  sedgey(i,js,k,n) = s(i,js-1,k,n)
               elseif (phys_bc(2,1) .eq. OUTLET) then
                  sedgey(i,js,k,n) = s_t(js)
               endif
               
               if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 2) then
                     sedgey(i,je+1,k,n) = ZERO
                  else if (is_vel .and. n .ne. 2) then
                     sedgey(i,je+1,k,n) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
                  else
                     sedgey(i,je+1,k,n) = s_b(je+1)
                  endif
               elseif (phys_bc(2,2) .eq. INLET) then
                  sedgey(i,je+1,k,n) = s(i,je+1,k,n)
               elseif (phys_bc(2,2) .eq. OUTLET) then
                  sedgey(i,je+1,k,n) = s_b(je+1)
               endif
            enddo
         enddo
      enddo
      
!     ******************************************************************
!     ******************************************************************
!     ******************************************************************

!     loop for z fluxes
      do n = 1,ncomp
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
                  
                  sprgt = merge(s(ie+1,j,k,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  splft = merge(s(ie+1,j,k,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                  
                  if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 1) then
                        sprgt = ZERO
                        splft = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                        splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                     else
                        sprgt = splft
                     endif
                  endif
                  
                  splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                  savg  = HALF * (splft + sprgt)
                  splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                  
                  smrgt = s(i  ,j,k,n) - (HALF + dth*u(i  ,j,k,1)/hx)*slopex(i  ,j,k,n)
                  !    $            + dth * force(i  ,j,k,n)
                  smlft = s(i-1,j,k,n) + (HALF - dth*u(i-1,j,k,1)/hx)*slopex(i-1,j,k,n)
                  !    $            + dth * force(i-1,j,k,n)
                  
                  smrgt = merge(s(is-1,j,k,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  smlft = merge(s(is-1,j,k,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                  
                  if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 1) then
                        smrgt = ZERO
                        smlft = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                        smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                     else
                        smlft = smrgt
                     endif
                  endif
                  
                  sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smlft + smrgt)
                  sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                  
                  str = HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                  
!        ******************************************************************
!        MAKE TRANSVERSE DERIVATIVES IN Y-DIRECTION
!        ******************************************************************

                  spbot = s(i,j  ,k,n) + (HALF - dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
                  !    $            + dth * force(i,  j,k,n)
                  sptop = s(i,j+1,k,n) - (HALF + dth*u(i,j+1,k,2)/hy)*slopey(i,j+1,k,n)
                  !    $            + dth * force(i,j+1,k,n)
                  
                  sptop = merge(s(i,je+1,k,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  spbot = merge(s(i,je+1,k,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                  
                  if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 2) then
                        sptop = ZERO
                        spbot = ZERO
                     else if (is_vel .and. n .ne. 2) then
                        sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                        spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                     else
                        sptop = spbot
                     endif
                  endif
                  
                  splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                  savg  = HALF * (spbot + sptop)
                  splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                  
                  smtop = s(i,j  ,k,n) - (HALF + dth*u(i,j  ,k,2)/hy)*slopey(i,j  ,k,n)
                  !    $            + dth * force(i,j  ,k,n)
                  smbot = s(i,j-1,k,n) + (HALF - dth*u(i,j-1,k,2)/hy)*slopey(i,j-1,k,n)
                  !    $            + dth * force(i,j-1,k,n)
                  
                  smtop = merge(s(i,js-1,k,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  smbot = merge(s(i,js-1,k,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                  
                  if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                     if (is_vel .and. n .eq. 2) then
                        smtop = ZERO
                        smbot = ZERO
                     else if (is_vel .and. n .ne. 2) then
                        smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                        smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                     else
                        smbot = smtop
                     endif
                  endif
                  
                  sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                  savg   = HALF * (smbot + smtop)
                  sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                  
                  str =  str + HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy

!        ******************************************************************
!        MAKE DOWN AND UP STATES
!        ******************************************************************

                  st = force(i,j,k,n) - str
                  
                  s_d(k+1)= s(i,j,k,n) + (HALF-dth*wadv(i,j,k+1))*slopez(i,j,k,1) + dth*st
                  s_u(k  )= s(i,j,k,n) - (HALF+dth*wadv(i,j,k  ))*slopez(i,j,k,1) + dth*st

               enddo
               
               do k = ks, ke+1 
                  sedgez(i,j,k,n)=merge(s_d(k),s_u(k),wadv(i,j,k).gt.ZERO)
                  savg = HALF*(s_d(k) + s_u(k))
                  sedgez(i,j,k,n)=merge(savg,sedgez(i,j,k,n),abs(wadv(i,j,k)) .lt. eps)
               enddo
               
               if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 3) then
                     sedgez(i,j,ks,n) = ZERO
                  else if (is_vel .and. n .ne. 3) then
                     sedgez(i,j,ks,n) = merge(ZERO,s_u(ks),phys_bc(3,1).eq.NO_SLIP_WALL)
                  else
                     sedgez(i,j,ks,n) = s_u(ks)
                  endif
               elseif (phys_bc(3,1) .eq. INLET) then
                  sedgez(i,j,ks,n) = s(i,j,ks-1,n)
               elseif (phys_bc(3,1) .eq. OUTLET) then
                  sedgez(i,j,ks,n) = s_u(ks)
               endif
               
               if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 3) then
                     sedgez(i,j,ke+1,n) = ZERO
                  else if (is_vel .and. n .ne. 3) then
                     sedgez(i,j,ke+1,n) = merge(ZERO,s_d(ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
                  else
                     sedgez(i,j,ke+1,n) = s_d(ke+1)
                  endif
               elseif (phys_bc(3,2) .eq. INLET) then
                  sedgez(i,j,ke+1,n) = s(i,j,ke+1,n)
               elseif (phys_bc(3,2) .eq. OUTLET) then
                  sedgez(i,j,ke+1,n) = s_d(ke+1)
               endif
               
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

      end subroutine mkflux_macvel_3d

end module mkflux_macvel_module
