module mkflux_lowmemory_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module

  implicit none

contains

      subroutine mkflux_lowmemory_2d(s,u,sedgex,sedgey,umac,vmac, &
                                     force,divu,lo,dx,dt,is_vel, &
                                     phys_bc,adv_bc,ng,use_minion,is_conservative)

      integer, intent(in) :: lo(:),ng

      real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) :: adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel, use_minion, is_conservative(:)

      ! Local variables
      real(kind=dp_t), allocatable:: slopex(:,:,:)
      real(kind=dp_t), allocatable:: slopey(:,:,:)

      real(kind=dp_t) hx, hy, dt2, dt4, savg
      real(kind=dp_t) :: abs_eps, eps, umax

      integer :: hi(2)
      integer :: i,j,is,js,ie,je,n
      integer :: jc,jp
      integer :: slope_order = 4
      integer :: ncomp

      ! these correspond to s_L^x, etc.
      real(kind=dp_t), allocatable:: slx(:,:),srx(:,:),simhx(:,:)
      real(kind=dp_t), allocatable:: sly(:,:),sry(:,:),simhy(:,:)

      ! these correspond to \mathrm{sedge}_L^x, etc.
      real(kind=dp_t), allocatable:: sedgelx(:),sedgerx(:)
      real(kind=dp_t), allocatable:: sedgely(:),sedgery(:)
 
      ncomp = size(s,dim=3)

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))

      call slopex_2d(s,slopex,lo,ng,ncomp,adv_bc,slope_order)
      call slopey_2d(s,slopey,lo,ng,ncomp,adv_bc,slope_order)

      ! Note: All of these arrays are allocated to exactly the 
      ! size they need to be in order to compute edge states on 
      ! a domain with faces from lo(1):hi(1)+1,lo(2):hi(2)+1

      !***********************
      ! Normal predictor terms
      !***********************

      ! lo(1):hi(1)+1 in the x-direction
      ! 2 rows needed in y-direction
      allocate(slx  (lo(1):hi(1)+1,2))
      allocate(srx  (lo(1):hi(1)+1,2))
      allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

      ! lo(1)-1:hi(1)+1 in the x-direction
      ! 2 rows needed in y-direction
      allocate(sly  (lo(1)-1:hi(1)+1,2))
      allocate(sry  (lo(1)-1:hi(1)+1,2))
      allocate(simhy(lo(1)-1:hi(1)+1,2))

      !************
      ! Edge states
      !************

      ! lo(1):hi(1)+1 in x-direction
      allocate(sedgelx(lo(1):hi(1)+1))
      allocate(sedgerx(lo(1):hi(1)+1))

      ! lo(1):hi(1) in x-direction
      allocate(sedgely(lo(1):hi(1)))
      allocate(sedgery(lo(1):hi(1)))

      abs_eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)

      dt2 = HALF*dt
      dt4 = dt/4.0d0

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
      if(umax .eq. 0.d0) then
         eps = abs_eps
      else
         eps = abs_eps * umax
      endif

      !*************************************
      ! Pseudo code
      !*************************************
      !
      !  do j=js-1,je+1
      !     1. Compute simhx(is:ie+1,j)
      !     if(j .gt. js-1) then
      !        2. Compute simhy(is-1:ie+1,j)
      !        3. Compute sedgey(is:ie,j)
      !     endif
      !     if(j .gt. js) then
      !        4. Compute sedgex(is:ie+1,j-1)
      !     endif
      !     5. Cycle indeces
      !  enddo
      !
      !*************************************
      ! End pseudo code
      !*************************************

      ! loop over components
      do n=1,ncomp

         jc = 1
         jp = 2

         do j=js-1,je+1

!******************************************************************
! 1. Compute simhx(is:ie+1,j)
!******************************************************************

            do i=is,ie+1
               ! make slx, srx with 1D extrapolation
               slx(i,jc) = s(i-1,j,n) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,n)
               srx(i,jc) = s(i  ,j,n) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,n)

               ! add s*u_x term where u_x = divu - v_y
               if(is_conservative(n)) then
                  slx(i,jc) = slx(i,jc) - dt2*s(i-1,j,n)*(divu(i-1,j) &
                       - (vmac(i-1,j+1)-vmac(i-1,j))/hy)
                  srx(i,jc) = srx(i,jc) - dt2*s(i  ,j,n)*(divu(i  ,j) &
                       - (vmac(i  ,j+1)-vmac(i  ,j))/hy)
               endif

               ! add source terms
               if(use_minion) then
                  slx(i,jc) = slx(i,jc) + dt2*force(i-1,j,n)
                  srx(i,jc) = srx(i,jc) + dt2*force(i  ,j,n)
               endif

               ! impose lo side bc's
               if(i .eq. is) then
                  slx(i,jc) = merge(s(is-1,j,n),slx(i,jc),phys_bc(1,1) .eq. INLET)
                  srx(i,jc) = merge(s(is-1,j,n),srx(i,jc),phys_bc(1,1) .eq. INLET)
                  if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                     if(is_vel .and. n .eq. 1) then
                        slx(i,jc) = ZERO
                        srx(i,jc) = ZERO
                     else if(is_vel .and. n .ne. 1) then
                        slx(i,jc) = merge(ZERO,srx(i,jc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        srx(i,jc) = merge(ZERO,srx(i,jc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                     else
                        slx(i,jc) = srx(i,jc)
                     endif
                  endif
               endif
               
               ! impose hi side bc's
               if(i .eq. ie+1) then
                  slx(i,jc) = merge(s(ie+1,j,n),slx(i,jc),phys_bc(1,2) .eq. INLET)
                  srx(i,jc) = merge(s(ie+1,j,n),srx(i,jc),phys_bc(1,2) .eq. INLET)
                  if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                     if (is_vel .and. n .eq. 1) then
                        slx(i,jc) = ZERO
                        srx(i,jc) = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        slx(i,jc) = merge(ZERO,slx(i,jc),phys_bc(1,2).eq.NO_SLIP_WALL)
                        srx(i,jc) = merge(ZERO,slx(i,jc),phys_bc(1,2).eq.NO_SLIP_WALL)
                     else
                        srx(i,jc) = slx(i,jc)
                     endif
                  endif
               endif
               
               ! make simhx by solving Riemann problem
               simhx(i,jc) = merge(slx(i,jc),srx(i,jc),umac(i,j) .gt. ZERO)
               savg = HALF*(slx(i,jc)+srx(i,jc))
               simhx(i,jc) = merge(simhx(i,jc),savg,abs(umac(i,j)) .gt. eps)
            enddo

            if(j .gt. js-1) then

!******************************************************************
! 2. Compute simhy(is-1:ie+1,j)
!******************************************************************

               do i=is-1,ie+1
                  ! make sly, sry with 1D extrapolation
                  sly(i,jc) = s(i,j-1,n) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,n)
                  sry(i,jc) = s(i,j  ,n) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,n)

                  ! add s*v_y term where v_y = divu - u_x
                  if(is_conservative(n)) then
                     sly(i,jc) = sly(i,jc) - dt2*s(i,j-1,n)*(divu(i,j-1) &
                          - (umac(i+1,j-1)-umac(i,j-1))/hx)
                     sry(i,jc) = sry(i,jc) - dt2*s(i,j  ,n)*(divu(i,j  ) &
                          - (umac(i+1,j  )-umac(i,j  ))/hx)
                  endif
               
                  ! add source terms
                  if(use_minion) then
                     sly(i,jc) = sly(i,jc) + dt2*force(1,j-1,n)
                     sry(i,jc) = sry(i,jc) + dt2*force(i,j  ,n)
                  endif
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     sly(i,jc) = merge(s(is,j-1,n),sly(i,jc),phys_bc(2,1) .eq. INLET)
                     sry(i,jc) = merge(s(is,j-1,n),sry(i,jc),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           sly(i,jc) = ZERO
                           sry(i,jc) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           sly(i,jc) = merge(ZERO,sry(i,jc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sry(i,jc) = merge(ZERO,sry(i,jc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           sly(i,jc) = sry(i,jc)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     sly(i,jc) = merge(s(i,je+1,n),sly(i,jc),phys_bc(2,2) .eq. INLET)
                     sry(i,jc) = merge(s(i,je+1,n),sry(i,jc),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           sly(i,jc) = ZERO
                           sry(i,jc) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           sly(i,jc) = merge(ZERO,sly(i,jc),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sry(i,jc) = merge(ZERO,sly(i,jc),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sry(i,jc) = sly(i,jc)
                        endif
                     endif
                  endif
                  
                  ! make simhy by solving Riemann problem
                  simhy(i,jc) = merge(sly(i,jc),sry(i,jc),vmac(i,j) .gt. ZERO)
                  savg = HALF*(sly(i,jc)+sry(i,jc))
                  simhy(i,jc) = merge(simhy(i,jc),savg,abs(vmac(i,j)) .gt. eps)
               enddo

!******************************************************************
! 3. Compute sedgey(is:ie,j)
!******************************************************************
               
               do i=is,ie
                  ! make sedgely, sedgery
                  if(is_conservative(n)) then
                     sedgely(i) = sly(i,jc) &
                          - (dt2/hx)*(simhx(i+1,jp)*umac(i+1,j-1) - simhx(i,jp)*umac(i,j-1))
                     sedgery(i) = sry(i,jc) &
                          - (dt2/hx)*(simhx(i+1,jc)*umac(i+1,j  ) - simhx(i,jc)*umac(i,j  ))
                  else
                     sedgely(i) = sly(i,jc) &
                          - (dt4/hx)*(umac(i+1,j-1)+umac(i,j-1))*(simhx(i+1,jp)-simhx(i,jp))
                     sedgery(i) = sry(i,jc) &
                          - (dt4/hx)*(umac(i+1,j  )+umac(i,j  ))*(simhx(i+1,jc)-simhx(i,jc))
                  endif
                  
                  ! if use_minion is true, we have already accounted for source terms
                  ! in sly and sry; otherwise, we need to account for them here.
                  if(.not. use_minion) then
                     sedgely(i) = sedgely(i) + dt2*force(i,j-1,n)
                     sedgery(i) = sedgery(i) + dt2*force(i,j  ,n)
                  endif
                  
                  ! make sedgey by solving Riemann problem
                  ! boundary conditions enforced outside of i,j loop
                  sedgey(i,j,n) = merge(sedgely(i),sedgery(i),vmac(i,j) .gt. ZERO)
                  savg = HALF*(sedgely(i)+sedgery(i))
                  sedgey(i,j,n) = merge(sedgey(i,j,n),savg,abs(vmac(i,j)) .gt. eps)
         
                  ! sedgey boundary conditions
                  if(j .eq. js) then
                     ! lo side
                     if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           sedgey(i,js,n) = ZERO
                        elseif (is_vel .and. n .ne. 2) then
                           sedgey(i,js,n) = merge(ZERO,sedgery(i),phys_bc(2,1).eq.NO_SLIP_WALL)
                        else 
                           sedgey(i,js,n) = sedgery(i)
                        endif
                     elseif (phys_bc(2,1) .eq. INLET) then
                        sedgey(i,js,n) = s(i,js-1,n)
                     elseif (phys_bc(2,1) .eq. OUTLET) then
                        if (is_vel .and. n.eq.2) then
                           sedgey(i,js,n) = MIN(sedgery(i),ZERO)
                        else
                           sedgey(i,js,n) = sedgery(i)
                        end if
                     endif
                  else if(j .eq. je+1) then
                     ! hi side
                     if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           sedgey(i,je+1,n) = ZERO
                        elseif (is_vel .and. n .ne. 2) then
                           sedgey(i,je+1,n) = merge(ZERO,sedgely(i),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else 
                           sedgey(i,je+1,n) = sedgely(i)
                        endif
                     elseif (phys_bc(2,2) .eq. INLET) then
                        sedgey(i,je+1,n) = s(i,je+1,n)
                     elseif (phys_bc(2,2) .eq. OUTLET) then
                        if (is_vel .and. n.eq.2) then
                           sedgey(i,je+1,n) = MAX(sedgely(i),ZERO)
                        else
                           sedgey(i,je+1,n) = sedgely(i)
                        end if
                     endif
                  endif
               enddo

            endif ! end if(j .gt. js-1)

            if(j .gt. js) then

!******************************************************************
! 4. Compute sedgex(is:ie+1,j-1)
!******************************************************************

               do i=is,ie+1
                  ! make sedgelx, sedgerx
                  if(is_conservative(n)) then
                     sedgelx(i) = slx(i,jp) &
                          - (dt2/hy)*(simhy(i-1,jc)*vmac(i-1,j) - simhy(i-1,jp)*vmac(i-1,j-1))
                     sedgerx(i) = srx(i,jp) &
                          - (dt2/hy)*(simhy(i,  jc)*vmac(i,  j) - simhy(i,  jp)*vmac(i,  j-1))
                  else
                     sedgelx(i) = slx(i,jp) &
                          - (dt4/hy)*(vmac(i-1,j)+vmac(i-1,j-1))*(simhy(i-1,jc)-simhy(i-1,jp))
                     sedgerx(i) = srx(i,jp) &
                          - (dt4/hy)*(vmac(i,  j)+vmac(i,  j-1))*(simhy(i,  jc)-simhy(i,  jp))
                  endif
                  
                  ! if use_minion is true, we have already accounted for source terms
                  ! in slx and srx; otherwise, we need to account for them here.
                  if(.not. use_minion) then
                     sedgelx(i) = sedgelx(i) + dt2*force(i-1,j-1,n)
                     sedgerx(i) = sedgerx(i) + dt2*force(i  ,j-1,n)
                  endif
                  
                  ! make sedgex by solving Riemann problem
                  ! boundary conditions enforced outside of i,j loop
                  sedgex(i,j-1,n) = merge(sedgelx(i),sedgerx(i),umac(i,j-1) .gt. ZERO)
                  savg = HALF*(sedgelx(i)+sedgerx(i))
                  sedgex(i,j-1,n) = merge(sedgex(i,j-1,n),savg,abs(umac(i,j-1)) .gt. eps)
               enddo
               
               ! sedgex boundary conditions
               ! lo side
               if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(is,j-1,n) = ZERO
                  elseif (is_vel .and. n .ne. 1) then
                     sedgex(is,j-1,n) = merge(ZERO,sedgerx(is),phys_bc(1,1).eq.NO_SLIP_WALL)
                  else
                     sedgex(is,j-1,n) = sedgerx(is)
                  endif
               elseif (phys_bc(1,1) .eq. INLET) then
                  sedgex(is,j-1,n) = s(is-1,j-1,n)
               elseif (phys_bc(1,1) .eq. OUTLET) then
                  if (is_vel .and. n.eq.1) then
                     sedgex(is,j-1,n) = MIN(sedgerx(is),ZERO)
                  else
                     sedgex(is,j-1,n) = sedgerx(is)
                  end if
               endif
               
               ! hi side
               if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                  if (is_vel .and. n .eq. 1) then
                     sedgex(ie+1,j-1,n) = ZERO
                  else if (is_vel .and. n .ne. 1) then
                     sedgex(ie+1,j-1,n) = merge(ZERO,sedgelx(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
                  else 
                     sedgex(ie+1,j-1,n) = sedgelx(ie+1)
                  endif
               elseif (phys_bc(1,2) .eq. INLET) then
                  sedgex(ie+1,j-1,n) = s(ie+1,j-1,n)
               elseif (phys_bc(1,2) .eq. OUTLET) then
                  if (is_vel .and. n.eq.1) then
                     sedgex(ie+1,j-1,n) = MAX(sedgelx(ie+1),ZERO)
                  else
                     sedgex(ie+1,j-1,n) = sedgelx(ie+1)
                  end if
               endif

            endif ! end if(j .gt. js)
            
!******************************************************************
! 5. Cycle indeces
!******************************************************************

            jc = 3 - jc
            jp = 3 - jp

         enddo ! end loop over j
      enddo ! end loop over components

      deallocate(slopex)
      deallocate(slopey)

      deallocate(slx)
      deallocate(srx)
      deallocate(sly)
      deallocate(sry)

      deallocate(simhx)
      deallocate(simhy)

      deallocate(sedgelx)
      deallocate(sedgerx)
      deallocate(sedgely)
      deallocate(sedgery)

      end subroutine mkflux_lowmemory_2d

      subroutine mkflux_lowmemory_3d(s,u,sedgex,sedgey,sedgez,&
                                     umac,vmac,wmac, &
                                     force,divu,lo,dx,dt,is_vel, &
                                     phys_bc,adv_bc,ng,use_minion,is_conservative)

      integer, intent(in) :: lo(:),ng

      real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
      real(kind=dp_t),intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t),intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
      real(kind=dp_t),intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)
      real(kind=dp_t),intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) :: adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel, use_minion, is_conservative(:)

      ! Local variables
      real(kind=dp_t), allocatable:: slopex(:,:,:,:)
      real(kind=dp_t), allocatable:: slopey(:,:,:,:)
      real(kind=dp_t), allocatable:: slopez(:,:,:,:)

      real(kind=dp_t) hx, hy, hz, dt2, dt3, dt4, dt6, savg
      real(kind=dp_t) :: abs_eps, eps, umax

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke,n
      integer :: kc,kp
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
      real(kind=dp_t), allocatable:: sedgelx(:,:),sedgerx(:,:)
      real(kind=dp_t), allocatable:: sedgely(:,:),sedgery(:,:)
      real(kind=dp_t), allocatable:: sedgelz(:,:),sedgerz(:,:)

      ncomp = size(s,dim=4)

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))
      allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))

      do k = lo(3)-1,hi(3)+1
         call slopex_2d(s(:,:,k,:),slopex(:,:,k,:),lo,ng,ncomp,adv_bc,slope_order)
         call slopey_2d(s(:,:,k,:),slopey(:,:,k,:),lo,ng,ncomp,adv_bc,slope_order)
      end do
      call slopez_3d(s,slopez,lo,ng,ncomp,adv_bc,slope_order)

      ! Note: All of these arrays are allocated to exactly the 
      ! size they need to be in order to compute edge states on 
      ! a domain with faces from lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1

      !***********************
      ! Normal predictor terms
      !***********************

      ! lo(1):hi(1)+1 in the x-direction
      ! lo(2)-1:hi(2)+1 in the y-direction
      ! 2 rows needed in z-direction
      allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

      ! lo(1)-1:hi(1)+1 in the x-direction
      ! lo(2):hi(2)+1 in the y-direction
      ! 2 rows needed in z-direction
      allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
      allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
      allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

      ! lo(1)-1:hi(1)+1 in the x-direction
      ! lo(2)-1:hi(2)+1 in the y-direction
      ! 2 rows needed in z-direction
      allocate(slz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(srz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(simhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

      !*****************
      ! Transverse terms
      !*****************

      ! lo(1):hi(1)+1 in the x-direction
      ! lo(2):hi(2) in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slxy  (lo(1):hi(1)+1,lo(2):hi(2),2))
      allocate(srxy  (lo(1):hi(1)+1,lo(2):hi(2),2))
      allocate(simhxy(lo(1):hi(1)+1,lo(2):hi(2),2))

      ! lo(1):hi(1)+1 in the x-direction
      ! lo(2)-1:hi(2)+1 in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(srxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(simhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

      ! lo(1):hi(1) in the x-direction
      ! lo(2):hi(2)+1 in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slyx  (lo(1):hi(1),lo(2):hi(2)+1,2))
      allocate(sryx  (lo(1):hi(1),lo(2):hi(2)+1,2))
      allocate(simhyx(lo(1):hi(1),lo(2):hi(2)+1,2))

      ! lo(1)-1:hi(1)+1 in the x-direction
      ! lo(2):hi(2)+1 in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slyz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
      allocate(sryz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
      allocate(simhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

      ! lo(1):hi(1) in the x-direction
      ! lo(2)-1:hi(2)+1 in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,2))
      allocate(srzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,2))
      allocate(simhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,2))

      ! lo(1)-1:hi(1)+1 in the x-direction
      ! lo(2):hi(2) in the y-direction
      ! 2 rows needed in the z-direction
      allocate(slzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),2))
      allocate(srzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),2))
      allocate(simhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),2))

      !************
      ! Edge states
      !************

      ! lo(1):hi(1)+1 in x-direction
      ! lo(2):hi(2) in y-direction
      allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2)))
      allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2)))

      ! lo(1):hi(1) in x-direction
      ! lo(2):hi(2)+1 in y-direction
      allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1))
      allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1))

      ! lo(1):hi(1) in x-direction
      ! lo(2):hi(2) in y-direction
      allocate(sedgelz(lo(1):hi(1),lo(2):hi(2)))
      allocate(sedgerz(lo(1):hi(1),lo(2):hi(2)))

      abs_eps = 1.0e-8

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

      dt2 = HALF*dt
      dt3 = dt/3.0d0
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
      if(umax .eq. 0.d0) then
         eps = abs_eps
      else
         eps = abs_eps * umax
      endif

      !*************************************
      ! Pseudo code
      !*************************************
      !
      !  do j=ks-1,ke+1
      !     1. Compute simhx (is  :ie+1,js-1:je+1,k) 
      !     2. Compute simhy (is-1:ie+1,js  :je+1,k)
      !     3. Compute simhxy(is  :ie+1,js  :je  ,k)
      !     4. Compute simhyx(is  :ie  ,js  :je+1,k)
      !     if(k .gt. ks-1) then
      !        5. Compute simhz (is-1:ie+1,js-1:je+1,k)
      !        6. Compute sedgez(is  :ie  ,js  :je  ,k)
      !        7. Compute simhzx(is  :ie  ,js-1:je+1,k)
      !        8. Compute simhzy(is-1:ie+1,js  :je  ,k)
      !     endif
      !     if(k .gt. ks) then
      !        9. Compute simhxz(is  :ie+1,js-1:je+1,k-1)
      !        10.Compute simhyz(is-1:ie+1,js  :je+1,k-1)
      !        11.Compute sedgex(is  :ie+1,js  :je,  k-1)
      !        12.Compute sedgey(is  :ie  ,js  :je+1,k-1)
      !     endif
      !     13. Cycle indeces
      !  enddo
      !
      !*************************************
      ! End pseudo code
      !*************************************

      ! loop over components
      do n=1,ncomp

         kc = 1
         kp = 2

         do k=ks-1,ke+1

!******************************************************************
! 1. Compute simhx (is  :ie+1,js-1:je+1,k)
!******************************************************************

            do j=js-1,je+1
               do i=is,ie+1
                  ! make slx, srx with 1D extrapolation
                  slx(i,j,kc) = s(i-1,j,k,n) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,n)
                  srx(i,j,kc) = s(i  ,j,k,n) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i,  j,k,n)
                  
                  ! add s*u_x term where u_x = divu - v_y - w_z
                  if(is_conservative(n)) then
                     slx(i,j,kc) = slx(i,j,kc) - dt2*s(i-1,j,k,n)*(divu(i-1,j,k) &
                          - (vmac(i-1,j+1,k)-vmac(i-1,j,k))/hy &
                          - (wmac(i-1,j,k+1)-wmac(i-1,j,k))/hz)
                     srx(i,j,kc) = srx(i,j,kc) - dt2*s(i  ,j,k,n)*(divu(i  ,j,k) &
                          - (vmac(i  ,j+1,k)-vmac(i  ,j,k))/hy &
                          - (wmac(i  ,j,k+1)-wmac(i  ,j,k))/hz)
                  endif
                  
                  ! add source terms
                  if(use_minion) then
                     slx(i,j,kc) = slx(i,j,kc) + dt2*force(i-1,j,k,n)
                     srx(i,j,kc) = srx(i,j,kc) + dt2*force(i  ,j,k,n)
                  endif
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slx(i,j,kc) = merge(s(is-1,j,k,n),slx(i,j,kc),phys_bc(1,1) .eq. INLET)
                     srx(i,j,kc) = merge(s(is-1,j,k,n),srx(i,j,kc),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 1) then
                           slx(i,j,kc) = ZERO
                           srx(i,j,kc) = ZERO
                        else if(is_vel .and. n .ne. 1) then
                           slx(i,j,kc) = merge(ZERO,srx(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srx(i,j,kc) = merge(ZERO,srx(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        else
                           slx(i,j,kc) = srx(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slx(i,j,kc) = merge(s(ie+1,j,k,n),slx(i,j,kc),phys_bc(1,2) .eq. INLET)
                     srx(i,j,kc) = merge(s(ie+1,j,k,n),srx(i,j,kc),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 1) then
                           slx(i,j,kc) = ZERO
                           srx(i,j,kc) = ZERO
                        else if (is_vel .and. n .ne. 1) then
                           slx(i,j,kc) = merge(ZERO,slx(i,j,kc),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srx(i,j,kc) = merge(ZERO,slx(i,j,kc),phys_bc(1,2).eq.NO_SLIP_WALL)
                        else
                           srx(i,j,kc) = slx(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! make simhx by solving Riemann problem
                  simhx(i,j,kc) = merge(slx(i,j,kc),srx(i,j,kc),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(slx(i,j,kc)+srx(i,j,kc))
                  simhx(i,j,kc) = merge(simhx(i,j,kc),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo

!******************************************************************
! 2. Compute simhy (is-1:ie+1,js  :je+1,k)
!******************************************************************

            do j=js,je+1
               do i=is-1,ie+1
                  ! make sly, sry with 1D extrapolation
                  sly(i,j,kc) = s(i,j-1,k,n) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,n)
                  sry(i,j,kc) = s(i,j,  k,n) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j,  k,n)
                  
                  ! add s*v_y term where v_y = divu - u_x - w_z
                  if(is_conservative(n)) then
                     sly(i,j,kc) = sly(i,j,kc) - dt2*s(i,j-1,k,n)*(divu(i,j-1,k) &
                          - (umac(i+1,j-1,k)-umac(i,j-1,k))/hx &
                          - (wmac(i,j-1,k+1)-wmac(i,j-1,k))/hz)
                     sry(i,j,kc) = sry(i,j,kc) - dt2*s(i,j  ,k,n)*(divu(i,j  ,k) &
                          - (umac(i+1,j  ,k)-umac(i,j  ,k))/hx &
                          - (wmac(i,j  ,k+1)-wmac(i,j  ,k))/hz) 
                  endif

                  ! add source terms
                  if(use_minion) then
                     sly(i,j,kc) = sly(i,j,kc) + dt2*force(1,j-1,k,n)
                     sry(i,j,kc) = sry(i,j,kc) + dt2*force(i,j  ,k,n)
                  endif

                  ! impose lo side bc's
                  if(j .eq. js) then
                     sly(i,j,kc) = merge(s(is,j-1,k,n),sly(i,j,kc),phys_bc(2,1) .eq. INLET)
                     sry(i,j,kc) = merge(s(is,j-1,k,n),sry(i,j,kc),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           sly(i,j,kc) = ZERO
                           sry(i,j,kc) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           sly(i,j,kc) = merge(ZERO,sry(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sry(i,j,kc) = merge(ZERO,sry(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           sly(i,j,kc) = sry(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     sly(i,j,kc) = merge(s(i,je+1,k,n),sly(i,j,kc),phys_bc(2,2) .eq. INLET)
                     sry(i,j,kc) = merge(s(i,je+1,k,n),sry(i,j,kc),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           sly(i,j,kc) = ZERO
                           sry(i,j,kc) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           sly(i,j,kc) = merge(ZERO,sly(i,j,kc),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sry(i,j,kc) = merge(ZERO,sly(i,j,kc),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sry(i,j,kc) = sly(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! make simhy by solving Riemann problem
                  simhy(i,j,kc) = merge(sly(i,j,kc),sry(i,j,kc),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(sly(i,j,kc)+sry(i,j,kc))
                  simhy(i,j,kc) = merge(simhy(i,j,kc),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo

!******************************************************************
! 3. Compute simhxy(is  :ie+1,js  :je  ,k)
!******************************************************************

            do j=js,je
               do i=is,ie+1
                  ! make slxy, srxy by updating 1D extrapolation
                  if(is_conservative(n)) then
                     slxy(i,j,kc) = slx(i,j,kc) - (dt3/hy)*(simhy(i-1,j+1,kc)*vmac(i-1,j+1,k) - simhy(i-1,j,kc)*vmac(i-1,j,k))
                     srxy(i,j,kc) = srx(i,j,kc) - (dt3/hy)*(simhy(i,  j+1,kc)*vmac(i,  j+1,k) - simhy(i,  j,kc)*vmac(i,  j,k))
                  else
                     slxy(i,j,kc) = slx(i,j,kc) - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhy(i-1,j+1,kc)-simhy(i-1,j,kc))
                     srxy(i,j,kc) = srx(i,j,kc) - (dt6/hy)*(vmac(i,  j+1,k)+vmac(i,  j,k))*(simhy(i,  j+1,kc)-simhy(i,  j,kc))
                  endif
                  
                  ! impose lo side bc's
                  if(i .eq. is) then
                     slxy(i,j,kc) = merge(s(is-1,j,k,n),slxy(i,j,kc),phys_bc(1,1) .eq. INLET)
                     srxy(i,j,kc) = merge(s(is-1,j,k,n),srxy(i,j,kc),phys_bc(1,1) .eq. INLET)
                     if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 1) then
                           slxy(i,j,kc) = ZERO
                           srxy(i,j,kc) = ZERO
                        else if(is_vel .and. n .ne. 1) then
                           slxy(i,j,kc) = merge(ZERO,srxy(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           srxy(i,j,kc) = merge(ZERO,srxy(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                        else
                           slxy(i,j,kc) = srxy(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(i .eq. ie+1) then
                     slxy(i,j,kc) = merge(s(ie+1,j,k,n),slxy(i,j,kc),phys_bc(1,2) .eq. INLET)
                     srxy(i,j,kc) = merge(s(ie+1,j,k,n),srxy(i,j,kc),phys_bc(1,2) .eq. INLET)
                     if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 1) then
                           slxy(i,j,kc) = ZERO
                           srxy(i,j,kc) = ZERO
                        else if (is_vel .and. n .ne. 1) then
                           slxy(i,j,kc) = merge(ZERO,slxy(i,j,kc),phys_bc(1,2).eq.NO_SLIP_WALL)
                           srxy(i,j,kc) = merge(ZERO,slxy(i,j,kc),phys_bc(1,2).eq.NO_SLIP_WALL)
                        else
                           srxy(i,j,kc) = slxy(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! make simhxy by solving Riemann problem
                  simhxy(i,j,kc) = merge(slxy(i,j,kc),srxy(i,j,kc),umac(i,j,k) .gt. ZERO)
                  savg = HALF*(slxy(i,j,kc)+srxy(i,j,kc))
                  simhxy(i,j,kc) = merge(simhxy(i,j,kc),savg,abs(umac(i,j,k)) .gt. eps)
               enddo
            enddo

!******************************************************************
! 4. Compute simhyx(is  :ie  ,js  :je+1,k)
!******************************************************************

            do j=js,je+1
               do i=is,ie
                  ! make slyx, sryx by updating 1D extrapolation
                  if(is_conservative(n)) then
                     slyx(i,j,kc) = sly(i,j,kc) - (dt3/hx)*(simhx(i+1,j-1,kc)*umac(i+1,j-1,k) - simhx(i,j-1,kc)*umac(i,j-1,k))
                     sryx(i,j,kc) = sry(i,j,kc) - (dt3/hx)*(simhx(i+1,j,  kc)*umac(i+1,j,  k) - simhx(i,j,  kc)*umac(i,j,  k))
                  else
                     slyx(i,j,kc) = sly(i,j,kc) - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhx(i+1,j-1,kc)-simhx(i,j-1,kc))
                     sryx(i,j,kc) = sry(i,j,kc) - (dt6/hx)*(umac(i+1,j,  k)+umac(i,j,  k))*(simhx(i+1,j,  kc)-simhx(i,j,  kc))
                  endif
                  
                  ! impose lo side bc's
                  if(j .eq. js) then
                     slyx(i,j,kc) = merge(s(is,j-1,k,n),slyx(i,j,kc),phys_bc(2,1) .eq. INLET)
                     sryx(i,j,kc) = merge(s(is,j-1,k,n),sryx(i,j,kc),phys_bc(2,1) .eq. INLET)
                     if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                        if(is_vel .and. n .eq. 2) then
                           slyx(i,j,kc) = ZERO
                           sryx(i,j,kc) = ZERO
                        else if(is_vel .and. n .ne. 2) then
                           slyx(i,j,kc) = merge(ZERO,sryx(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           sryx(i,j,kc) = merge(ZERO,sryx(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                        else
                           slyx(i,j,kc) = sryx(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! impose hi side bc's
                  if(j .eq. je+1) then
                     slyx(i,j,kc) = merge(s(i,je+1,k,n),slyx(i,j,kc),phys_bc(2,2) .eq. INLET)
                     sryx(i,j,kc) = merge(s(i,je+1,k,n),sryx(i,j,kc),phys_bc(2,2) .eq. INLET)
                     if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                        if (is_vel .and. n .eq. 2) then
                           slyx(i,j,kc) = ZERO
                           sryx(i,j,kc) = ZERO
                        else if (is_vel .and. n .ne. 2) then
                           slyx(i,j,kc) = merge(ZERO,slyx(i,j,kc),phys_bc(2,2).eq.NO_SLIP_WALL)
                           sryx(i,j,kc) = merge(ZERO,slyx(i,j,kc),phys_bc(2,2).eq.NO_SLIP_WALL)
                        else
                           sryx(i,j,kc) = slyx(i,j,kc)
                        endif
                     endif
                  endif
                  
                  ! make simhyx by solving Riemann problem
                  simhyx(i,j,kc) = merge(slyx(i,j,kc),sryx(i,j,kc),vmac(i,j,k) .gt. ZERO)
                  savg = HALF*(slyx(i,j,kc)+sryx(i,j,kc))
                  simhyx(i,j,kc) = merge(simhyx(i,j,kc),savg,abs(vmac(i,j,k)) .gt. eps)
               enddo
            enddo

            if(k .gt. ks-1) then

!******************************************************************
! 5. Compute simhz (is-1:ie+1,js-1:je+1,k)
!******************************************************************

               do j=js-1,je+1
                  do i=is-1,ie+1
                     ! make slz, srz with 1D extrapolation
                     slz(i,j,kc) = s(i,j,k-1,n) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,n)
                     srz(i,j,kc) = s(i,j,k,  n) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k,  n)
                     
                     ! add s*w_z term where w_z = divu - u_x - w_y
                     if(is_conservative(n)) then
                        slz(i,j,kc) = slz(i,j,kc) - dt2*s(i,j,k-1,n)*(divu(i,j,k-1) &
                             - (umac(i+1,j,k-1)-umac(i,j,k-1))/hx &
                             - (vmac(i,j+1,k-1)-vmac(i,j,k-1))/hy)
                        srz(i,j,kc) = srz(i,j,kc) - dt2*s(i,j,k,n)*(divu(i,j,k) &
                             - (umac(i+1,j,k)-umac(i,j,k))/hx &
                             - (vmac(i,j+1,k)-vmac(i,j,k))/hy)     
                     endif
                     
                     ! add source terms
                     if(use_minion) then
                        slz(i,j,kc) = slz(i,j,kc) + dt2*force(i,j,k-1,n)
                        srz(i,j,kc) = srz(i,j,kc) + dt2*force(i,j,k,n)
                     endif
                     
                     ! impose lo side bc's
                     if(k .eq. ks) then
                        slz(i,j,kc) = merge(s(is,j,k-1,n),slz(i,j,kc),phys_bc(3,1) .eq. INLET)
                        srz(i,j,kc) = merge(s(is,j,k-1,n),srz(i,j,kc),phys_bc(3,1) .eq. INLET)
                        if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                           if(is_vel .and. n .eq. 3) then
                              slz(i,j,kc) = ZERO
                              srz(i,j,kc) = ZERO
                           else if(is_vel .and. n .ne. 3) then
                              slz(i,j,kc) = merge(ZERO,srz(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                              srz(i,j,kc) = merge(ZERO,srz(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           else
                              slz(i,j,kc) = srz(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! impose hi side bc's
                     if(k .eq. ke+1) then
                        slz(i,j,kc) = merge(s(i,j,ke+1,n),slz(i,j,kc),phys_bc(3,2) .eq. INLET)
                        srz(i,j,kc) = merge(s(i,j,ke+1,n),srz(i,j,kc),phys_bc(3,2) .eq. INLET)
                        if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 3) then
                              slz(i,j,kc) = ZERO
                              srz(i,j,kc) = ZERO
                           else if (is_vel .and. n .ne. 3) then
                              slz(i,j,kc) = merge(ZERO,slz(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                              srz(i,j,kc) = merge(ZERO,slz(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                           else
                              srz(i,j,kc) = slz(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! make simhz by solving Riemann problem
                     simhz(i,j,kc) = merge(slz(i,j,kc),srz(i,j,kc),wmac(i,j,k) .gt. ZERO)
                     savg = HALF*(slz(i,j,kc)+srz(i,j,kc))
                     simhz(i,j,kc) = merge(simhz(i,j,kc),savg,abs(wmac(i,j,k)) .gt. eps)
                  enddo
               enddo

!******************************************************************
! 6. Compute sedgez(is  :ie  ,js  :je  ,k)
!******************************************************************
               
               do j=js,je
                  do i=is,ie
                     ! make sedgelz, sedgerz
                     if(is_conservative(n)) then
                        sedgelz(i,j) = slz(i,j,kc) &
                             - (dt2/hx)*(simhxy(i+1,j,kp)*umac(i+1,j,k-1) - simhxy(i,j,kp)*umac(i,j,k-1)) &
                             - (dt2/hy)*(simhyx(i,j+1,kp)*vmac(i,j+1,k-1) - simhyx(i,j,kp)*vmac(i,j,k-1))
                        sedgerz(i,j) = srz(i,j,kc) &
                             - (dt2/hx)*(simhxy(i+1,j,kc)*umac(i+1,j,k  ) - simhxy(i,j,kc)*umac(i,j,k  )) &
                             - (dt2/hy)*(simhyx(i,j+1,kc)*vmac(i,j+1,k  ) - simhyx(i,j,kc)*vmac(i,j,k  ))
                     else
                        sedgelz(i,j) = slz(i,j,kc) &
                             - (dt4/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhxy(i+1,j,kp)-simhxy(i,j,kp)) &
                             - (dt4/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhyx(i,j+1,kp)-simhyx(i,j,kp))
                        sedgerz(i,j) = srz(i,j,kc) &
                             - (dt4/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhxy(i+1,j,kc)-simhxy(i,j,kc)) &
                             - (dt4/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhyx(i,j+1,kc)-simhyx(i,j,kc))
                     endif
                     
                     ! if use_minion is true, we have already accounted for source terms
                     ! in slz and srz; otherwise, we need to account for them here.
                     if(.not. use_minion) then
                        sedgelz(i,j) = sedgelz(i,j) + dt2*force(i,j,k-1,n)
                        sedgerz(i,j) = sedgerz(i,j) + dt2*force(i,j,k  ,n)
                     endif
                     
                     ! make sedgez by solving Riemann problem
                     ! boundary conditions enforced outside of i,j,k loop
                     sedgez(i,j,k,n) = merge(sedgelz(i,j),sedgerz(i,j),wmac(i,j,k) .gt. ZERO)
                     savg = HALF*(sedgelz(i,j)+sedgerz(i,j))
                     sedgez(i,j,k,n) = merge(sedgez(i,j,k,n),savg,abs(wmac(i,j,k)) .gt. eps)
                     
                     ! sedgez boundary conditions
                     if(k .eq. ks) then
                        ! lo side
                        if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 2) then
                              sedgez(i,j,ks,n) = ZERO
                           elseif (is_vel .and. n .ne. 2) then
                              sedgez(i,j,ks,n) = merge(ZERO,sedgerz(i,j),phys_bc(3,1).eq.NO_SLIP_WALL)
                           else 
                              sedgez(i,j,ks,n) = sedgerz(i,j)
                           endif
                        elseif (phys_bc(3,1) .eq. INLET) then
                           sedgez(i,j,ks,n) = s(i,j,ks-1,n)
                        elseif (phys_bc(3,1) .eq. OUTLET) then
                           if (is_vel .and. n.eq.3) then
                              sedgez(i,j,ks,n) = MIN(sedgerz(i,j),ZERO)
                           else
                              sedgez(i,j,ks,n) = sedgerz(i,j)
                           end if
                        endif
                     else if(k .eq. ke+1) then
                        ! hi side
                        if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 2) then
                              sedgez(i,j,ke+1,n) = ZERO
                           elseif (is_vel .and. n .ne. 2) then
                              sedgez(i,j,ke+1,n) = merge(ZERO,sedgelz(i,j),phys_bc(3,2).eq.NO_SLIP_WALL)
                           else 
                              sedgez(i,j,ke+1,n) = sedgelz(i,j)
                           endif
                        elseif (phys_bc(3,2) .eq. INLET) then
                           sedgez(i,j,ke+1,n) = s(i,j,ke+1,n)
                        elseif (phys_bc(3,2) .eq. OUTLET) then
                           if (is_vel .and. n.eq.3) then
                              sedgez(i,j,ke+1,n) = MAX(sedgelz(i,j),ZERO)
                           else
                              sedgez(i,j,ke+1,n) = sedgelz(i,j)
                           end if
                        endif
                     endif
                  enddo
               enddo
               
!******************************************************************
! 7. Compute simhzx(is  :ie  ,js-1:je+1,k)
!******************************************************************

               do j=js-1,je+1
                  do i=is,ie
                     ! make slzx, srzx by updating 1D extrapolation
                     if(is_conservative(n)) then
                        slzx(i,j,kc) = slz(i,j,kc) - (dt3/hx)*(simhx(i+1,j,kp)*umac(i+1,j,k-1) - simhx(i,j,kp)*umac(i,j,k-1))
                        srzx(i,j,kc) = srz(i,j,kc) - (dt3/hx)*(simhx(i+1,j,kc)*umac(i+1,j,k  ) - simhx(i,j,kc)*umac(i,j,k  ))
                     else
                        slzx(i,j,kc) = slz(i,j,kc) - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhx(i+1,j,kp)-simhx(i,j,kp))
                        srzx(i,j,kc) = srz(i,j,kc) - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhx(i+1,j,kc)-simhx(i,j,kc))
                     endif
                     
                     ! impose lo side bc's
                     if(k .eq. ks) then
                        slzx(i,j,kc) = merge(s(is,j,k-1,n),slzx(i,j,kc),phys_bc(3,1) .eq. INLET)
                        srzx(i,j,kc) = merge(s(is,j,k-1,n),srzx(i,j,kc),phys_bc(3,1) .eq. INLET)
                        if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                           if(is_vel .and. n .eq. 3) then
                              slzx(i,j,kc) = ZERO
                              srzx(i,j,kc) = ZERO
                           else if(is_vel .and. n .ne. 3) then
                              slzx(i,j,kc) = merge(ZERO,srzx(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                              srzx(i,j,kc) = merge(ZERO,srzx(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           else
                              slzx(i,j,kc) = srzx(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! impose hi side bc's
                     if(k .eq. ke+1) then
                        slzx(i,j,kc) = merge(s(i,j,ke+1,n),slzx(i,j,kc),phys_bc(3,2) .eq. INLET)
                        srzx(i,j,kc) = merge(s(i,j,ke+1,n),srzx(i,j,kc),phys_bc(3,2) .eq. INLET)
                        if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 3) then
                              slzx(i,j,kc) = ZERO
                              srzx(i,j,kc) = ZERO
                           else if (is_vel .and. n .ne. 3) then
                              slzx(i,j,kc) = merge(ZERO,slzx(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                              srzx(i,j,kc) = merge(ZERO,slzx(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                           else
                              srzx(i,j,kc) = slzx(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! make simhzx by solving Riemann problem
                     simhzx(i,j,kc) = merge(slzx(i,j,kc),srzx(i,j,kc),wmac(i,j,k) .gt. ZERO)
                     savg = HALF*(slzx(i,j,kc)+srzx(i,j,kc))
                     simhzx(i,j,kc) = merge(simhzx(i,j,kc),savg,abs(wmac(i,j,k)) .gt. eps)
                  enddo
               enddo

!******************************************************************
! 8. Compute simhzy(is-1:ie+1,js  :je  ,k)
!******************************************************************

               do j=js,je
                  do i=is-1,ie+1
                     ! make slzy, srzy by updating 1D extrapolation
                     if(is_conservative(n)) then
                        slzy(i,j,kc) = slz(i,j,kc) - (dt3/hy)*(simhy(i,j+1,kp)*vmac(i,j+1,k-1) - simhy(i,j,kp)*vmac(i,j,k-1))
                        srzy(i,j,kc) = srz(i,j,kc) - (dt3/hy)*(simhy(i,j+1,kc)*vmac(i,j+1,k  ) - simhy(i,j,kc)*vmac(i,j,k  ))
                     else
                        slzy(i,j,kc) = slz(i,j,kc) - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhy(i,j+1,kp)-simhy(i,j,kp))
                        srzy(i,j,kc) = srz(i,j,kc) - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhy(i,j+1,kc)-simhy(i,j,kc))
                     endif
                     
                     ! impose lo side bc's
                     if(k .eq. ks) then
                        slzy(i,j,kc) = merge(s(is,j,k-1,n),slzy(i,j,kc),phys_bc(3,1) .eq. INLET)
                        srzy(i,j,kc) = merge(s(is,j,k-1,n),srzy(i,j,kc),phys_bc(3,1) .eq. INLET)
                        if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                           if(is_vel .and. n .eq. 3) then
                              slzy(i,j,kc) = ZERO
                              srzy(i,j,kc) = ZERO
                           else if(is_vel .and. n .ne. 3) then
                              slzy(i,j,kc) = merge(ZERO,srzy(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                              srzy(i,j,kc) = merge(ZERO,srzy(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                           else
                              slzy(i,j,kc) = srzy(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! impose hi side bc's
                     if(k .eq. ke+1) then
                        slzy(i,j,kc) = merge(s(i,j,ke+1,n),slzy(i,j,kc),phys_bc(3,2) .eq. INLET)
                        srzy(i,j,kc) = merge(s(i,j,ke+1,n),srzy(i,j,kc),phys_bc(3,2) .eq. INLET)
                        if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 3) then
                              slzy(i,j,kc) = ZERO
                              srzy(i,j,kc) = ZERO
                           else if (is_vel .and. n .ne. 3) then
                              slzy(i,j,kc) = merge(ZERO,slzy(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                              srzy(i,j,kc) = merge(ZERO,slzy(i,j,kc),phys_bc(3,2).eq.NO_SLIP_WALL)
                           else
                              srzy(i,j,kc) = slzy(i,j,kc)
                           endif
                        endif
                     endif
                     
                     ! make simhzy by solving Riemann problem
                     simhzy(i,j,kc) = merge(slzy(i,j,kc),srzy(i,j,kc),wmac(i,j,k) .gt. ZERO)
                     savg = HALF*(slzy(i,j,kc)+srzy(i,j,kc))
                     simhzy(i,j,kc) = merge(simhzy(i,j,kc),savg,abs(wmac(i,j,k)) .gt. eps)
                  enddo
               enddo

            endif ! end if(k .gt. ks-1)

            if(k .gt. ks) then

!******************************************************************
! 9. Compute simhxz(is  :ie+1,js-1:je+1,k-1)
!******************************************************************

               do j=js-1,je+1
                  do i=is,ie+1
                     ! make slxz, srxz by updating 1D extrapolation
                     if(is_conservative(n)) then
                        slxz(i,j,kp) = slx(i,j,kp) - (dt3/hz)*(simhz(i-1,j,kc)*wmac(i-1,j,k) - simhz(i-1,j,kp)*wmac(i-1,j,k-1))
                        srxz(i,j,kp) = srx(i,j,kp) - (dt3/hz)*(simhz(i,  j,kc)*wmac(i,  j,k) - simhz(i,  j,kp)*wmac(i,  j,k-1))
                     else
                        slxz(i,j,kp) = slx(i,j,kp) - (dt6/hz)*(wmac(i-1,j,k)+wmac(i-1,j,k-1))*(simhz(i-1,j,kc)-simhz(i-1,j,kp))
                        srxz(i,j,kp) = srx(i,j,kp) - (dt6/hz)*(wmac(i,  j,k)+wmac(i,  j,k-1))*(simhz(i,  j,kc)-simhz(i,  j,kp))
                     endif
                     
                     ! impose lo side bc's
                     if(i .eq. is) then
                        slxz(i,j,kp) = merge(s(is-1,j,k-1,n),slxz(i,j,kp),phys_bc(1,1) .eq. INLET)
                        srxz(i,j,kp) = merge(s(is-1,j,k-1,n),srxz(i,j,kp),phys_bc(1,1) .eq. INLET)
                        if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                           if(is_vel .and. n .eq. 1) then
                              slxz(i,j,kp) = ZERO
                              srxz(i,j,kp) = ZERO
                           else if(is_vel .and. n .ne. 1) then
                              slxz(i,j,kp) = merge(ZERO,srxz(i,j,kp),phys_bc(1,1) .eq. NO_SLIP_WALL)
                              srxz(i,j,kp) = merge(ZERO,srxz(i,j,kp),phys_bc(1,1) .eq. NO_SLIP_WALL)
                           else
                              slxz(i,j,kp) = srxz(i,j,kp)
                           endif
                        endif
                     endif
                     
                     ! impose hi side bc's
                     if(i .eq. ie+1) then
                        slxz(i,j,kp) = merge(s(ie+1,j,k-1,n),slxz(i,j,kp),phys_bc(1,2) .eq. INLET)
                        srxz(i,j,kp) = merge(s(ie+1,j,k-1,n),srxz(i,j,kp),phys_bc(1,2) .eq. INLET)
                        if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 1) then
                              slxz(i,j,kp) = ZERO
                              srxz(i,j,kp) = ZERO
                           else if (is_vel .and. n .ne. 1) then
                              slxz(i,j,kp) = merge(ZERO,slxz(i,j,kp),phys_bc(1,2).eq.NO_SLIP_WALL)
                              srxz(i,j,kp) = merge(ZERO,slxz(i,j,kp),phys_bc(1,2).eq.NO_SLIP_WALL)
                           else
                              srxz(i,j,kp) = slxz(i,j,kp)
                           endif
                        endif
                     endif
                     
                     ! make simhxz by solving Riemann problem
                     simhxz(i,j,kp) = merge(slxz(i,j,kp),srxz(i,j,kp),umac(i,j,k-1) .gt. ZERO)
                     savg = HALF*(slxz(i,j,kp)+srxz(i,j,kp))
                     simhxz(i,j,kp) = merge(simhxz(i,j,kp),savg,abs(umac(i,j,k-1)) .gt. eps)
                  enddo
               enddo

!******************************************************************
! 10.Compute simhyz(is-1:ie+1,js  :je+1,k-1)
!******************************************************************

               do j=js,je+1
                  do i=is-1,ie+1
                     ! make slyz, sryz by updating 1D extrapolation
                     if(is_conservative(n)) then
                        slyz(i,j,kp) = sly(i,j,kp) - (dt3/hz)*(simhz(i,j-1,kc)*wmac(i,j-1,k) - simhz(i,j-1,kp)*wmac(i,j-1,k-1))
                        sryz(i,j,kp) = sry(i,j,kp) - (dt3/hz)*(simhz(i,j,  kc)*wmac(i,j,  k) - simhz(i,j,  kp)*wmac(i,j,  k-1))
                     else
                        slyz(i,j,kp) = sly(i,j,kp) - (dt6/hz)*(wmac(i,j-1,k)+wmac(i,j-1,k-1))*(simhz(i,j-1,kc)-simhz(i,j-1,kp))
                        sryz(i,j,kp) = sry(i,j,kp) - (dt6/hz)*(wmac(i,j,  k)+wmac(i,j,  k-1))*(simhz(i,j,  kc)-simhz(i,j,  kp))
                     endif
                     
                     ! impose lo side bc's
                     if(j .eq. js) then
                        slyz(i,j,kp) = merge(s(is,j-1,k-1,n),slyz(i,j,kp),phys_bc(2,1) .eq. INLET)
                        sryz(i,j,kp) = merge(s(is,j-1,k-1,n),sryz(i,j,kp),phys_bc(2,1) .eq. INLET)
                        if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                           if(is_vel .and. n .eq. 2) then
                              slyz(i,j,kp) = ZERO
                              sryz(i,j,kp) = ZERO
                           else if(is_vel .and. n .ne. 2) then
                              slyz(i,j,kp) = merge(ZERO,sryz(i,j,kp),phys_bc(2,1) .eq. NO_SLIP_WALL)
                              sryz(i,j,kp) = merge(ZERO,sryz(i,j,kp),phys_bc(2,1) .eq. NO_SLIP_WALL)
                           else
                              slyz(i,j,kp) = sryz(i,j,kp)
                           endif
                        endif
                     endif
                     
                     ! impose hi side bc's
                     if(j .eq. je+1) then
                        slyz(i,j,kp) = merge(s(i,je+1,k-1,n),slyz(i,j,kp),phys_bc(2,2) .eq. INLET)
                        sryz(i,j,kp) = merge(s(i,je+1,k-1,n),sryz(i,j,kp),phys_bc(2,2) .eq. INLET)
                        if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                           if (is_vel .and. n .eq. 2) then
                              slyz(i,j,kp) = ZERO
                              sryz(i,j,kp) = ZERO
                           else if (is_vel .and. n .ne. 2) then
                              slyz(i,j,kp) = merge(ZERO,slyz(i,j,kp),phys_bc(2,2).eq.NO_SLIP_WALL)
                              sryz(i,j,kp) = merge(ZERO,slyz(i,j,kp),phys_bc(2,2).eq.NO_SLIP_WALL)
                           else
                              sryz(i,j,kp) = slyz(i,j,kp)
                           endif
                        endif
                     endif
                     
                     ! make simhyz by solving Riemann problem
                     simhyz(i,j,kp) = merge(slyz(i,j,kp),sryz(i,j,kp),vmac(i,j,k-1) .gt. ZERO)
                     savg = HALF*(slyz(i,j,kp)+sryz(i,j,kp))
                     simhyz(i,j,kp) = merge(simhyz(i,j,kp),savg,abs(vmac(i,j,k-1)) .gt. eps)
                  enddo
               enddo

!******************************************************************
! 11.Compute sedgex(is  :ie+1,js  :je,  k-1)
!******************************************************************
               
               do j=js,je
                  do i=is,ie+1
                     ! make sedgelx, sedgerx
                     if(is_conservative(n)) then
                        sedgelx(i,j) = slx(i,j,kp) &
                             - (dt2/hy)*(simhyz(i-1,j+1,kp)*vmac(i-1,j+1,k-1) - simhyz(i-1,j,kp)*vmac(i-1,j,k-1)) &
                             - (dt2/hz)*(simhzy(i-1,j  ,kc)*wmac(i-1,j  ,k  ) - simhzy(i-1,j,kp)*wmac(i-1,j,k-1))
                        sedgerx(i,j) = srx(i,j,kp) &
                             - (dt2/hy)*(simhyz(i,  j+1,kp)*vmac(i,  j+1,k-1) - simhyz(i,  j,kp)*vmac(i,  j,k-1)) &
                             - (dt2/hz)*(simhzy(i,  j  ,kc)*wmac(i,  j  ,k  ) - simhzy(i,  j,kp)*wmac(i,  j,k-1))
                     else
                        sedgelx(i,j) = slx(i,j,kp) &
                             - (dt4/hy)*(vmac(i-1,j+1,k-1)+vmac(i-1,j,k-1))*(simhyz(i-1,j+1,kp)-simhyz(i-1,j,kp)) &
                             - (dt4/hz)*(wmac(i-1,j  ,k  )+wmac(i-1,j,k-1))*(simhzy(i-1,j  ,kc)-simhzy(i-1,j,kp))
                        sedgerx(i,j) = srx(i,j,kp) &
                             - (dt4/hy)*(vmac(i,  j+1,k-1)+vmac(i,  j,k-1))*(simhyz(i,  j+1,kp)-simhyz(i,  j,kp)) &
                             - (dt4/hz)*(wmac(i,  j  ,k  )+wmac(i,  j,k-1))*(simhzy(i,  j  ,kc)-simhzy(i,  j,kp))
                     endif
                     
                     ! if use_minion is true, we have already accounted for source terms
                     ! in slx and srx; otherwise, we need to account for them here.
                     if(.not. use_minion) then
                        sedgelx(i,j) = sedgelx(i,j) + dt2*force(i-1,j,k-1,n)
                        sedgerx(i,j) = sedgerx(i,j) + dt2*force(i  ,j,k-1,n)
                     endif
                     
                     ! make sedgex by solving Riemann problem
                     ! boundary conditions enforced outside of i,j,k loop
                     sedgex(i,j,k-1,n) = merge(sedgelx(i,j),sedgerx(i,j),umac(i,j,k-1) .gt. ZERO)
                     savg = HALF*(sedgelx(i,j)+sedgerx(i,j))
                     sedgex(i,j,k-1,n) = merge(sedgex(i,j,k-1,n),savg,abs(umac(i,j,k-1)) .gt. eps)
                  enddo
               enddo

               ! sedgex boundary conditions
               if(i .eq. is) then
                  ! lo side
                  if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                     if (is_vel .and. n .eq. 1) then
                        sedgex(is,j,k-1,n) = ZERO
                     elseif (is_vel .and. n .ne. 1) then
                        sedgex(is,j,k-1,n) = merge(ZERO,sedgerx(is,j),phys_bc(1,1).eq.NO_SLIP_WALL)
                     else
                        sedgex(is,j,k-1,n) = sedgerx(is,j)
                     endif
                  elseif (phys_bc(1,1) .eq. INLET) then
                     sedgex(is,j,k-1,n) = s(is-1,j,k-1,n)
                  elseif (phys_bc(1,1) .eq. OUTLET) then
                     if (is_vel .and. n.eq.1) then
                        sedgex(is,j,k-1,n) = MIN(sedgerx(is,j),ZERO)
                     else
                        sedgex(is,j,k-1,n) = sedgerx(is,j)
                     end if
                  endif
               else if(i .eq. ie+1) then
                  ! hi side
                  if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                     if (is_vel .and. n .eq. 1) then
                        sedgex(ie+1,j,k-1,n) = ZERO
                     else if (is_vel .and. n .ne. 1) then
                        sedgex(ie+1,j,k-1,n) = merge(ZERO,sedgelx(ie+1,j),phys_bc(1,2).eq.NO_SLIP_WALL)
                     else 
                        sedgex(ie+1,j,k-1,n) = sedgelx(ie+1,j)
                     endif
                  elseif (phys_bc(1,2) .eq. INLET) then
                     sedgex(ie+1,j,k-1,n) = s(ie+1,j,k-1,n)
                  elseif (phys_bc(1,2) .eq. OUTLET) then
                     if (is_vel .and. n.eq.1) then
                        sedgex(ie+1,j,k-1,n) = MAX(sedgelx(ie+1,j),ZERO)
                     else
                        sedgex(ie+1,j,k-1,n) = sedgelx(ie+1,j)
                     end if
                  endif
               endif


!******************************************************************
! 12.Compute sedgey(is  :ie  ,js  :je+1,k-1)
!******************************************************************
               
               do j=js,je+1
                  do i=is,ie
                     ! make sedgely, sedgery
                     if(is_conservative(n)) then
                        sedgely(i,j) = sly(i,j,kp) &
                             - (dt2/hx)*(simhxz(i+1,j-1,kp)*umac(i+1,j-1,k-1) - simhxz(i,j-1,kp)*umac(i,j-1,k-1)) &
                             - (dt2/hz)*(simhzx(i  ,j-1,kc)*wmac(i  ,j-1,k  ) - simhzx(i,j-1,kp)*wmac(i,j-1,k-1))
                        sedgery(i,j) = sry(i,j,kp) &
                             - (dt2/hx)*(simhxz(i+1,j,  kp)*umac(i+1,j,  k-1) - simhxz(i,j,  kp)*umac(i,j,  k-1)) &
                             - (dt2/hz)*(simhzx(i  ,j,  kc)*wmac(i  ,j,  k  ) - simhzx(i,j,  kp)*wmac(i,j,  k-1))
                     else
                        sedgely(i,j) = sly(i,j,kp) &
                             - (dt4/hx)*(umac(i+1,j-1,k-1)+umac(i,j-1,k-1))*(simhxz(i+1,j-1,kp)-simhxz(i,j-1,kp)) &
                             - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k  ))*(simhzx(i  ,j-1,kc)-simhzx(i,j-1,kp))
                        sedgery(i,j) = sry(i,j,kp) &
                             - (dt4/hx)*(umac(i+1,j,  k-1)+umac(i,j,  k-1))*(simhxz(i+1,j,  kp)-simhxz(i,j,  kp)) &
                             - (dt4/hz)*(wmac(i  ,j,  k  )+wmac(i,j,  k-1))*(simhzx(i  ,j,  kc)-simhzx(i,j,  kp))
                     endif
                     
                     ! if use_minion is true, we have already accounted for source terms
                     ! in sly and sry; otherwise, we need to account for them here.
                     if(.not. use_minion) then
                        sedgely(i,j) = sedgely(i,j) + dt2*force(i,j-1,k-1,n)
                        sedgery(i,j) = sedgery(i,j) + dt2*force(i,j  ,k-1,n)
                     endif
                     
                     ! make sedgey by solving Riemann problem
                     ! boundary conditions enforced outside of i,j,k loop
                     sedgey(i,j,k-1,n) = merge(sedgely(i,j),sedgery(i,j),vmac(i,j,k-1) .gt. ZERO)
                     savg = HALF*(sedgely(i,j)+sedgery(i,j))
                     sedgey(i,j,k-1,n) = merge(sedgey(i,j,k-1,n),savg,abs(vmac(i,j,k-1)) .gt. eps)
                  enddo
               enddo
               
               ! sedgey boundary conditions
               if(j .eq. js) then
                  ! lo side
                  if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                     if (is_vel .and. n .eq. 2) then
                        sedgey(i,js,k-1,n) = ZERO
                     elseif (is_vel .and. n .ne. 2) then
                        sedgey(i,js,k-1,n) = merge(ZERO,sedgery(i,js),phys_bc(2,1).eq.NO_SLIP_WALL)
                     else 
                        sedgey(i,js,k-1,n) = sedgery(i,js)
                     endif
                  elseif (phys_bc(2,1) .eq. INLET) then
                     sedgey(i,js,k-1,n) = s(i,js-1,k-1,n)
                  elseif (phys_bc(2,1) .eq. OUTLET) then
                     if (is_vel .and. n.eq.2) then
                        sedgey(i,js,k-1,n) = MIN(sedgery(i,js),ZERO)
                     else
                        sedgey(i,js,k-1,n) = sedgery(i,js)
                     end if
                  endif
               else if(j .eq. je+1) then
                  ! hi side
                  if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                     if (is_vel .and. n .eq. 2) then
                        sedgey(i,je+1,k-1,n) = ZERO
                     elseif (is_vel .and. n .ne. 2) then
                        sedgey(i,je+1,k-1,n) = merge(ZERO,sedgely(i,je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
                     else 
                        sedgey(i,je+1,k-1,n) = sedgely(i,je+1)
                     endif
                  elseif (phys_bc(2,2) .eq. INLET) then
                     sedgey(i,je+1,k-1,n) = s(i,je+1,k-1,n)
                  elseif (phys_bc(2,2) .eq. OUTLET) then
                     if (is_vel .and. n.eq.2) then
                        sedgey(i,je+1,k-1,n) = MAX(sedgely(i,je+1),ZERO)
                     else
                        sedgey(i,je+1,k-1,n) = sedgely(i,je+1)
                     end if
                  endif
               endif

            endif ! end if(k .gt. ks)

!******************************************************************
! 13. Cycle indeces
!******************************************************************

            kc = 3 - kc
            kp = 3 - kp

         enddo ! end loop over k
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

      end subroutine mkflux_lowmemory_3d

end module mkflux_lowmemory_module
