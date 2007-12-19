module velpred_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: velpred

contains

  subroutine velpred(nlevs,u,umac,force,dx,dt,the_bc_level,use_minion, &
                     use_godunov_debug,mla)

    use ml_restriction_module, only: ml_edge_restriction

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical        , intent(in)    :: use_minion
    logical        , intent(in)    :: use_godunov_debug
    type(ml_layout), intent(inout) :: mla

    ! local
    integer :: i,n,dm,ng,comp
    integer :: lo(u(1)%dim)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: ump(:,:,:,:)
    real(kind=dp_t), pointer:: vmp(:,:,:,:)
    real(kind=dp_t), pointer:: wmp(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)

    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       ! Create the edge states to be used for the MAC velocity 
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          up  => dataptr(u(n),i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          fp  => dataptr(force(n),i)
          lo =  lwb(get_box(u(n),i))
          select case (dm)
          case (2)
             if(use_godunov_debug) then
                call velpred_debug_2d(up(:,:,1,:), &
                                      ump(:,:,1,1),  vmp(:,:,1,1), &
                                      fp(:,:,1,:), &
                                      lo, dx(n,:), dt, &
                                      the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                      ng, use_minion)
             else
                call velpred_2d(up(:,:,1,:), &
                                ump(:,:,1,1),  vmp(:,:,1,1), &
                                fp(:,:,1,:), &
                                lo, dx(n,:), dt, &
                                the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                ng, use_minion)
             endif
          case (3)
             wmp  => dataptr(umac(n,3), i)
             if(use_godunov_debug) then
                call velpred_debug_3d(up(:,:,:,:), &
                                      ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                                      fp(:,:,:,:), &
                                      lo, dx(n,:), dt, &
                                      the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                      ng, use_minion)
                else
                call velpred_3d(up(:,:,:,:), &
                                ump(:,:,:,1),  vmp(:,:,:,1), wmp(:,:,:,1), &
                                fp(:,:,:,:), &
                                lo, dx(n,:), dt, &
                                the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                ng, use_minion)
             endif
          end select
       end do

       do comp = 1, dm
          call multifab_fill_boundary(umac(n,comp))
       enddo
       
    enddo ! end loop over levels

    do n = nlevs,2,-1
       do comp = 1,dm
          call ml_edge_restriction(umac(n-1,comp),umac(n,comp),mla%mba%rr(n-1,:),comp)
       end do
    end do
    
  end subroutine velpred


  subroutine velpred_2d(u,umac,vmac,force,lo,dx,dt,phys_bc,adv_bc,ng,use_minion)

    use bc_module
    use slope_module
    use bl_constants_module

    integer         ,intent(in) :: lo(2)
    integer         ,intent(in) :: ng

    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)

    real(kind=dp_t) ,intent(in) :: dx(:),dt
    integer         ,intent(in) :: phys_bc(:,:)
    integer         ,intent(in) :: adv_bc(:,:,:)
    logical         ,intent(in) :: use_minion

    ! Local variables
    real(kind=dp_t), allocatable::  slopex(:,:,:)
    real(kind=dp_t), allocatable::  slopey(:,:,:)

    real(kind=dp_t) hx, hy, dt2, dt4, uavg

    integer :: hi(2)
    integer :: slope_order = 4
    logical :: test

    real(kind=dp_t) :: abs_eps, eps, umax
    integer :: i,j,is,js,ie,je
    integer :: jc,jp ! "current" and "previous" j

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable:: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    real(kind=dp_t), allocatable:: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable:: umacl(:),umacr(:)
    real(kind=dp_t), allocatable:: vmacl(:),vmacr(:)

    hi(1) = lo(1) + size(u,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(u,dim=2) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    call slopex_2d(u,slopex,lo,ng,2,adv_bc,slope_order)
    call slopey_2d(u,slopey,lo,ng,2,adv_bc,slope_order)

    ! Note: All of these arrays are allocated to exactly the 
    ! size they need to be in order to compute MAC velocities on
    ! a domain with faces from lo(1):hi(1)+1,lo(2):hi(2)+1

    !***********************
    ! Normal predictor terms
    !***********************

    ! lo(1):hi(1)+1 in the x-direction
    ! 2 rows needed in y-direction
    allocate(ulx  (lo(1):hi(1)+1,2,2))
    allocate(urx  (lo(1):hi(1)+1,2,2))
    allocate(uimhx(lo(1):hi(1)+1,2,2))

    ! lo(1)-1:hi(1)+1 in the x-direction
    ! 2 rows needed in y-direction
    allocate(uly  (lo(1)-1:hi(1)+1,2,2))
    allocate(ury  (lo(1)-1:hi(1)+1,2,2))
    allocate(uimhy(lo(1)-1:hi(1)+1,2,2))

    !***************
    ! MAC velocities
    !***************

    ! lo(1):hi(1)+1 in x-direction
    allocate(umacl(lo(1):hi(1)+1))
    allocate(umacr(lo(1):hi(1)+1))

    ! lo(1):hi(1) in x-direction
    allocate(vmacl(lo(1):hi(1)))
    allocate(vmacr(lo(1):hi(1)))

    abs_eps = 1.0d-8

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    ! Compute eps, which is relative to the max velocity
    umax = abs(u(is,js,1))
    do j = js,je
       do i = is,ie
          umax = max(umax,abs(u(i,j,1)))
          umax = max(umax,abs(u(i,j,2)))
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
    !     1. Compute uimhx(is:ie+1,j)
    !     if(j .gt. js-1) then
    !        2. Compute uimhy(is-1:ie+1,j)
    !        3. Compute vmac(is:ie,j)
    !     endif
    !     if(j .gt. js) then
    !        4. Compute umac(is:ie+1,j-1)
    !     endif
    !     5. Cycle indices
    !  enddo
    !
    !*************************************
    ! End pseudo code
    !*************************************

    jc = 1
    jp = 2

    do j=js-1,je+1

       !******************************************************************
       ! 1. Compute uimhx(is:ie+1,j)
       !******************************************************************

       do i=is,ie+1
          ! extrapolate both components of velocity to left face
          ulx(i,jc,1) = u(i-1,j,1) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,1)
          ulx(i,jc,2) = u(i-1,j,2) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,2)

          ! extrapolate both components of velocity to right face
          urx(i,jc,1) = u(i  ,j,1) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,1)
          urx(i,jc,2) = u(i  ,j,2) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,2)

          ! add source terms
          if(use_minion) then
             ulx(i,jc,1) = ulx(i,jc,1) + dt2*force(i-1,j,1)
             ulx(i,jc,2) = ulx(i,jc,2) + dt2*force(i-1,j,2)
             urx(i,jc,1) = urx(i,jc,1) + dt2*force(i  ,j,1)
             urx(i,jc,2) = urx(i,jc,2) + dt2*force(i  ,j,2)
          endif

          ! impose lo side bc's
          if(i .eq. is) then
             ulx(i,jc,1) = merge(u(is-1,j,1),ulx(i,jc,1),phys_bc(1,1) .eq. INLET)
             urx(i,jc,1) = merge(u(is-1,j,1),urx(i,jc,1),phys_bc(1,1) .eq. INLET)
             ulx(i,jc,2) = merge(u(is-1,j,2),ulx(i,jc,2),phys_bc(1,1) .eq. INLET)
             urx(i,jc,2) = merge(u(is-1,j,2),urx(i,jc,2),phys_bc(1,1) .eq. INLET)
             if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                ulx(i,jc,1) = ZERO
                urx(i,jc,1) = ZERO
                ulx(i,jc,2) = merge(ZERO,urx(i,jc,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                urx(i,jc,2) = merge(ZERO,urx(i,jc,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
             endif
          endif

          ! impose hi side bc's
          if(i .eq. ie+1) then
             ulx(i,jc,1) = merge(u(ie+1,j,1),ulx(i,jc,1),phys_bc(1,2) .eq. INLET)
             urx(i,jc,1) = merge(u(ie+1,j,1),urx(i,jc,1),phys_bc(1,2) .eq. INLET)
             ulx(i,jc,2) = merge(u(ie+1,j,2),ulx(i,jc,2),phys_bc(1,2) .eq. INLET)
             urx(i,jc,2) = merge(u(ie+1,j,2),urx(i,jc,2),phys_bc(1,2) .eq. INLET)
             if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                ulx(i,jc,1) = ZERO
                urx(i,jc,1) = ZERO
                ulx(i,jc,2) = merge(ZERO,ulx(i,jc,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                urx(i,jc,2) = merge(ZERO,ulx(i,jc,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
             endif
          endif

          ! make normal component of uimhx by first solving a normal Riemann problem
          uavg = HALF*(ulx(i,jc,1)+urx(i,jc,1))
          test = ((ulx(i,jc,1) .le. ZERO .and. urx(i,jc,1) .ge. ZERO) .or. &
               (abs(ulx(i,jc,1)+urx(i,jc,1)) .lt. eps))
          uimhx(i,jc,1) = merge(ulx(i,jc,1),urx(i,jc,1),uavg .gt. ZERO)
          uimhx(i,jc,1) = merge(ZERO,uimhx(i,jc,1),test)

          ! now upwind to get transverse component of uimhx
          uimhx(i,jc,2) = merge(ulx(i,jc,2),urx(i,jc,2),uimhx(i,jc,1).gt.ZERO)
          uavg = HALF*(ulx(i,jc,2)+urx(i,jc,2))
          uimhx(i,jc,2) = merge(uavg,uimhx(i,jc,2),abs(uimhx(i,jc,1)).lt.eps)
       enddo

       if(j .gt. js-1) then

          !******************************************************************
          ! 2. Compute uimhy(is-1:ie+1,j)
          !******************************************************************

          do i=is-1,ie+1
             ! extrapolate both components of velocity to left face
             uly(i,jc,1) = u(i,j-1,1) + (HALF - dt2*max(ZERO,u(i,j-1,2)/hy))*slopey(i,j-1,1)
             uly(i,jc,2) = u(i,j-1,2) + (HALF - dt2*max(ZERO,u(i,j-1,2)/hy))*slopey(i,j-1,2)

             ! extrapolate both components of velocity to right face
             ury(i,jc,1) = u(i,j  ,1) - (HALF + dt2*min(ZERO,u(i,j  ,2)/hy))*slopey(i,j  ,1)
             ury(i,jc,2) = u(i,j  ,2) - (HALF + dt2*min(ZERO,u(i,j  ,2)/hy))*slopey(i,j  ,2)

             ! add source terms
             if(use_minion) then
                uly(i,jc,1) = uly(i,jc,1) + dt2*force(i,j-1,1)
                uly(i,jc,2) = uly(i,jc,2) + dt2*force(i,j-1,2)
                ury(i,jc,1) = ury(i,jc,1) + dt2*force(i,j  ,1)
                ury(i,jc,2) = ury(i,jc,2) + dt2*force(i,j  ,2)
             endif

             ! impose lo side bc's
             if(j .eq. js) then
                uly(i,jc,1) = merge(u(i,js-1,1),uly(i,jc,1),phys_bc(2,1) .eq. INLET)
                ury(i,jc,1) = merge(u(i,js-1,1),ury(i,jc,1),phys_bc(2,1) .eq. INLET)
                uly(i,jc,2) = merge(u(i,js-1,2),uly(i,jc,2),phys_bc(2,1) .eq. INLET)
                ury(i,jc,2) = merge(u(i,js-1,2),ury(i,jc,2),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   uly(i,jc,1) = merge(ZERO,ury(i,jc,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,jc,1) = merge(ZERO,ury(i,jc,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uly(i,jc,2) = ZERO
                   ury(i,jc,2) = ZERO
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                uly(i,jc,1) = merge(u(i,je+1,1),uly(i,jc,1),phys_bc(2,2) .eq. INLET)
                ury(i,jc,1) = merge(u(i,je+1,1),ury(i,jc,1),phys_bc(2,2) .eq. INLET)
                uly(i,jc,2) = merge(u(i,je+1,2),uly(i,jc,2),phys_bc(2,2) .eq. INLET)
                ury(i,jc,2) = merge(u(i,je+1,2),ury(i,jc,2),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   uly(i,jc,1) = merge(ZERO,uly(i,jc,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,jc,1) = merge(ZERO,uly(i,jc,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uly(i,jc,2) = ZERO
                   ury(i,jc,2) = ZERO
                endif
             endif

             ! make normal component of uimhx by first solving a normal Riemann problem
             uavg = HALF*(uly(i,jc,2)+ury(i,jc,2))
             test = ((uly(i,jc,2) .le. ZERO .and. ury(i,jc,2) .ge. ZERO) .or. &
                  (abs(uly(i,jc,2)+ury(i,jc,2)) .lt. eps))
             uimhy(i,jc,2) = merge(uly(i,jc,2),ury(i,jc,2),uavg .gt. ZERO)
             uimhy(i,jc,2) = merge(ZERO,uimhy(i,jc,2),test)

             ! now upwind to get transverse component of uimhy
             uimhy(i,jc,1) = merge(uly(i,jc,1),ury(i,jc,1),uimhy(i,jc,2).gt.ZERO)
             uavg = HALF*(uly(i,jc,1)+ury(i,jc,1))
             uimhy(i,jc,1) = merge(uavg,uimhy(i,jc,1),abs(uimhy(i,jc,2)).lt.eps)
          enddo

          !******************************************************************
          ! 3. Compute vmac(is:ie,j)
          !******************************************************************

          do i=is,ie
             ! extrapolate to edges
             vmacl(i) = uly(i,jc,2) &
                  - (dt4/hx)*(uimhx(i+1,jp,1)+uimhx(i,jp,1))*(uimhx(i+1,jp,2)-uimhx(i,jp,2))
             vmacr(i) = ury(i,jc,2) &
                  - (dt4/hx)*(uimhx(i+1,jc,1)+uimhx(i,jc,1))*(uimhx(i+1,jc,2)-uimhx(i,jc,2))

             ! if use_minion is true, we have already accounted for source terms
             ! in uly and ury; otherwise, we need to account for them here.
             if(.not. use_minion) then
                vmacl(i) = vmacl(i) + dt2*force(i,j-1,2)
                vmacr(i) = vmacr(i) + dt2*force(i,j  ,2)
             endif

             ! solve Riemann problem
             uavg = HALF*(vmacl(i)+vmacr(i))
             test = ((vmacl(i) .le. ZERO .and. vmacr(i) .ge. ZERO) .or. &
                  (abs(vmacl(i)+vmacr(i)) .lt. eps))
             vmac(i,j) = merge(vmacl(i),vmacr(i),uavg .gt. ZERO)
             vmac(i,j) = merge(ZERO,vmac(i,j),test)

             ! Apply boundary conditions
             if(j .eq. js) then
                ! lo side
                if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   vmac(i,js) = ZERO
                elseif (phys_bc(2,1) .eq. INLET) then
                   vmac(i,js) = u(i,js-1,2)
                elseif (phys_bc(2,1) .eq. OUTLET) then
                   vmac(i,js) = min(vmacr(i),ZERO)
                endif
             else if(j .eq. je+1) then
                ! hi side
                if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   vmac(i,je+1) = ZERO
                elseif (phys_bc(2,2) .eq. INLET) then
                   vmac(i,je+1) = u(i,je+1,2)
                elseif (phys_bc(2,2) .eq. OUTLET) then
                   vmac(i,je+1) = max(vmacl(i),ZERO)
                endif
             endif
          enddo

       endif ! end if(j .gt. js-1)

       if(j .gt. js) then

          !******************************************************************
          ! 4. Compute umac(is:ie+1,j-1)
          !******************************************************************

          do i=is,ie+1
             ! extrapolate to edges
             umacl(i) = ulx(i,jp,1) &
                  - (dt4/hy)*(uimhy(i-1,jc,2)+uimhy(i-1,jp,2))*(uimhy(i-1,jc,1)-uimhy(i-1,jp,1))
             umacr(i) = urx(i,jp,1) &
                  - (dt4/hy)*(uimhy(i  ,jc,2)+uimhy(i  ,jp,2))*(uimhy(i  ,jc,1)-uimhy(i  ,jp,1))

             ! if use_minion is true, we have already accounted for source terms
             ! in ulx and urx; otherwise, we need to account for them here.
             if(.not. use_minion) then
                umacl(i) = umacl(i) + dt2*force(i-1,j-1,1)
                umacr(i) = umacr(i) + dt2*force(i  ,j-1,1)
             endif

             ! solve Riemann problem
             uavg = HALF*(umacl(i)+umacr(i))
             test = ((umacl(i) .le. ZERO .and. umacr(i) .ge. ZERO) .or. &
                  (abs(umacl(i)+umacr(i)) .lt. eps))
             umac(i,j-1) = merge(umacl(i),umacr(i),uavg .gt. ZERO)
             umac(i,j-1) = merge(ZERO,umac(i,j-1),test)

             ! Apply boundary conditions
             if(i .eq. is) then
                ! lo side
                if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   umac(is,j-1) = ZERO
                elseif (phys_bc(1,1) .eq. INLET) then
                   umac(is,j-1) = u(is-1,j-1,1)
                elseif (phys_bc(1,1) .eq. OUTLET) then
                   umac(is,j-1) = min(umacr(is),ZERO)
                endif
             else if(i .eq. ie+1) then
                ! hi side
                if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   umac(ie+1,j-1) = ZERO
                elseif (phys_bc(1,2) .eq. INLET) then
                   umac(ie+1,j-1) = u(ie+1,j-1,1)
                elseif (phys_bc(1,2) .eq. OUTLET) then
                   umac(ie+1,j-1) = max(umacl(ie+1),ZERO)
                endif
             endif
          enddo

       endif ! end if(j .gt. js)

       !******************************************************************
       ! 5. Cycle indices
       !******************************************************************

       jc = 3 - jc
       jp = 3 - jp

    enddo ! end loop over j

    deallocate(slopex)
    deallocate(slopey)

    deallocate(ulx)
    deallocate(urx)
    deallocate(uly)
    deallocate(ury)
    deallocate(uimhx)
    deallocate(uimhy)

    deallocate(umacl)
    deallocate(umacr)
    deallocate(vmacl)
    deallocate(vmacr)

  end subroutine velpred_2d


  subroutine velpred_debug_2d(u,umac,vmac,force,lo,dx,dt,phys_bc,adv_bc,ng,use_minion)

    use bc_module
    use slope_module
    use bl_constants_module

    integer         ,intent(in) :: lo(2)
    integer         ,intent(in) :: ng

    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)

    real(kind=dp_t) ,intent(in) :: dx(:),dt
    integer         ,intent(in) :: phys_bc(:,:)
    integer         ,intent(in) :: adv_bc(:,:,:)
    logical         ,intent(in) :: use_minion

    ! Local variables
    real(kind=dp_t), allocatable::  slopex(:,:,:)
    real(kind=dp_t), allocatable::  slopey(:,:,:)

    real(kind=dp_t) hx, hy, dt2, dt4, uavg

    integer :: hi(2)
    integer :: slope_order = 4
    logical :: test

    real(kind=dp_t) :: abs_eps, eps, umax
    integer :: i,j,is,js,ie,je

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable:: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    real(kind=dp_t), allocatable:: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable:: umacl(:,:),umacr(:,:)
    real(kind=dp_t), allocatable:: vmacl(:,:),vmacr(:,:)

    hi(1) = lo(1) + size(u,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(u,dim=2) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse direction
    allocate(ulx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(urx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

    allocate(uly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(ury  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2)))

    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1))

    call slopex_2d(u,slopex,lo,ng,2,adv_bc,slope_order)
    call slopey_2d(u,slopey,lo,ng,2,adv_bc,slope_order)

    abs_eps = 1.0d-8

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    ! Compute eps, which is relative to the max velocity
    umax = abs(u(is,js,1))
    do j = js,je
       do i = is,ie
          umax = max(umax,abs(u(i,j,1)))
          umax = max(umax,abs(u(i,j,2)))
       end do
    end do
    if(umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    do j=js-1,je+1
       do i=is,ie+1
          ! extrapolate both components of velocity to left face
          ulx(i,j,1) = u(i-1,j,1) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,1)
          ulx(i,j,2) = u(i-1,j,2) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,2)

          ! extrapolate both components of velocity to right face
          urx(i,j,1) = u(i  ,j,1) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,1)
          urx(i,j,2) = u(i  ,j,2) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,2)

          ! add source terms
          if(use_minion) then
             ulx(i,j,1) = ulx(i,j,1) + dt2*force(i-1,j,1)
             ulx(i,j,2) = ulx(i,j,2) + dt2*force(i-1,j,2)
             urx(i,j,1) = urx(i,j,1) + dt2*force(i  ,j,1)
             urx(i,j,2) = urx(i,j,2) + dt2*force(i  ,j,2)
          endif

          ! impose lo side bc's
          if(i .eq. is) then
             ulx(i,j,1) = merge(u(is-1,j,1),ulx(i,j,1),phys_bc(1,1) .eq. INLET)
             urx(i,j,1) = merge(u(is-1,j,1),urx(i,j,1),phys_bc(1,1) .eq. INLET)
             ulx(i,j,2) = merge(u(is-1,j,2),ulx(i,j,2),phys_bc(1,1) .eq. INLET)
             urx(i,j,2) = merge(u(is-1,j,2),urx(i,j,2),phys_bc(1,1) .eq. INLET)
             if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                ulx(i,j,1) = ZERO
                urx(i,j,1) = ZERO
                ulx(i,j,2) = merge(ZERO,urx(i,j,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                urx(i,j,2) = merge(ZERO,urx(i,j,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
             endif
          endif

          ! impose hi side bc's
          if(i .eq. ie+1) then
             ulx(i,j,1) = merge(u(ie+1,j,1),ulx(i,j,1),phys_bc(1,2) .eq. INLET)
             urx(i,j,1) = merge(u(ie+1,j,1),urx(i,j,1),phys_bc(1,2) .eq. INLET)
             ulx(i,j,2) = merge(u(ie+1,j,2),ulx(i,j,2),phys_bc(1,2) .eq. INLET)
             urx(i,j,2) = merge(u(ie+1,j,2),urx(i,j,2),phys_bc(1,2) .eq. INLET)
             if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                ulx(i,j,1) = ZERO
                urx(i,j,1) = ZERO
                ulx(i,j,2) = merge(ZERO,ulx(i,j,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                urx(i,j,2) = merge(ZERO,ulx(i,j,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
             endif
          endif

          ! make normal component of uimhx by first solving a normal Riemann problem
          uavg = HALF*(ulx(i,j,1)+urx(i,j,1))
          test = ((ulx(i,j,1) .le. ZERO .and. urx(i,j,1) .ge. ZERO) .or. &
               (abs(ulx(i,j,1)+urx(i,j,1)) .lt. eps))
          uimhx(i,j,1) = merge(ulx(i,j,1),urx(i,j,1),uavg .gt. ZERO)
          uimhx(i,j,1) = merge(ZERO,uimhx(i,j,1),test)

          ! now upwind to get transverse component of uimhx
          uimhx(i,j,2) = merge(ulx(i,j,2),urx(i,j,2),uimhx(i,j,1).gt.ZERO)
          uavg = HALF*(ulx(i,j,2)+urx(i,j,2))
          uimhx(i,j,2) = merge(uavg,uimhx(i,j,2),abs(uimhx(i,j,1)).lt.eps)
       enddo
    enddo

    do j=js,je+1
       do i=is-1,ie+1
          ! extrapolate both components of velocity to left face
          uly(i,j,1) = u(i,j-1,1) + (HALF - dt2*max(ZERO,u(i,j-1,2)/hy))*slopey(i,j-1,1)
          uly(i,j,2) = u(i,j-1,2) + (HALF - dt2*max(ZERO,u(i,j-1,2)/hy))*slopey(i,j-1,2)

          ! extrapolate both components of velocity to right face
          ury(i,j,1) = u(i,j  ,1) - (HALF + dt2*min(ZERO,u(i,j  ,2)/hy))*slopey(i,j  ,1)
          ury(i,j,2) = u(i,j  ,2) - (HALF + dt2*min(ZERO,u(i,j  ,2)/hy))*slopey(i,j  ,2)

          ! add source terms
          if(use_minion) then
             uly(i,j,1) = uly(i,j,1) + dt2*force(i,j-1,1)
             uly(i,j,2) = uly(i,j,2) + dt2*force(i,j-1,2)
             ury(i,j,1) = ury(i,j,1) + dt2*force(i,j  ,1)
             ury(i,j,2) = ury(i,j,2) + dt2*force(i,j  ,2)
          endif

          ! impose lo side bc's
          if(j .eq. js) then
             uly(i,j,1) = merge(u(i,js-1,1),uly(i,j,1),phys_bc(2,1) .eq. INLET)
             ury(i,j,1) = merge(u(i,js-1,1),ury(i,j,1),phys_bc(2,1) .eq. INLET)
             uly(i,j,2) = merge(u(i,js-1,2),uly(i,j,2),phys_bc(2,1) .eq. INLET)
             ury(i,j,2) = merge(u(i,js-1,2),ury(i,j,2),phys_bc(2,1) .eq. INLET)
             if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                uly(i,j,1) = merge(ZERO,ury(i,j,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                ury(i,j,1) = merge(ZERO,ury(i,j,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                uly(i,j,2) = ZERO
                ury(i,j,2) = ZERO
             endif
          endif

          ! impose hi side bc's
          if(j .eq. je+1) then
             uly(i,j,1) = merge(u(i,je+1,1),uly(i,j,1),phys_bc(2,2) .eq. INLET)
             ury(i,j,1) = merge(u(i,je+1,1),ury(i,j,1),phys_bc(2,2) .eq. INLET)
             uly(i,j,2) = merge(u(i,je+1,2),uly(i,j,2),phys_bc(2,2) .eq. INLET)
             ury(i,j,2) = merge(u(i,je+1,2),ury(i,j,2),phys_bc(2,2) .eq. INLET)
             if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                uly(i,j,1) = merge(ZERO,uly(i,j,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                ury(i,j,1) = merge(ZERO,uly(i,j,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                uly(i,j,2) = ZERO
                ury(i,j,2) = ZERO
             endif
          endif

          ! make normal component of uimhx by first solving a normal Riemann problem
          uavg = HALF*(uly(i,j,2)+ury(i,j,2))
          test = ((uly(i,j,2) .le. ZERO .and. ury(i,j,2) .ge. ZERO) .or. &
               (abs(uly(i,j,2)+ury(i,j,2)) .lt. eps))
          uimhy(i,j,2) = merge(uly(i,j,2),ury(i,j,2),uavg .gt. ZERO)
          uimhy(i,j,2) = merge(ZERO,uimhy(i,j,2),test)

          ! now upwind to get transverse component of uimhy
          uimhy(i,j,1) = merge(uly(i,j,1),ury(i,j,1),uimhy(i,j,2).gt.ZERO)
          uavg = HALF*(uly(i,j,1)+ury(i,j,1))
          uimhy(i,j,1) = merge(uavg,uimhy(i,j,1),abs(uimhy(i,j,2)).lt.eps)
       enddo
    enddo

    !******************************************************************
    ! Create umac and vmac
    !******************************************************************

    do j=js,je
       do i=is,ie+1
          ! extrapolate to edges
          umacl(i,j) = ulx(i,j,1) &
               - (dt4/hy)*(uimhy(i-1,j+1,2)+uimhy(i-1,j,2))*(uimhy(i-1,j+1,1)-uimhy(i-1,j,1))
          umacr(i,j) = urx(i,j,1) &
               - (dt4/hy)*(uimhy(i  ,j+1,2)+uimhy(i  ,j,2))*(uimhy(i  ,j+1,1)-uimhy(i  ,j,1))

          ! if use_minion is true, we have already accounted for source terms
          ! in ulx and urx; otherwise, we need to account for them here.
          if(.not. use_minion) then
             umacl(i,j) = umacl(i,j) + dt2*force(i-1,j,1)
             umacr(i,j) = umacr(i,j) + dt2*force(i  ,j,1)
          endif

          ! solve Riemann problem
          uavg = HALF*(umacl(i,j)+umacr(i,j))
          test = ((umacl(i,j) .le. ZERO .and. umacr(i,j) .ge. ZERO) .or. &
               (abs(umacl(i,j)+umacr(i,j)) .lt. eps))
          umac(i,j) = merge(umacl(i,j),umacr(i,j),uavg .gt. ZERO)
          umac(i,j) = merge(ZERO,umac(i,j),test)
       enddo
    enddo

    ! Apply boundary conditions
    do j=js,je
       ! lo side
       if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
          umac(is,j) = ZERO
       elseif (phys_bc(1,1) .eq. INLET) then
          umac(is,j) = u(is-1,j,1)
       elseif (phys_bc(1,1) .eq. OUTLET) then
          umac(is,j) = min(umacr(is,j),ZERO)
       endif

       ! hi side
       if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
          umac(ie+1,j) = ZERO
       elseif (phys_bc(1,2) .eq. INLET) then
          umac(ie+1,j) = u(ie+1,j,1)
       elseif (phys_bc(1,2) .eq. OUTLET) then
          umac(ie+1,j) = max(umacl(ie+1,j),ZERO)
       endif
    enddo

    do j=js,je+1
       do i=is,ie
          ! extrapolate to edges
          vmacl(i,j) = uly(i,j,2) &
               - (dt4/hx)*(uimhx(i+1,j-1,1)+uimhx(i,j-1,1))*(uimhx(i+1,j-1,2)-uimhx(i,j-1,2))
          vmacr(i,j) = ury(i,j,2) &
               - (dt4/hx)*(uimhx(i+1,j  ,1)+uimhx(i,j  ,1))*(uimhx(i+1,j  ,2)-uimhx(i,j  ,2))

          ! if use_minion is true, we have already accounted for source terms
          ! in uly and ury; otherwise, we need to account for them here.
          if(.not. use_minion) then
             vmacl(i,j) = vmacl(i,j) + dt2*force(i,j-1,2)
             vmacr(i,j) = vmacr(i,j) + dt2*force(i,j  ,2)
          endif

          ! solve Riemann problem
          uavg = HALF*(vmacl(i,j)+vmacr(i,j))
          test = ((vmacl(i,j) .le. ZERO .and. vmacr(i,j) .ge. ZERO) .or. &
               (abs(vmacl(i,j)+vmacr(i,j)) .lt. eps))
          vmac(i,j) = merge(vmacl(i,j),vmacr(i,j),uavg .gt. ZERO)
          vmac(i,j) = merge(ZERO,vmac(i,j),test)
       enddo
    enddo

    ! Apply boundary conditions
    do i=is,ie
       ! lo side
       if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          vmac(i,js) = ZERO
       elseif (phys_bc(2,1) .eq. INLET) then
          vmac(i,js) = u(i,js-1,2)
       elseif (phys_bc(2,1) .eq. OUTLET) then
          vmac(i,js) = min(vmacr(i,js),ZERO)
       endif

       ! hi side
       if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          vmac(i,je+1) = ZERO
       elseif (phys_bc(2,2) .eq. INLET) then
          vmac(i,je+1) = u(i,je+1,2)
       elseif (phys_bc(2,2) .eq. OUTLET) then
          vmac(i,je+1) = max(vmacl(i,je+1),ZERO)
       endif
    enddo

    deallocate(slopex)
    deallocate(slopey)

    deallocate(ulx)
    deallocate(urx)
    deallocate(uly)
    deallocate(ury)
    deallocate(uimhx)
    deallocate(uimhy)

    deallocate(umacl)
    deallocate(umacr)
    deallocate(vmacl)
    deallocate(vmacr)

  end subroutine velpred_debug_2d


  subroutine velpred_debug_3d(u, umac,vmac,wmac,force,lo,dx,dt,phys_bc,adv_bc,ng,use_minion)

    use bc_module
    use slope_module
    use bl_constants_module

    integer         ,intent(in) :: lo(3)
    integer         ,intent(in) :: ng

    real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
    real(kind=dp_t),intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)

    real(kind=dp_t) ,intent(in) :: dx(:),dt
    integer         ,intent(in) :: phys_bc(:,:)
    integer         ,intent(in) :: adv_bc(:,:,:)
    logical         ,intent(in) :: use_minion

    ! Local variables
    real(kind=dp_t), allocatable::  slopex(:,:,:,:)
    real(kind=dp_t), allocatable::  slopey(:,:,:,:)
    real(kind=dp_t), allocatable::  slopez(:,:,:,:)

    real(kind=dp_t) hx, hy, hz, dt2, dt4, dt6, uavg

    integer :: hi(3)
    integer :: slope_order = 4
    logical :: test

    real(kind=dp_t) :: abs_eps, eps, umax
    integer :: i,j,k,is,js,ks,ie,je,ke

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable:: ulx(:,:,:,:),urx(:,:,:,:),uimhx(:,:,:,:)
    real(kind=dp_t), allocatable:: uly(:,:,:,:),ury(:,:,:,:),uimhy(:,:,:,:)
    real(kind=dp_t), allocatable:: ulz(:,:,:,:),urz(:,:,:,:),uimhz(:,:,:,:)

    ! these correspond to u_L^{y|z}, etc.
    real(kind=dp_t), allocatable:: ulyz(:,:,:)
    real(kind=dp_t), allocatable:: uryz(:,:,:)
    real(kind=dp_t), allocatable:: uimhyz(:,:,:)

    real(kind=dp_t), allocatable:: ulzy(:,:,:)
    real(kind=dp_t), allocatable:: urzy(:,:,:)
    real(kind=dp_t), allocatable:: uimhzy(:,:,:)

    real(kind=dp_t), allocatable:: vlxz(:,:,:)
    real(kind=dp_t), allocatable:: vrxz(:,:,:)
    real(kind=dp_t), allocatable:: vimhxz(:,:,:)

    real(kind=dp_t), allocatable:: vlzx(:,:,:)
    real(kind=dp_t), allocatable:: vrzx(:,:,:)
    real(kind=dp_t), allocatable:: vimhzx(:,:,:)

    real(kind=dp_t), allocatable:: wlxy(:,:,:)
    real(kind=dp_t), allocatable:: wrxy(:,:,:)
    real(kind=dp_t), allocatable:: wimhxy(:,:,:)

    real(kind=dp_t), allocatable:: wlyx(:,:,:)
    real(kind=dp_t), allocatable:: wryx(:,:,:)
    real(kind=dp_t), allocatable:: wimhyx(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable:: umacl(:,:,:),umacr(:,:,:)
    real(kind=dp_t), allocatable:: vmacl(:,:,:),vmacr(:,:,:)
    real(kind=dp_t), allocatable:: wmacl(:,:,:),wmacr(:,:,:)

    hi(1) = lo(1) + size(u,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(u,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(u,dim=3) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(ulz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(urz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(uimhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uryz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uimhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    allocate(ulzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(urzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(uimhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    allocate(vlxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vrxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vimhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    allocate(vlzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vrzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vimhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    allocate(wlxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wrxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wimhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    allocate(wlyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wryx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wimhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(wmacl(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(wmacr(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

    do k = lo(3)-1,hi(3)+1
       call slopex_2d(u(:,:,k,:),slopex(:,:,k,:),lo,ng,3,adv_bc,slope_order)
       call slopey_2d(u(:,:,k,:),slopey(:,:,k,:),lo,ng,3,adv_bc,slope_order)
    end do
    call slopez_3d(u,slopez,lo,ng,3,adv_bc,slope_order)

    abs_eps = 1.0d-8

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

    ! Compute eps, which is relative to the max velocity
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
    if(umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate all components of velocity to left face
             ulx(i,j,k,1) = u(i-1,j,k,1) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,1)
             ulx(i,j,k,2) = u(i-1,j,k,2) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,2)
             ulx(i,j,k,3) = u(i-1,j,k,3) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,3)

             ! extrapolate all components of velocity to right face
             urx(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,1)
             urx(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,2)
             urx(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,3)

             ! add source terms
             if(use_minion) then
                ulx(i,j,k,1) = ulx(i,j,k,1) + dt2*force(i-1,j,k,1)
                ulx(i,j,k,2) = ulx(i,j,k,2) + dt2*force(i-1,j,k,2)
                ulx(i,j,k,3) = ulx(i,j,k,3) + dt2*force(i-1,j,k,3)
                urx(i,j,k,1) = urx(i,j,k,1) + dt2*force(i  ,j,k,1)
                urx(i,j,k,2) = urx(i,j,k,2) + dt2*force(i  ,j,k,2)
                urx(i,j,k,3) = urx(i,j,k,3) + dt2*force(i  ,j,k,3)
             endif

             ! impose lo side bc's
             if(i .eq. is) then
                ulx(i,j,k,1) = merge(u(is-1,j,k,1),ulx(i,j,k,1),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,1) = merge(u(is-1,j,k,1),urx(i,j,k,1),phys_bc(1,1) .eq. INLET)
                ulx(i,j,k,2) = merge(u(is-1,j,k,2),ulx(i,j,k,2),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,2) = merge(u(is-1,j,k,2),urx(i,j,k,2),phys_bc(1,1) .eq. INLET)
                ulx(i,j,k,3) = merge(u(is-1,j,k,3),ulx(i,j,k,3),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,3) = merge(u(is-1,j,k,3),urx(i,j,k,3),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   ulx(i,j,k,1) = ZERO
                   urx(i,j,k,1) = ZERO
                   ulx(i,j,k,2) = merge(ZERO,urx(i,j,k,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,k,2) = merge(ZERO,urx(i,j,k,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   ulx(i,j,k,3) = merge(ZERO,urx(i,j,k,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,k,3) = merge(ZERO,urx(i,j,k,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                ulx(i,j,k,1) = merge(u(ie+1,j,k,1),ulx(i,j,k,1),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,1) = merge(u(ie+1,j,k,1),urx(i,j,k,1),phys_bc(1,2) .eq. INLET)
                ulx(i,j,k,2) = merge(u(ie+1,j,k,2),ulx(i,j,k,2),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,2) = merge(u(ie+1,j,k,2),urx(i,j,k,2),phys_bc(1,2) .eq. INLET)
                ulx(i,j,k,3) = merge(u(ie+1,j,k,3),ulx(i,j,k,3),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,3) = merge(u(ie+1,j,k,3),urx(i,j,k,3),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   ulx(i,j,k,1) = ZERO
                   urx(i,j,k,1) = ZERO
                   ulx(i,j,k,2) = merge(ZERO,ulx(i,j,k,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,k,2) = merge(ZERO,ulx(i,j,k,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   ulx(i,j,k,3) = merge(ZERO,ulx(i,j,k,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,k,3) = merge(ZERO,ulx(i,j,k,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhx by first solving a normal Riemann problem
             uavg = HALF*(ulx(i,j,k,1)+urx(i,j,k,1))
             test = ((ulx(i,j,k,1) .le. ZERO .and. urx(i,j,k,1) .ge. ZERO) .or. &
                  (abs(ulx(i,j,k,1)+urx(i,j,k,1)) .lt. eps))
             uimhx(i,j,k,1) = merge(ulx(i,j,k,1),urx(i,j,k,1),uavg .gt. ZERO)
             uimhx(i,j,k,1) = merge(ZERO,uimhx(i,j,k,1),test)

             ! now upwind to get transverse components of uimhx
             uimhx(i,j,k,2) = merge(ulx(i,j,k,2),urx(i,j,k,2),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,2)+urx(i,j,k,2))
             uimhx(i,j,k,2) = merge(uavg,uimhx(i,j,k,2),abs(uimhx(i,j,k,1)).lt.eps)

             uimhx(i,j,k,3) = merge(ulx(i,j,k,3),urx(i,j,k,3),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,3)+urx(i,j,k,3))
             uimhx(i,j,k,3) = merge(uavg,uimhx(i,j,k,3),abs(uimhx(i,j,k,1)).lt.eps)
          enddo
       enddo
    enddo

    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate all components of velocity to left face
             uly(i,j,k,1) = u(i,j-1,k,1) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,1)
             uly(i,j,k,2) = u(i,j-1,k,2) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,2)
             uly(i,j,k,3) = u(i,j-1,k,3) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,3)

             ! extrapolate all components of velocity to right face
             ury(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,1)
             ury(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,2)
             ury(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,3)

             ! add source terms
             if(use_minion) then
                uly(i,j,k,1) = uly(i,j,k,1) + dt2*force(i,j-1,k,1)
                uly(i,j,k,2) = uly(i,j,k,2) + dt2*force(i,j-1,k,2)
                uly(i,j,k,3) = uly(i,j,k,3) + dt2*force(i,j-1,k,3)
                ury(i,j,k,1) = ury(i,j,k,1) + dt2*force(i,j  ,k,1)
                ury(i,j,k,2) = ury(i,j,k,2) + dt2*force(i,j  ,k,2)
                ury(i,j,k,3) = ury(i,j,k,3) + dt2*force(i,j  ,k,3)
             endif

             ! impose lo side bc's
             if(j .eq. js) then
                uly(i,j,k,1) = merge(u(i,js-1,k,1),uly(i,j,k,1),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,1) = merge(u(i,js-1,k,1),ury(i,j,k,1),phys_bc(2,1) .eq. INLET)
                uly(i,j,k,2) = merge(u(i,js-1,k,2),uly(i,j,k,2),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,2) = merge(u(i,js-1,k,2),ury(i,j,k,2),phys_bc(2,1) .eq. INLET)
                uly(i,j,k,3) = merge(u(i,js-1,k,3),uly(i,j,k,3),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,3) = merge(u(i,js-1,k,3),ury(i,j,k,3),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   uly(i,j,k,1) = merge(ZERO,ury(i,j,k,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,k,1) = merge(ZERO,ury(i,j,k,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uly(i,j,k,2) = ZERO
                   ury(i,j,k,2) = ZERO
                   uly(i,j,k,3) = merge(ZERO,ury(i,j,k,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,k,3) = merge(ZERO,ury(i,j,k,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                uly(i,j,k,1) = merge(u(i,je+1,k,1),uly(i,j,k,1),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,1) = merge(u(i,je+1,k,1),ury(i,j,k,1),phys_bc(2,2) .eq. INLET)
                uly(i,j,k,2) = merge(u(i,je+1,k,2),uly(i,j,k,2),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,2) = merge(u(i,je+1,k,2),ury(i,j,k,2),phys_bc(2,2) .eq. INLET)
                uly(i,j,k,3) = merge(u(i,je+1,k,3),uly(i,j,k,3),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,3) = merge(u(i,je+1,k,3),ury(i,j,k,3),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   uly(i,j,k,1) = merge(ZERO,uly(i,j,k,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,k,1) = merge(ZERO,uly(i,j,k,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uly(i,j,k,2) = ZERO
                   ury(i,j,k,2) = ZERO
                   uly(i,j,k,3) = merge(ZERO,uly(i,j,k,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,k,3) = merge(ZERO,uly(i,j,k,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhy by first solving a normal Riemann problem
             uavg = HALF*(uly(i,j,k,2)+ury(i,j,k,2))
             test = ((uly(i,j,k,2) .le. ZERO .and. ury(i,j,k,2) .ge. ZERO) .or. &
                  (abs(uly(i,j,k,2)+ury(i,j,k,2)) .lt. eps))
             uimhy(i,j,k,2) = merge(uly(i,j,k,2),ury(i,j,k,2),uavg .gt. ZERO)
             uimhy(i,j,k,2) = merge(ZERO,uimhy(i,j,k,2),test)

             ! now upwind to get transverse components of uimhy
             uimhy(i,j,k,1) = merge(uly(i,j,k,1),ury(i,j,k,1),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(uly(i,j,k,1)+ury(i,j,k,1))
             uimhy(i,j,k,1) = merge(uavg,uimhy(i,j,k,1),abs(uimhy(i,j,k,2)).lt.eps)

             uimhy(i,j,k,3) = merge(uly(i,j,k,3),ury(i,j,k,3),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(uly(i,j,k,3)+ury(i,j,k,3))
             uimhy(i,j,k,3) = merge(uavg,uimhy(i,j,k,3),abs(uimhy(i,j,k,2)).lt.eps)
          enddo
       enddo
    enddo

    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! extrapolate all components of velocity to left face
             ulz(i,j,k,1) = u(i,j,k-1,1) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,1)
             ulz(i,j,k,2) = u(i,j,k-1,2) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,2)
             ulz(i,j,k,3) = u(i,j,k-1,3) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,3)

             ! extrapolate all components of velocity to right face
             urz(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,1)
             urz(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,2)
             urz(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,3)

             ! add source terms
             if(use_minion) then
                ulz(i,j,k,1) = ulz(i,j,k,1) + dt2*force(i,j,k-1,1)
                ulz(i,j,k,2) = ulz(i,j,k,2) + dt2*force(i,j,k-1,2)
                ulz(i,j,k,3) = ulz(i,j,k,3) + dt2*force(i,j,k-1,3)
                urz(i,j,k,1) = urz(i,j,k,1) + dt2*force(i,j,k  ,1)
                urz(i,j,k,2) = urz(i,j,k,2) + dt2*force(i,j,k  ,2)
                urz(i,j,k,3) = urz(i,j,k,3) + dt2*force(i,j,k  ,3)
             endif

             ! impose lo side bc's
             if(k .eq. ks) then
                ulz(i,j,k,1) = merge(u(i,j,ks-1,1),ulz(i,j,k,1),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,1) = merge(u(i,j,ks-1,1),urz(i,j,k,1),phys_bc(3,1) .eq. INLET)
                ulz(i,j,k,2) = merge(u(i,j,ks-1,2),ulz(i,j,k,2),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,2) = merge(u(i,j,ks-1,2),urz(i,j,k,2),phys_bc(3,1) .eq. INLET)
                ulz(i,j,k,3) = merge(u(i,j,ks-1,3),ulz(i,j,k,3),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,3) = merge(u(i,j,ks-1,3),urz(i,j,k,3),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   ulz(i,j,k,1) = merge(ZERO,urz(i,j,k,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urz(i,j,k,1) = merge(ZERO,urz(i,j,k,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,2) = merge(ZERO,urz(i,j,k,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urz(i,j,k,2) = merge(ZERO,urz(i,j,k,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,3) = ZERO
                   urz(i,j,k,3) = ZERO
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                ulz(i,j,k,1) = merge(u(i,j,ke+1,1),ulz(i,j,k,1),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,1) = merge(u(i,j,ke+1,1),urz(i,j,k,1),phys_bc(3,2) .eq. INLET)
                ulz(i,j,k,2) = merge(u(i,j,ke+1,2),ulz(i,j,k,2),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,2) = merge(u(i,j,ke+1,2),urz(i,j,k,2),phys_bc(3,2) .eq. INLET)
                ulz(i,j,k,3) = merge(u(i,j,ke+1,3),ulz(i,j,k,3),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,3) = merge(u(i,j,ke+1,3),urz(i,j,k,3),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   ulz(i,j,k,1) = merge(ZERO,ulz(i,j,k,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urz(i,j,k,1) = merge(ZERO,ulz(i,j,k,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,2) = merge(ZERO,ulz(i,j,k,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urz(i,j,k,2) = merge(ZERO,ulz(i,j,k,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,3) = ZERO
                   urz(i,j,k,3) = ZERO
                endif
             endif

             ! make normal component of uimhz by first solving a normal Riemann problem
             uavg = HALF*(ulz(i,j,k,3)+urz(i,j,k,3))
             test = ((ulz(i,j,k,3) .le. ZERO .and. urz(i,j,k,3) .ge. ZERO) .or. &
                  (abs(ulz(i,j,k,3)+urz(i,j,k,3)) .lt. eps))
             uimhz(i,j,k,3) = merge(ulz(i,j,k,3),urz(i,j,k,3),uavg .gt. ZERO)
             uimhz(i,j,k,3) = merge(ZERO,uimhz(i,j,k,3),test)

             ! now upwind to get transverse components of uimhz
             uimhz(i,j,k,1) = merge(ulz(i,j,k,1),urz(i,j,k,1),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,1)+urz(i,j,k,1))
             uimhz(i,j,k,1) = merge(uavg,uimhz(i,j,k,1),abs(uimhz(i,j,k,3)).lt.eps)

             uimhz(i,j,k,2) = merge(ulz(i,j,k,2),urz(i,j,k,2),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,2)+urz(i,j,k,2))
             uimhz(i,j,k,2) = merge(uavg,uimhz(i,j,k,2),abs(uimhz(i,j,k,3)).lt.eps)
          enddo
       enddo
    enddo

    !******************************************************************
    ! Create u_{\i-\half\e_y}^{y|z}, etc.
    !******************************************************************

    ! uimhyz loop
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate to faces
             ulyz(i,j,k) = uly(i,j,k,1) &
                  - (dt6/hz)*(uimhz(i,j-1,k+1,3)+uimhz(i,j-1,k,3))*(uimhz(i,j-1,k+1,1)-uimhz(i,j-1,k,1))               
             uryz(i,j,k) = ury(i,j,k,1) &
                  - (dt6/hz)*(uimhz(i,j  ,k+1,3)+uimhz(i,j  ,k,3))*(uimhz(i,j  ,k+1,1)-uimhz(i,j  ,k,1)) 

             ! impose lo side bc's
             if(j .eq. js) then
                ulyz(i,j,k) = merge(u(i,js-1,k,1),ulyz(i,j,k),phys_bc(2,1) .eq. INLET)
                uryz(i,j,k) = merge(u(i,js-1,k,1),uryz(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   ulyz(i,j,k) = merge(ZERO,uryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uryz(i,j,k) = merge(ZERO,uryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                ulyz(i,j,k) = merge(u(i,je+1,k,1),ulyz(i,j,k),phys_bc(2,2) .eq. INLET)
                uryz(i,j,k) = merge(u(i,je+1,k,1),uryz(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   ulyz(i,j,k) = merge(ZERO,ulyz(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uryz(i,j,k) = merge(ZERO,ulyz(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             uimhyz(i,j,k) = merge(ulyz(i,j,k),uryz(i,j,k),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(ulyz(i,j,k)+uryz(i,j,k))
             uimhyz(i,j,k) = merge(uavg,uimhyz(i,j,k),abs(uimhy(i,j,k,2)).lt.eps)
          enddo
       enddo
    enddo

    ! uimhzy loop
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! extrapolate to faces
             ulzy(i,j,k) = ulz(i,j,k,1) &
                  - (dt6/hy)*(uimhy(i,j+1,k-1,2)+uimhy(i,j,k-1,2))*(uimhy(i,j+1,k-1,1)-uimhy(i,j,k-1,1))
             urzy(i,j,k) = urz(i,j,k,1) &
                  - (dt6/hy)*(uimhy(i,j+1,k  ,2)+uimhy(i,j,k  ,2))*(uimhy(i,j+1,k  ,1)-uimhy(i,j,k  ,1))

             ! impose lo side bc's
             if(k .eq. ks) then
                ulzy(i,j,k) = merge(u(i,j,ks-1,1),ulzy(i,j,k),phys_bc(3,1) .eq. INLET)
                urzy(i,j,k) = merge(u(i,j,ks-1,1),urzy(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   ulzy(i,j,k) = merge(ZERO,urzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urzy(i,j,k) = merge(ZERO,urzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                ulzy(i,j,k) = merge(u(i,j,ke+1,1),ulzy(i,j,k),phys_bc(3,2) .eq. INLET)
                urzy(i,j,k) = merge(u(i,j,ke+1,1),urzy(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   ulzy(i,j,k) = merge(ZERO,ulzy(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urzy(i,j,k) = merge(ZERO,ulzy(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             uimhzy(i,j,k) = merge(ulzy(i,j,k),urzy(i,j,k),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulzy(i,j,k)+urzy(i,j,k))
             uimhzy(i,j,k) = merge(uavg,uimhzy(i,j,k),abs(uimhz(i,j,k,3)).lt.eps)
          enddo
       enddo
    enddo

    ! vimhxz loop
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate to faces
             vlxz(i,j,k) = ulx(i,j,k,2) &
                  - (dt6/hz)*(uimhz(i-1,j,k+1,3)+uimhz(i-1,j,k,3))*(uimhz(i-1,j,k+1,2)-uimhz(i-1,j,k,2))
             vrxz(i,j,k) = urx(i,j,k,2) &
                  - (dt6/hz)*(uimhz(i  ,j,k+1,3)+uimhz(i  ,j,k,3))*(uimhz(i  ,j,k+1,2)-uimhz(i  ,j,k,2))

             ! impose lo side bc's
             if(i .eq. is) then
                vlxz(i,j,k) = merge(u(is-1,j,k,2),vlxz(i,j,k),phys_bc(1,1) .eq. INLET)
                vrxz(i,j,k) = merge(u(is-1,j,k,2),vrxz(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   vlxz(i,j,k) = merge(ZERO,vrxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   vrxz(i,j,k) = merge(ZERO,vrxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                vlxz(i,j,k) = merge(u(ie+1,j,k,2),vlxz(i,j,k),phys_bc(1,2) .eq. INLET)
                vrxz(i,j,k) = merge(u(ie+1,j,k,2),vrxz(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   vlxz(i,j,k) = merge(ZERO,vlxz(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   vrxz(i,j,k) = merge(ZERO,vlxz(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             vimhxz(i,j,k) = merge(vlxz(i,j,k),vrxz(i,j,k),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(vlxz(i,j,k)+vrxz(i,j,k))
             vimhxz(i,j,k) = merge(uavg,vimhxz(i,j,k),abs(uimhx(i,j,k,1)).lt.eps)
          enddo
       enddo
    enddo

    ! vimhzx loop
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! extrapolate to faces
             vlzx(i,j,k) = ulz(i,j,k,2) &
                  - (dt6/hx)*(uimhx(i+1,j,k-1,1)+uimhx(i,j,k-1,1))*(uimhx(i+1,j,k-1,2)-uimhx(i,j,k-1,2))
             vrzx(i,j,k) = urz(i,j,k,2) &
                  - (dt6/hx)*(uimhx(i+1,j,k  ,1)+uimhx(i,j,k  ,1))*(uimhx(i+1,j,k  ,2)-uimhx(i,j,k  ,2))

             ! impose lo side bc's
             if(k .eq. ks) then
                vlzx(i,j,k) = merge(u(i,j,ks-1,1),vlzx(i,j,k),phys_bc(3,1) .eq. INLET)
                vrzx(i,j,k) = merge(u(i,j,ks-1,1),vrzx(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   vlzx(i,j,k) = merge(ZERO,vrzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   vrzx(i,j,k) = merge(ZERO,vrzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                vlzx(i,j,k) = merge(u(i,j,ke+1,1),vlzx(i,j,k),phys_bc(3,2) .eq. INLET)
                vrzx(i,j,k) = merge(u(i,j,ke+1,1),vrzx(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   vlzx(i,j,k) = merge(ZERO,vlzx(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   vrzx(i,j,k) = merge(ZERO,vlzx(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             vimhzx(i,j,k) = merge(vlzx(i,j,k),vrzx(i,j,k),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(vlzx(i,j,k)+vrzx(i,j,k))
             vimhzx(i,j,k) = merge(uavg,vimhzx(i,j,k),abs(uimhz(i,j,k,3)).lt.eps)
          enddo
       enddo
    enddo

    ! wimhxy loop
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! extrapolate to faces
             wlxy(i,j,k) = ulx(i,j,k,3) &
                  - (dt6/hy)*(uimhy(i-1,j+1,k,2)+uimhy(i-1,j,k,2))*(uimhy(i-1,j+1,k,3)-uimhy(i-1,j,k,3))
             wrxy(i,j,k) = urx(i,j,k,3) &
                  - (dt6/hy)*(uimhy(i  ,j+1,k,2)+uimhy(i  ,j,k,2))*(uimhy(i  ,j+1,k,3)-uimhy(i  ,j,k,3))

             ! impose lo side bc's
             if(i .eq. is) then
                wlxy(i,j,k) = merge(u(is-1,j,k,3),wlxy(i,j,k),phys_bc(1,1) .eq. INLET)
                wrxy(i,j,k) = merge(u(is-1,j,k,3),wrxy(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,k) = merge(ZERO,wrxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   wrxy(i,j,k) = merge(ZERO,wrxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                wlxy(i,j,k) = merge(u(ie+1,j,k,3),wlxy(i,j,k),phys_bc(1,2) .eq. INLET)
                wrxy(i,j,k) = merge(u(ie+1,j,k,3),wrxy(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,k) = merge(ZERO,wlxy(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   wrxy(i,j,k) = merge(ZERO,wlxy(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhxy(i,j,k) = merge(wlxy(i,j,k),wrxy(i,j,k),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(wlxy(i,j,k)+wrxy(i,j,k))
             wimhxy(i,j,k) = merge(uavg,wimhxy(i,j,k),abs(uimhx(i,j,k,1)).lt.eps)
          enddo
       enddo
    enddo

    ! wimhyx loop
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! extrapolate to faces
             wlyx(i,j,k) = uly(i,j,k,3) &
                  - (dt6/hx)*(uimhx(i+1,j-1,k,1)+uimhx(i,j-1,k,1))*(uimhx(i+1,j-1,k,3)-uimhx(i,j-1,k,3))
             wryx(i,j,k) = ury(i,j,k,3) &
                  - (dt6/hx)*(uimhx(i+1,j  ,k,1)+uimhx(i,j  ,k,1))*(uimhx(i+1,j  ,k,3)-uimhx(i,j  ,k,3))

             ! impose lo side bc's
             if(j .eq. js) then
                wlyx(i,j,k) = merge(u(i,js-1,k,3),wlyx(i,j,k),phys_bc(2,1) .eq. INLET)
                wryx(i,j,k) = merge(u(i,js-1,k,3),wryx(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,k) = merge(ZERO,wryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   wryx(i,j,k) = merge(ZERO,wryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                wlyx(i,j,k) = merge(u(i,je+1,k,3),wlyx(i,j,k),phys_bc(2,2) .eq. INLET)
                wryx(i,j,k) = merge(u(i,je+1,k,3),wryx(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,k) = merge(ZERO,wlyx(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   wryx(i,j,k) = merge(ZERO,wlyx(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhyx(i,j,k) = merge(wlyx(i,j,k),wryx(i,j,k),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(wlyx(i,j,k)+wryx(i,j,k))
             wimhyx(i,j,k) = merge(uavg,wimhyx(i,j,k),abs(uimhy(i,j,k,2)).lt.eps)
          enddo
       enddo
    enddo


    !******************************************************************
    ! Create umac, etc.
    !******************************************************************

    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             umacl(i,j,k) = ulx(i,j,k,1) &
                  - (dt4/hy)*(uimhy(i-1,j+1,k  ,2)+uimhy(i-1,j,k,2))*(uimhyz(i-1,j+1,k  )-uimhyz(i-1,j,k)) &
                  - (dt4/hz)*(uimhz(i-1,j  ,k+1,3)+uimhz(i-1,j,k,3))*(uimhzy(i-1,j  ,k+1)-uimhzy(i-1,j,k))
             umacr(i,j,k) = urx(i,j,k,1) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k  ,2)+uimhy(i  ,j,k,2))*(uimhyz(i  ,j+1,k  )-uimhyz(i  ,j,k)) &
                  - (dt4/hz)*(uimhz(i  ,j  ,k+1,3)+uimhz(i  ,j,k,3))*(uimhzy(i  ,j  ,k+1)-uimhzy(i  ,j,k))

             ! if use_minion is true, we have already accounted for source terms
             ! in ulx and urx; otherwise, we need to account for them here.
             if(.not. use_minion) then
                umacl(i,j,k) = umacl(i,j,k) + dt2*force(i-1,j,k,1)
                umacr(i,j,k) = umacr(i,j,k) + dt2*force(i  ,j,k,1)
             endif

             ! solve Riemann problem
             uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
             test = ((umacl(i,j,k) .le. ZERO .and. umacr(i,j,k) .ge. ZERO) .or. &
                  (abs(umacl(i,j,k)+umacr(i,j,k)) .lt. eps))
             umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg .gt. ZERO)
             umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do k=ks,ke
       do j=js,je
          ! lo side
          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             umac(is,j,k) = ZERO
          elseif (phys_bc(1,1) .eq. INLET) then
             umac(is,j,k) = u(is-1,j,k,1)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             umac(is,j,k) = min(umacr(is,j,k),ZERO)
          endif

          ! hi side
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             umac(ie+1,j,k) = ZERO
          elseif (phys_bc(1,2) .eq. INLET) then
             umac(ie+1,j,k) = u(ie+1,j,k,1)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             umac(ie+1,j,k) = max(umacl(ie+1,j,k),ZERO)
          endif
       enddo
    enddo

    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vmacl(i,j,k) = uly(i,j,k,2) &
                  - (dt4/hx)*(uimhx(i+1,j-1,k  ,1)+uimhx(i,j-1,k,1))*(vimhxz(i+1,j-1,k  )-vimhxz(i,j-1,k)) &
                  - (dt4/hz)*(uimhz(i  ,j-1,k+1,3)+uimhz(i,j-1,k,3))*(vimhzx(i  ,j-1,k+1)-vimhzx(i,j-1,k))
             vmacr(i,j,k) = ury(i,j,k,2) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k  ,1)+uimhx(i,j  ,k,1))*(vimhxz(i+1,j  ,k  )-vimhxz(i,j  ,k)) &
                  - (dt4/hz)*(uimhz(i  ,j  ,k+1,3)+uimhz(i,j  ,k,3))*(vimhzx(i  ,j  ,k+1)-vimhzx(i,j  ,k))

             ! if use_minion is true, we have already accounted for source terms
             ! in uly and ury; otherwise, we need to account for them here.
             if(.not. use_minion) then
                vmacl(i,j,k) = vmacl(i,j,k) + dt2*force(i,j-1,k,2)
                vmacr(i,j,k) = vmacr(i,j,k) + dt2*force(i,j  ,k,2)
             endif

             ! solve Riemann problem
             uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
             test = ((vmacl(i,j,k) .le. ZERO .and. vmacr(i,j,k) .ge. ZERO) .or. &
                  (abs(vmacl(i,j,k)+vmacr(i,j,k)) .lt. eps))
             vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg .gt. ZERO)
             vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do k=ks,ke
       do i=is,ie
          ! lo side
          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             vmac(i,js,k) = ZERO
          elseif (phys_bc(2,1) .eq. INLET) then
             vmac(i,js,k) = u(i,js-1,k,2)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             vmac(i,js,k) = min(vmacr(i,js,k),ZERO)
          endif

          ! hi side
          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             vmac(i,je+1,k) = ZERO
          elseif (phys_bc(2,2) .eq. INLET) then
             vmac(i,je+1,k) = u(i,je+1,k,2)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             vmac(i,je+1,k) = max(vmacl(i,je+1,k),ZERO)
          endif
       enddo
    enddo

    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! extrapolate to edges
             wmacl(i,j,k) = ulz(i,j,k,3) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k-1,1)+uimhx(i,j,k-1,1))*(wimhxy(i+1,j  ,k-1)-wimhxy(i,j,k-1)) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k-1,2)+uimhy(i,j,k-1,2))*(wimhyx(i  ,j+1,k-1)-wimhyx(i,j,k-1))
             wmacr(i,j,k) = urz(i,j,k,3) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k  ,1)+uimhx(i,j,k  ,1))*(wimhxy(i+1,j  ,k  )-wimhxy(i,j,k  )) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k  ,2)+uimhy(i,j,k  ,2))*(wimhyx(i  ,j+1,k  )-wimhyx(i,j,k  ))

             ! if use_minion is true, we have already accounted for source terms
             ! in uly and ury; otherwise, we need to account for them here.
             if(.not. use_minion) then
                wmacl(i,j,k) = wmacl(i,j,k) + dt2*force(i,j,k-1,3)
                wmacr(i,j,k) = wmacr(i,j,k) + dt2*force(i,j,k  ,3)
             endif

             ! solve Riemann problem
             uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
             test = ((wmacl(i,j,k) .le. ZERO .and. wmacr(i,j,k) .ge. ZERO) .or. &
                  (abs(wmacl(i,j,k)+wmacr(i,j,k)) .lt. eps))
             wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg .gt. ZERO)
             wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do j=js,je
       do i=is,ie
          ! lo side
          if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
             wmac(i,j,ks) = ZERO
          elseif (phys_bc(3,1) .eq. INLET) then
             wmac(i,j,ks) = u(i,j,ks-1,3)
          elseif (phys_bc(3,1) .eq. OUTLET) then
             wmac(i,j,ks) = min(wmacr(i,j,ks),ZERO)
          endif

          ! hi side
          if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
             wmac(i,j,ke+1) = ZERO
          elseif (phys_bc(3,2) .eq. INLET) then
             wmac(i,j,ke+1) = u(i,j,ke+1,3)
          elseif (phys_bc(3,2) .eq. OUTLET) then
             wmac(i,j,ke+1) = max(wmacl(i,j,ke+1),ZERO)
          endif
       enddo
    enddo

    deallocate(slopex)
    deallocate(slopey)
    deallocate(slopez)

    deallocate(ulx)
    deallocate(urx)
    deallocate(uimhx)
    deallocate(uly)
    deallocate(ury)
    deallocate(uimhy)
    deallocate(ulz)
    deallocate(urz)
    deallocate(uimhz)

    deallocate(ulyz)
    deallocate(uryz)
    deallocate(uimhyz)

    deallocate(ulzy)
    deallocate(urzy)
    deallocate(uimhzy)

    deallocate(vlxz)
    deallocate(vrxz)
    deallocate(vimhxz)

    deallocate(vlzx)
    deallocate(vrzx)
    deallocate(vimhzx)

    deallocate(wlxy)
    deallocate(wrxy)
    deallocate(wimhxy)

    deallocate(wlyx)
    deallocate(wryx)
    deallocate(wimhyx)

    deallocate(umacl)
    deallocate(umacr)
    deallocate(vmacl)
    deallocate(vmacr)
    deallocate(wmacl)
    deallocate(wmacr)

  end subroutine velpred_debug_3d


  subroutine velpred_3d(u,umac,vmac,wmac,force,lo,dx,dt,phys_bc,adv_bc,ng,use_minion)

    use bc_module
    use slope_module
    use bl_constants_module

    integer         ,intent(in) :: lo(3)
    integer         ,intent(in) :: ng

    real(kind=dp_t),intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
    real(kind=dp_t),intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)

    real(kind=dp_t) ,intent(in) :: dx(:),dt
    integer         ,intent(in) :: phys_bc(:,:)
    integer         ,intent(in) :: adv_bc(:,:,:)
    logical         ,intent(in) :: use_minion

    ! Local variables
    real(kind=dp_t), allocatable::  slopex(:,:,:,:)
    real(kind=dp_t), allocatable::  slopey(:,:,:,:)
    real(kind=dp_t), allocatable::  slopez(:,:,:,:)

    real(kind=dp_t) hx, hy, hz, dt2, dt4, dt6, uavg

    integer :: hi(3)
    integer :: slope_order = 4
    logical :: test

    real(kind=dp_t) :: abs_eps, eps, umax
    integer :: i,j,k,is,js,ks,ie,je,ke
    integer :: kc,kp ! "current" and "previous" k

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable:: ulx(:,:,:,:),urx(:,:,:,:),uimhx(:,:,:,:)
    real(kind=dp_t), allocatable:: uly(:,:,:,:),ury(:,:,:,:),uimhy(:,:,:,:)
    real(kind=dp_t), allocatable:: ulz(:,:,:,:),urz(:,:,:,:),uimhz(:,:,:,:)

    ! these correspond to u_L^{y|z}, etc.
    real(kind=dp_t), allocatable:: ulyz(:,:,:)
    real(kind=dp_t), allocatable:: uryz(:,:,:)
    real(kind=dp_t), allocatable:: uimhyz(:,:,:)

    real(kind=dp_t), allocatable:: ulzy(:,:,:)
    real(kind=dp_t), allocatable:: urzy(:,:,:)
    real(kind=dp_t), allocatable:: uimhzy(:,:,:)

    real(kind=dp_t), allocatable:: vlxz(:,:,:)
    real(kind=dp_t), allocatable:: vrxz(:,:,:)
    real(kind=dp_t), allocatable:: vimhxz(:,:,:)

    real(kind=dp_t), allocatable:: vlzx(:,:,:)
    real(kind=dp_t), allocatable:: vrzx(:,:,:)
    real(kind=dp_t), allocatable:: vimhzx(:,:,:)

    real(kind=dp_t), allocatable:: wlxy(:,:,:)
    real(kind=dp_t), allocatable:: wrxy(:,:,:)
    real(kind=dp_t), allocatable:: wimhxy(:,:,:)

    real(kind=dp_t), allocatable:: wlyx(:,:,:)
    real(kind=dp_t), allocatable:: wryx(:,:,:)
    real(kind=dp_t), allocatable:: wimhyx(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable:: umacl(:,:),umacr(:,:)
    real(kind=dp_t), allocatable:: vmacl(:,:),vmacr(:,:)
    real(kind=dp_t), allocatable:: wmacl(:,:),wmacr(:,:)

    hi(1) = lo(1) + size(u,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(u,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(u,dim=3) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    do k = lo(3)-1,hi(3)+1
       call slopex_2d(u(:,:,k,:),slopex(:,:,k,:),lo,ng,3,adv_bc,slope_order)
       call slopey_2d(u(:,:,k,:),slopey(:,:,k,:),lo,ng,3,adv_bc,slope_order)
    end do
    call slopez_3d(u,slopez,lo,ng,3,adv_bc,slope_order)

    ! Note: All of these arrays are allocated to exactly the 
    ! size they need to be in order to compute edge states on 
    ! a domain with faces from lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1

    !***********************
    ! Normal predictor terms
    !***********************

    ! lo(1):hi(1)+1 in the x-direction
    ! lo(2)-1:hi(2)+1 in the y-direction
    ! 2 rows needed in z-direction
    allocate(ulx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2,3))
    allocate(urx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2,3))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2,3))

    ! lo(1)-1:hi(1)+1 in the x-direction
    ! lo(2):hi(2)+1 in the y-direction
    ! 2 rows needed in z-direction
    allocate(uly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2,3))
    allocate(ury  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2,3))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2,3))

    ! lo(1)-1:hi(1)+1 in the x-direction
    ! lo(2)-1:hi(2)+1 in the y-direction
    ! 2 rows needed in z-direction
    allocate(ulz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2,3))
    allocate(urz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2,3))
    allocate(uimhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2,3))

    !*****************
    ! Transverse terms
    !*****************

    ! lo(1)-1:hi(1)+1 in the x-direction
    ! lo(2):hi(2)+1 in the y-direction
    ! 2 rows needed in the z-direction
    allocate(ulyz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(uryz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(uimhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

    ! lo(1)-1:hi(1)+1 in the x-direction
    ! lo(2):hi(2) in the y-direction
    ! 2 rows needed in the z-direction
    allocate(ulzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),2))
    allocate(urzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),2))
    allocate(uimhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),2))

    ! lo(1):hi(1)+1 in the x-direction
    ! lo(2)-1:hi(2)+1 in the y-direction
    ! 2 rows needed in the z-direction
    allocate(vlxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(vrxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(vimhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

    ! lo(1):hi(1) in the x-direction
    ! lo(2)-1:hi(2)+1 in the y-direction
    ! 2 rows needed in the z-direction
    allocate(vlzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,2))
    allocate(vrzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,2))
    allocate(vimhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,2))

    ! lo(1):hi(1)+1 in the x-direction
    ! lo(2):hi(2) in the y-direction
    ! 2 rows needed in the z-direction
    allocate(wlxy  (lo(1):hi(1)+1,lo(2):hi(2),2))
    allocate(wrxy  (lo(1):hi(1)+1,lo(2):hi(2),2))
    allocate(wimhxy(lo(1):hi(1)+1,lo(2):hi(2),2))

    ! lo(1):hi(1) in the x-direction
    ! lo(2):hi(2)+1 in the y-direction
    ! 2 rows needed in the z-direction
    allocate(wlyx  (lo(1):hi(1),lo(2):hi(2)+1,2))
    allocate(wryx  (lo(1):hi(1),lo(2):hi(2)+1,2))
    allocate(wimhyx(lo(1):hi(1),lo(2):hi(2)+1,2))

    !************
    ! Edge states
    !************

    ! lo(1):hi(1)+1 in x-direction
    ! lo(2):hi(2) in y-direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2)))

    ! lo(1):hi(1) in x-direction
    ! lo(2):hi(2)+1 in y-direction
    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1))

    ! lo(1):hi(1) in x-direction
    ! lo(2):hi(2) in y-direction
    allocate(wmacl(lo(1):hi(1),lo(2):hi(2)))
    allocate(wmacr(lo(1):hi(1),lo(2):hi(2)))

    abs_eps = 1.0d-8

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

    ! Compute eps, which is relative to the max velocity
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
    !     1. Compute uimhx (is  :ie+1,js-1:je+1,k) 
    !     2. Compute uimhy (is-1:ie+1,js  :je+1,k)
    !     3. Compute wimhxy(is  :ie+1,js  :je  ,k)
    !     4. Compute wimhyx(is  :ie  ,js  :je+1,k)
    !     if(k .gt. ks-1) then
    !        5. Compute uimhz (is-1:ie+1,js-1:je+1,k)
    !        6. Compute wmac  (is  :ie  ,js  :je  ,k)
    !        7. Compute vimhzx(is  :ie  ,js-1:je+1,k)
    !        8. Compute uimhzy(is-1:ie+1,js  :je  ,k)
    !     endif
    !     if(k .gt. ks) then
    !        9. Compute vimhxz(is  :ie+1,js-1:je+1,k-1)
    !        10.Compute uimhyz(is-1:ie+1,js  :je+1,k-1)
    !        11.Compute umac  (is  :ie+1,js  :je,  k-1)
    !        12.Compute vmac  (is  :ie  ,js  :je+1,k-1)
    !     endif
    !     13. Cycle indices
    !  enddo
    !
    !*************************************
    ! End pseudo code
    !*************************************

    kc = 1
    kp = 2

    do k=ks-1,ke+1

       !******************************************************************
       ! 1. Compute uimhx (is  :ie+1,js-1:je+1,k) 
       !******************************************************************

       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate all components of velocity to left face
             ulx(i,j,kc,1) = u(i-1,j,k,1) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,1)
             ulx(i,j,kc,2) = u(i-1,j,k,2) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,2)
             ulx(i,j,kc,3) = u(i-1,j,k,3) + (HALF - dt2*max(ZERO,u(i-1,j,k,1))/hx)*slopex(i-1,j,k,3)

             ! extrapolate all components of velocity to right face
             urx(i,j,kc,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,1)
             urx(i,j,kc,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,2)
             urx(i,j,kc,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,1))/hx)*slopex(i,j,k,3)

             ! add source terms
             if(use_minion) then
                ulx(i,j,kc,1) = ulx(i,j,kc,1) + dt2*force(i-1,j,k,1)
                ulx(i,j,kc,2) = ulx(i,j,kc,2) + dt2*force(i-1,j,k,2)
                ulx(i,j,kc,3) = ulx(i,j,kc,3) + dt2*force(i-1,j,k,3)
                urx(i,j,kc,1) = urx(i,j,kc,1) + dt2*force(i  ,j,k,1)
                urx(i,j,kc,2) = urx(i,j,kc,2) + dt2*force(i  ,j,k,2)
                urx(i,j,kc,3) = urx(i,j,kc,3) + dt2*force(i  ,j,k,3)
             endif

             ! impose lo side bc's
             if(i .eq. is) then
                ulx(i,j,kc,1) = merge(u(is-1,j,k,1),ulx(i,j,kc,1),phys_bc(1,1) .eq. INLET)
                urx(i,j,kc,1) = merge(u(is-1,j,k,1),urx(i,j,kc,1),phys_bc(1,1) .eq. INLET)
                ulx(i,j,kc,2) = merge(u(is-1,j,k,2),ulx(i,j,kc,2),phys_bc(1,1) .eq. INLET)
                urx(i,j,kc,2) = merge(u(is-1,j,k,2),urx(i,j,kc,2),phys_bc(1,1) .eq. INLET)
                ulx(i,j,kc,3) = merge(u(is-1,j,k,3),ulx(i,j,kc,3),phys_bc(1,1) .eq. INLET)
                urx(i,j,kc,3) = merge(u(is-1,j,k,3),urx(i,j,kc,3),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   ulx(i,j,kc,1) = ZERO
                   urx(i,j,kc,1) = ZERO
                   ulx(i,j,kc,2) = merge(ZERO,urx(i,j,kc,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,kc,2) = merge(ZERO,urx(i,j,kc,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   ulx(i,j,kc,3) = merge(ZERO,urx(i,j,kc,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,kc,3) = merge(ZERO,urx(i,j,kc,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                ulx(i,j,kc,1) = merge(u(ie+1,j,k,1),ulx(i,j,kc,1),phys_bc(1,2) .eq. INLET)
                urx(i,j,kc,1) = merge(u(ie+1,j,k,1),urx(i,j,kc,1),phys_bc(1,2) .eq. INLET)
                ulx(i,j,kc,2) = merge(u(ie+1,j,k,2),ulx(i,j,kc,2),phys_bc(1,2) .eq. INLET)
                urx(i,j,kc,2) = merge(u(ie+1,j,k,2),urx(i,j,kc,2),phys_bc(1,2) .eq. INLET)
                ulx(i,j,kc,3) = merge(u(ie+1,j,k,3),ulx(i,j,kc,3),phys_bc(1,2) .eq. INLET)
                urx(i,j,kc,3) = merge(u(ie+1,j,k,3),urx(i,j,kc,3),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   ulx(i,j,kc,1) = ZERO
                   urx(i,j,kc,1) = ZERO
                   ulx(i,j,kc,2) = merge(ZERO,ulx(i,j,kc,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,kc,2) = merge(ZERO,ulx(i,j,kc,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   ulx(i,j,kc,3) = merge(ZERO,ulx(i,j,kc,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,kc,3) = merge(ZERO,ulx(i,j,kc,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhx by first solving a normal Riemann problem
             uavg = HALF*(ulx(i,j,kc,1)+urx(i,j,kc,1))
             test = ((ulx(i,j,kc,1) .le. ZERO .and. urx(i,j,kc,1) .ge. ZERO) .or. &
                  (abs(ulx(i,j,kc,1)+urx(i,j,kc,1)) .lt. eps))
             uimhx(i,j,kc,1) = merge(ulx(i,j,kc,1),urx(i,j,kc,1),uavg .gt. ZERO)
             uimhx(i,j,kc,1) = merge(ZERO,uimhx(i,j,kc,1),test)

             ! now upwind to get transverse components of uimhx
             uimhx(i,j,kc,2) = merge(ulx(i,j,kc,2),urx(i,j,kc,2),uimhx(i,j,kc,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,kc,2)+urx(i,j,kc,2))
             uimhx(i,j,kc,2) = merge(uavg,uimhx(i,j,kc,2),abs(uimhx(i,j,kc,1)).lt.eps)

             uimhx(i,j,kc,3) = merge(ulx(i,j,kc,3),urx(i,j,kc,3),uimhx(i,j,kc,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,kc,3)+urx(i,j,kc,3))
             uimhx(i,j,kc,3) = merge(uavg,uimhx(i,j,kc,3),abs(uimhx(i,j,kc,1)).lt.eps)
          enddo
       enddo

       !******************************************************************
       ! 2. Compute uimhy (is-1:ie+1,js  :je+1,k)
       !******************************************************************

       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate all components of velocity to left face
             uly(i,j,kc,1) = u(i,j-1,k,1) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,1)
             uly(i,j,kc,2) = u(i,j-1,k,2) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,2)
             uly(i,j,kc,3) = u(i,j-1,k,3) + (HALF - dt2*max(ZERO,u(i,j-1,k,2)/hy))*slopey(i,j-1,k,3)

             ! extrapolate all components of velocity to right face
             ury(i,j,kc,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,1)
             ury(i,j,kc,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,2)
             ury(i,j,kc,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,2))/hy)*slopey(i,j,k,3)

             ! add source terms
             if(use_minion) then
                uly(i,j,kc,1) = uly(i,j,kc,1) + dt2*force(i,j-1,k,1)
                uly(i,j,kc,2) = uly(i,j,kc,2) + dt2*force(i,j-1,k,2)
                uly(i,j,kc,3) = uly(i,j,kc,3) + dt2*force(i,j-1,k,3)
                ury(i,j,kc,1) = ury(i,j,kc,1) + dt2*force(i,j  ,k,1)
                ury(i,j,kc,2) = ury(i,j,kc,2) + dt2*force(i,j  ,k,2)
                ury(i,j,kc,3) = ury(i,j,kc,3) + dt2*force(i,j  ,k,3)
             endif

             ! impose lo side bc's
             if(j .eq. js) then
                uly(i,j,kc,1) = merge(u(i,js-1,k,1),uly(i,j,kc,1),phys_bc(2,1) .eq. INLET)
                ury(i,j,kc,1) = merge(u(i,js-1,k,1),ury(i,j,kc,1),phys_bc(2,1) .eq. INLET)
                uly(i,j,kc,2) = merge(u(i,js-1,k,2),uly(i,j,kc,2),phys_bc(2,1) .eq. INLET)
                ury(i,j,kc,2) = merge(u(i,js-1,k,2),ury(i,j,kc,2),phys_bc(2,1) .eq. INLET)
                uly(i,j,kc,3) = merge(u(i,js-1,k,3),uly(i,j,kc,3),phys_bc(2,1) .eq. INLET)
                ury(i,j,kc,3) = merge(u(i,js-1,k,3),ury(i,j,kc,3),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   uly(i,j,kc,1) = merge(ZERO,ury(i,j,kc,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,kc,1) = merge(ZERO,ury(i,j,kc,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uly(i,j,kc,2) = ZERO
                   ury(i,j,kc,2) = ZERO
                   uly(i,j,kc,3) = merge(ZERO,ury(i,j,kc,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,kc,3) = merge(ZERO,ury(i,j,kc,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                uly(i,j,kc,1) = merge(u(i,je+1,k,1),uly(i,j,kc,1),phys_bc(2,2) .eq. INLET)
                ury(i,j,kc,1) = merge(u(i,je+1,k,1),ury(i,j,kc,1),phys_bc(2,2) .eq. INLET)
                uly(i,j,kc,2) = merge(u(i,je+1,k,2),uly(i,j,kc,2),phys_bc(2,2) .eq. INLET)
                ury(i,j,kc,2) = merge(u(i,je+1,k,2),ury(i,j,kc,2),phys_bc(2,2) .eq. INLET)
                uly(i,j,kc,3) = merge(u(i,je+1,k,3),uly(i,j,kc,3),phys_bc(2,2) .eq. INLET)
                ury(i,j,kc,3) = merge(u(i,je+1,k,3),ury(i,j,kc,3),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   uly(i,j,kc,1) = merge(ZERO,uly(i,j,kc,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,kc,1) = merge(ZERO,uly(i,j,kc,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uly(i,j,kc,2) = ZERO
                   ury(i,j,kc,2) = ZERO
                   uly(i,j,kc,3) = merge(ZERO,uly(i,j,kc,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,kc,3) = merge(ZERO,uly(i,j,kc,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhy by first solving a normal Riemann problem
             uavg = HALF*(uly(i,j,kc,2)+ury(i,j,kc,2))
             test = ((uly(i,j,kc,2) .le. ZERO .and. ury(i,j,kc,2) .ge. ZERO) .or. &
                  (abs(uly(i,j,kc,2)+ury(i,j,kc,2)) .lt. eps))
             uimhy(i,j,kc,2) = merge(uly(i,j,kc,2),ury(i,j,kc,2),uavg .gt. ZERO)
             uimhy(i,j,kc,2) = merge(ZERO,uimhy(i,j,kc,2),test)

             ! now upwind to get transverse components of uimhy
             uimhy(i,j,kc,1) = merge(uly(i,j,kc,1),ury(i,j,kc,1),uimhy(i,j,kc,2).gt.ZERO)
             uavg = HALF*(uly(i,j,kc,1)+ury(i,j,kc,1))
             uimhy(i,j,kc,1) = merge(uavg,uimhy(i,j,kc,1),abs(uimhy(i,j,kc,2)).lt.eps)

             uimhy(i,j,kc,3) = merge(uly(i,j,kc,3),ury(i,j,kc,3),uimhy(i,j,kc,2).gt.ZERO)
             uavg = HALF*(uly(i,j,kc,3)+ury(i,j,kc,3))
             uimhy(i,j,kc,3) = merge(uavg,uimhy(i,j,kc,3),abs(uimhy(i,j,kc,2)).lt.eps)
          enddo
       enddo

       !******************************************************************
       ! 3. Compute wimhxy(is  :ie+1,js  :je  ,k)
       !******************************************************************

       do j=js,je
          do i=is,ie+1
             ! extrapolate to faces
             wlxy(i,j,kc) = ulx(i,j,kc,3) &
                  - (dt6/hy)*(uimhy(i-1,j+1,kc,2)+uimhy(i-1,j,kc,2))*(uimhy(i-1,j+1,kc,3)-uimhy(i-1,j,kc,3))
             wrxy(i,j,kc) = urx(i,j,kc,3) &
                  - (dt6/hy)*(uimhy(i  ,j+1,kc,2)+uimhy(i  ,j,kc,2))*(uimhy(i  ,j+1,kc,3)-uimhy(i  ,j,kc,3))

             ! impose lo side bc's
             if(i .eq. is) then
                wlxy(i,j,kc) = merge(u(is-1,j,k,3),wlxy(i,j,kc),phys_bc(1,1) .eq. INLET)
                wrxy(i,j,kc) = merge(u(is-1,j,k,3),wrxy(i,j,kc),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,kc) = merge(ZERO,wrxy(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   wrxy(i,j,kc) = merge(ZERO,wrxy(i,j,kc),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                wlxy(i,j,kc) = merge(u(ie+1,j,k,3),wlxy(i,j,kc),phys_bc(1,2) .eq. INLET)
                wrxy(i,j,kc) = merge(u(ie+1,j,k,3),wrxy(i,j,kc),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,kc) = merge(ZERO,wlxy(i,j,kc),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   wrxy(i,j,kc) = merge(ZERO,wlxy(i,j,kc),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhxy(i,j,kc) = merge(wlxy(i,j,kc),wrxy(i,j,kc),uimhx(i,j,kc,1).gt.ZERO)
             uavg = HALF*(wlxy(i,j,kc)+wrxy(i,j,kc))
             wimhxy(i,j,kc) = merge(uavg,wimhxy(i,j,kc),abs(uimhx(i,j,kc,1)).lt.eps)
          enddo
       enddo

       !******************************************************************
       ! 4. Compute wimhyx(is  :ie  ,js  :je+1,k)
       !******************************************************************

       do j=js,je+1
          do i=is,ie
             ! extrapolate to faces
             wlyx(i,j,kc) = uly(i,j,kc,3) &
                  - (dt6/hx)*(uimhx(i+1,j-1,kc,1)+uimhx(i,j-1,kc,1))*(uimhx(i+1,j-1,kc,3)-uimhx(i,j-1,kc,3))
             wryx(i,j,kc) = ury(i,j,kc,3) &
                  - (dt6/hx)*(uimhx(i+1,j  ,kc,1)+uimhx(i,j  ,kc,1))*(uimhx(i+1,j  ,kc,3)-uimhx(i,j  ,kc,3))

             ! impose lo side bc's
             if(j .eq. js) then
                wlyx(i,j,kc) = merge(u(i,js-1,k,3),wlyx(i,j,kc),phys_bc(2,1) .eq. INLET)
                wryx(i,j,kc) = merge(u(i,js-1,k,3),wryx(i,j,kc),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,kc) = merge(ZERO,wryx(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   wryx(i,j,kc) = merge(ZERO,wryx(i,j,kc),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                wlyx(i,j,kc) = merge(u(i,je+1,k,3),wlyx(i,j,kc),phys_bc(2,2) .eq. INLET)
                wryx(i,j,kc) = merge(u(i,je+1,k,3),wryx(i,j,kc),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,kc) = merge(ZERO,wlyx(i,j,kc),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   wryx(i,j,kc) = merge(ZERO,wlyx(i,j,kc),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhyx(i,j,kc) = merge(wlyx(i,j,kc),wryx(i,j,kc),uimhy(i,j,kc,2).gt.ZERO)
             uavg = HALF*(wlyx(i,j,kc)+wryx(i,j,kc))
             wimhyx(i,j,kc) = merge(uavg,wimhyx(i,j,kc),abs(uimhy(i,j,kc,2)).lt.eps)
          enddo
       enddo

       if(k .gt. ks-1) then

          !******************************************************************
          ! 5. Compute uimhz (is-1:ie+1,js-1:je+1,k)
          !******************************************************************

          do j=js-1,je+1
             do i=is-1,ie+1
                ! extrapolate all components of velocity to left face
                ulz(i,j,kc,1) = u(i,j,k-1,1) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,1)
                ulz(i,j,kc,2) = u(i,j,k-1,2) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,2)
                ulz(i,j,kc,3) = u(i,j,k-1,3) + (HALF - dt2*max(ZERO,u(i,j,k-1,3))/hz)*slopez(i,j,k-1,3)

                ! extrapolate all components of velocity to right face
                urz(i,j,kc,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,1)
                urz(i,j,kc,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,2)
                urz(i,j,kc,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,u(i,j,k,3))/hz)*slopez(i,j,k,3)

                ! add source terms
                if(use_minion) then
                   ulz(i,j,kc,1) = ulz(i,j,kc,1) + dt2*force(i,j,k-1,1)
                   ulz(i,j,kc,2) = ulz(i,j,kc,2) + dt2*force(i,j,k-1,2)
                   ulz(i,j,kc,3) = ulz(i,j,kc,3) + dt2*force(i,j,k-1,3)
                   urz(i,j,kc,1) = urz(i,j,kc,1) + dt2*force(i,j,k  ,1)
                   urz(i,j,kc,2) = urz(i,j,kc,2) + dt2*force(i,j,k  ,2)
                   urz(i,j,kc,3) = urz(i,j,kc,3) + dt2*force(i,j,k  ,3)
                endif

                ! impose lo side bc's
                if(k .eq. ks) then
                   ulz(i,j,kc,1) = merge(u(i,j,ks-1,1),ulz(i,j,kc,1),phys_bc(3,1) .eq. INLET)
                   urz(i,j,kc,1) = merge(u(i,j,ks-1,1),urz(i,j,kc,1),phys_bc(3,1) .eq. INLET)
                   ulz(i,j,kc,2) = merge(u(i,j,ks-1,2),ulz(i,j,kc,2),phys_bc(3,1) .eq. INLET)
                   urz(i,j,kc,2) = merge(u(i,j,ks-1,2),urz(i,j,kc,2),phys_bc(3,1) .eq. INLET)
                   ulz(i,j,kc,3) = merge(u(i,j,ks-1,3),ulz(i,j,kc,3),phys_bc(3,1) .eq. INLET)
                   urz(i,j,kc,3) = merge(u(i,j,ks-1,3),urz(i,j,kc,3),phys_bc(3,1) .eq. INLET)
                   if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                      ulz(i,j,kc,1) = merge(ZERO,urz(i,j,kc,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      urz(i,j,kc,1) = merge(ZERO,urz(i,j,kc,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      ulz(i,j,kc,2) = merge(ZERO,urz(i,j,kc,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      urz(i,j,kc,2) = merge(ZERO,urz(i,j,kc,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      ulz(i,j,kc,3) = ZERO
                      urz(i,j,kc,3) = ZERO
                   endif
                endif

                ! impose hi side bc's
                if(k .eq. ke+1) then
                   ulz(i,j,kc,1) = merge(u(i,j,ke+1,1),ulz(i,j,kc,1),phys_bc(3,2) .eq. INLET)
                   urz(i,j,kc,1) = merge(u(i,j,ke+1,1),urz(i,j,kc,1),phys_bc(3,2) .eq. INLET)
                   ulz(i,j,kc,2) = merge(u(i,j,ke+1,2),ulz(i,j,kc,2),phys_bc(3,2) .eq. INLET)
                   urz(i,j,kc,2) = merge(u(i,j,ke+1,2),urz(i,j,kc,2),phys_bc(3,2) .eq. INLET)
                   ulz(i,j,kc,3) = merge(u(i,j,ke+1,3),ulz(i,j,kc,3),phys_bc(3,2) .eq. INLET)
                   urz(i,j,kc,3) = merge(u(i,j,ke+1,3),urz(i,j,kc,3),phys_bc(3,2) .eq. INLET)
                   if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                      ulz(i,j,kc,1) = merge(ZERO,ulz(i,j,kc,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      urz(i,j,kc,1) = merge(ZERO,ulz(i,j,kc,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      ulz(i,j,kc,2) = merge(ZERO,ulz(i,j,kc,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      urz(i,j,kc,2) = merge(ZERO,ulz(i,j,kc,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      ulz(i,j,kc,3) = ZERO
                      urz(i,j,kc,3) = ZERO
                   endif
                endif

                ! make normal component of uimhz by first solving a normal Riemann problem
                uavg = HALF*(ulz(i,j,kc,3)+urz(i,j,kc,3))
                test = ((ulz(i,j,kc,3) .le. ZERO .and. urz(i,j,kc,3) .ge. ZERO) .or. &
                     (abs(ulz(i,j,kc,3)+urz(i,j,kc,3)) .lt. eps))
                uimhz(i,j,kc,3) = merge(ulz(i,j,kc,3),urz(i,j,kc,3),uavg .gt. ZERO)
                uimhz(i,j,kc,3) = merge(ZERO,uimhz(i,j,kc,3),test)

                ! now upwind to get transverse components of uimhz
                uimhz(i,j,kc,1) = merge(ulz(i,j,kc,1),urz(i,j,kc,1),uimhz(i,j,kc,3).gt.ZERO)
                uavg = HALF*(ulz(i,j,kc,1)+urz(i,j,kc,1))
                uimhz(i,j,kc,1) = merge(uavg,uimhz(i,j,kc,1),abs(uimhz(i,j,kc,3)).lt.eps)

                uimhz(i,j,kc,2) = merge(ulz(i,j,kc,2),urz(i,j,kc,2),uimhz(i,j,kc,3).gt.ZERO)
                uavg = HALF*(ulz(i,j,kc,2)+urz(i,j,kc,2))
                uimhz(i,j,kc,2) = merge(uavg,uimhz(i,j,kc,2),abs(uimhz(i,j,kc,3)).lt.eps)
             enddo
          enddo

          !******************************************************************
          ! 6. Compute wmac  (is  :ie  ,js  :je  ,k)
          !******************************************************************

          do j=js,je
             do i=is,ie
                ! extrapolate to edges
                wmacl(i,j) = ulz(i,j,kc,3) &
                     - (dt4/hx)*(uimhx(i+1,j,kp,1)+uimhx(i,j,kp,1))*(wimhxy(i+1,j,kp)-wimhxy(i,j,kp)) &
                     - (dt4/hy)*(uimhy(i,j+1,kp,2)+uimhy(i,j,kp,2))*(wimhyx(i,j+1,kp)-wimhyx(i,j,kp))
                wmacr(i,j) = urz(i,j,kc,3) &
                     - (dt4/hx)*(uimhx(i+1,j,kc,1)+uimhx(i,j,kc,1))*(wimhxy(i+1,j,kc)-wimhxy(i,j,kc)) &
                     - (dt4/hy)*(uimhy(i,j+1,kc,2)+uimhy(i,j,kc,2))*(wimhyx(i,j+1,kc)-wimhyx(i,j,kc))

                ! if use_minion is true, we have already accounted for source terms
                ! in uly and ury; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   wmacl(i,j) = wmacl(i,j) + dt2*force(i,j,k-1,3)
                   wmacr(i,j) = wmacr(i,j) + dt2*force(i,j,k  ,3)
                endif

                ! solve Riemann problem
                uavg = HALF*(wmacl(i,j)+wmacr(i,j))
                test = ((wmacl(i,j) .le. ZERO .and. wmacr(i,j) .ge. ZERO) .or. &
                     (abs(wmacl(i,j)+wmacr(i,j)) .lt. eps))
                wmac(i,j,k) = merge(wmacl(i,j),wmacr(i,j),uavg .gt. ZERO)
                wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)

                ! Apply boundary conditions
                if(k .eq. ks) then
                   ! lo side
                   if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                      wmac(i,j,ks) = ZERO
                   elseif (phys_bc(3,1) .eq. INLET) then
                      wmac(i,j,ks) = u(i,j,ks-1,3)
                   elseif (phys_bc(3,1) .eq. OUTLET) then
                      wmac(i,j,ks) = min(wmacr(i,j),ZERO)
                   endif
                else if(k .eq. ke+1) then
                   ! hi side
                   if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                      wmac(i,j,ke+1) = ZERO
                   elseif (phys_bc(3,2) .eq. INLET) then
                      wmac(i,j,ke+1) = u(i,j,ke+1,3)
                   elseif (phys_bc(3,2) .eq. OUTLET) then
                      wmac(i,j,ke+1) = max(wmacl(i,j),ZERO)
                   endif
                endif
             enddo
          enddo

          !******************************************************************
          ! 7. Compute vimhzx(is  :ie  ,js-1:je+1,k)
          !******************************************************************

          do j=js-1,je+1
             do i=is,ie
                ! extrapolate to faces
                vlzx(i,j,kc) = ulz(i,j,kc,2) &
                     - (dt6/hx)*(uimhx(i+1,j,kp,1)+uimhx(i,j,kp,1))*(uimhx(i+1,j,kp,2)-uimhx(i,j,kp,2))
                vrzx(i,j,kc) = urz(i,j,kc,2) &
                     - (dt6/hx)*(uimhx(i+1,j,kc,1)+uimhx(i,j,kc,1))*(uimhx(i+1,j,kc,2)-uimhx(i,j,kc,2))

                ! impose lo side bc's
                if(k .eq. ks) then
                   vlzx(i,j,kc) = merge(u(i,j,ks-1,1),vlzx(i,j,kc),phys_bc(3,1) .eq. INLET)
                   vrzx(i,j,kc) = merge(u(i,j,ks-1,1),vrzx(i,j,kc),phys_bc(3,1) .eq. INLET)
                   if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                      vlzx(i,j,kc) = merge(ZERO,vrzx(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      vrzx(i,j,kc) = merge(ZERO,vrzx(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! impose hi side bc's
                if(k .eq. ke+1) then
                   vlzx(i,j,kc) = merge(u(i,j,ke+1,1),vlzx(i,j,kc),phys_bc(3,2) .eq. INLET)
                   vrzx(i,j,kc) = merge(u(i,j,ke+1,1),vrzx(i,j,kc),phys_bc(3,2) .eq. INLET)
                   if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                      vlzx(i,j,kc) = merge(ZERO,vlzx(i,j,kc),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      vrzx(i,j,kc) = merge(ZERO,vlzx(i,j,kc),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! upwind
                vimhzx(i,j,kc) = merge(vlzx(i,j,kc),vrzx(i,j,kc),uimhz(i,j,kc,3).gt.ZERO)
                uavg = HALF*(vlzx(i,j,kc)+vrzx(i,j,kc))
                vimhzx(i,j,kc) = merge(uavg,vimhzx(i,j,kc),abs(uimhz(i,j,kc,3)).lt.eps)
             enddo
          enddo

          !******************************************************************
          ! 8. Compute uimhzy(is-1:ie+1,js  :je  ,k)
          !******************************************************************

          do j=js,je
             do i=is-1,ie+1
                ! extrapolate to faces
                ulzy(i,j,kc) = ulz(i,j,kc,1) &
                     - (dt6/hy)*(uimhy(i,j+1,kp,2)+uimhy(i,j,kp,2))*(uimhy(i,j+1,kp,1)-uimhy(i,j,kp,1))
                urzy(i,j,kc) = urz(i,j,kc,1) &
                     - (dt6/hy)*(uimhy(i,j+1,kc,2)+uimhy(i,j,kc,2))*(uimhy(i,j+1,kc,1)-uimhy(i,j,kc,1))

                ! impose lo side bc's
                if(k .eq. ks) then
                   ulzy(i,j,kc) = merge(u(i,j,ks-1,1),ulzy(i,j,kc),phys_bc(3,1) .eq. INLET)
                   urzy(i,j,kc) = merge(u(i,j,ks-1,1),urzy(i,j,kc),phys_bc(3,1) .eq. INLET)
                   if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                      ulzy(i,j,kc) = merge(ZERO,urzy(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      urzy(i,j,kc) = merge(ZERO,urzy(i,j,kc),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! impose hi side bc's
                if(k .eq. ke+1) then
                   ulzy(i,j,kc) = merge(u(i,j,ke+1,1),ulzy(i,j,kc),phys_bc(3,2) .eq. INLET)
                   urzy(i,j,kc) = merge(u(i,j,ke+1,1),urzy(i,j,kc),phys_bc(3,2) .eq. INLET)
                   if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                      ulzy(i,j,kc) = merge(ZERO,ulzy(i,j,kc),phys_bc(3,2) .eq. NO_SLIP_WALL)
                      urzy(i,j,kc) = merge(ZERO,ulzy(i,j,kc),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! upwind
                uimhzy(i,j,kc) = merge(ulzy(i,j,kc),urzy(i,j,kc),uimhz(i,j,kc,3).gt.ZERO)
                uavg = HALF*(ulzy(i,j,kc)+urzy(i,j,kc))
                uimhzy(i,j,kc) = merge(uavg,uimhzy(i,j,kc),abs(uimhz(i,j,kc,3)).lt.eps)
             enddo
          enddo

       endif ! end if(k .gt. ks-1)

       if(k .gt. ks) then

          !******************************************************************
          ! 9. Compute vimhxz(is  :ie+1,js-1:je+1,k-1)
          !******************************************************************

          do j=js-1,je+1
             do i=is,ie+1
                ! extrapolate to faces
                vlxz(i,j,kp) = ulx(i,j,kp,2) &
                     - (dt6/hz)*(uimhz(i-1,j,kc,3)+uimhz(i-1,j,kp,3))*(uimhz(i-1,j,kc,2)-uimhz(i-1,j,kp,2))
                vrxz(i,j,kp) = urx(i,j,kp,2) &
                     - (dt6/hz)*(uimhz(  i,j,kc,3)+uimhz(i  ,j,kp,3))*(uimhz(i  ,j,kc,2)-uimhz(i  ,j,kp,2))

                ! impose lo side bc's
                if(i .eq. is) then
                   vlxz(i,j,kp) = merge(u(is-1,j,k-1,2),vlxz(i,j,kp),phys_bc(1,1) .eq. INLET)
                   vrxz(i,j,kp) = merge(u(is-1,j,k-1,2),vrxz(i,j,kp),phys_bc(1,1) .eq. INLET)
                   if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                      vlxz(i,j,kp) = merge(ZERO,vrxz(i,j,kp),phys_bc(1,1) .eq. NO_SLIP_WALL)
                      vrxz(i,j,kp) = merge(ZERO,vrxz(i,j,kp),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! impose hi side bc's
                if(i .eq. ie+1) then
                   vlxz(i,j,kp) = merge(u(ie+1,j,k-1,2),vlxz(i,j,kp),phys_bc(1,2) .eq. INLET)
                   vrxz(i,j,kp) = merge(u(ie+1,j,k-1,2),vrxz(i,j,kp),phys_bc(1,2) .eq. INLET)
                   if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                      vlxz(i,j,kp) = merge(ZERO,vlxz(i,j,kp),phys_bc(1,2) .eq. NO_SLIP_WALL)
                      vrxz(i,j,kp) = merge(ZERO,vlxz(i,j,kp),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! upwind
                vimhxz(i,j,kp) = merge(vlxz(i,j,kp),vrxz(i,j,kp),uimhx(i,j,kp,1).gt.ZERO)
                uavg = HALF*(vlxz(i,j,kp)+vrxz(i,j,kp))
                vimhxz(i,j,kp) = merge(uavg,vimhxz(i,j,kp),abs(uimhx(i,j,kp,1)).lt.eps)
             enddo
          enddo

          !******************************************************************
          ! 10.Compute uimhyz(is-1:ie+1,js  :je+1,k-1)
          !******************************************************************

          do j=js,je+1
             do i=is-1,ie+1
                ! extrapolate to faces
                ulyz(i,j,kp) = uly(i,j,kp,1) &
                     - (dt6/hz)*(uimhz(i,j-1,kc,3)+uimhz(i,j-1,kp,3))*(uimhz(i,j-1,kc,1)-uimhz(i,j-1,kp,1))               
                uryz(i,j,kp) = ury(i,j,kp,1) &
                     - (dt6/hz)*(uimhz(i,j  ,kc,3)+uimhz(i,j  ,kp,3))*(uimhz(i,j  ,kc,1)-uimhz(i,j  ,kp,1)) 

                ! impose lo side bc's
                if(j .eq. js) then
                   ulyz(i,j,kp) = merge(u(i,js-1,k-1,1),ulyz(i,j,kp),phys_bc(2,1) .eq. INLET)
                   uryz(i,j,kp) = merge(u(i,js-1,k-1,1),uryz(i,j,kp),phys_bc(2,1) .eq. INLET)
                   if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                      ulyz(i,j,kp) = merge(ZERO,uryz(i,j,kp),phys_bc(2,1) .eq. NO_SLIP_WALL)
                      uryz(i,j,kp) = merge(ZERO,uryz(i,j,kp),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! impose hi side bc's
                if(j .eq. je+1) then
                   ulyz(i,j,kp) = merge(u(i,je+1,k-1,1),ulyz(i,j,kp),phys_bc(2,2) .eq. INLET)
                   uryz(i,j,kp) = merge(u(i,je+1,k-1,1),uryz(i,j,kp),phys_bc(2,2) .eq. INLET)
                   if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                      ulyz(i,j,kp) = merge(ZERO,ulyz(i,j,kp),phys_bc(2,2) .eq. NO_SLIP_WALL)
                      uryz(i,j,kp) = merge(ZERO,ulyz(i,j,kp),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   endif
                endif

                ! upwind
                uimhyz(i,j,kp) = merge(ulyz(i,j,kp),uryz(i,j,kp),uimhy(i,j,kp,2).gt.ZERO)
                uavg = HALF*(ulyz(i,j,kp)+uryz(i,j,kp))
                uimhyz(i,j,kp) = merge(uavg,uimhyz(i,j,kp),abs(uimhy(i,j,kp,2)).lt.eps)
             enddo
          enddo

          !******************************************************************
          ! 11.Compute umac  (is  :ie+1,js  :je,  k-1)
          !******************************************************************

          do j=js,je
             do i=is,ie+1
                ! extrapolate to edges
                umacl(i,j) = ulx(i,j,kp,1) &
                     - (dt4/hy)*(uimhy(i-1,j+1,kp,2)+uimhy(i-1,j,kp,2))*(uimhyz(i-1,j+1,kp)-uimhyz(i-1,j,kp)) &
                     - (dt4/hz)*(uimhz(i-1,j  ,kc,3)+uimhz(i-1,j,kp,3))*(uimhzy(i-1,j  ,kc)-uimhzy(i-1,j,kp))
                umacr(i,j) = urx(i,j,kp,1) &
                     - (dt4/hy)*(uimhy(i  ,j+1,kp,2)+uimhy(i  ,j,kp,2))*(uimhyz(i  ,j+1,kp)-uimhyz(i  ,j,kp)) &
                     - (dt4/hz)*(uimhz(i  ,j  ,kc,3)+uimhz(i  ,j,kp,3))*(uimhzy(i  ,j  ,kc)-uimhzy(i  ,j,kp))

                ! if use_minion is true, we have already accounted for source terms
                ! in ulx and urx; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   umacl(i,j) = umacl(i,j) + dt2*force(i-1,j,k-1,1)
                   umacr(i,j) = umacr(i,j) + dt2*force(i  ,j,k-1,1)
                endif

                ! solve Riemann problem
                uavg = HALF*(umacl(i,j)+umacr(i,j))
                test = ((umacl(i,j) .le. ZERO .and. umacr(i,j) .ge. ZERO) .or. &
                     (abs(umacl(i,j)+umacr(i,j)) .lt. eps))
                umac(i,j,k-1) = merge(umacl(i,j),umacr(i,j),uavg .gt. ZERO)
                umac(i,j,k-1) = merge(ZERO,umac(i,j,k-1),test)

                if(i .eq. is) then
                   ! lo side
                   if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                      umac(is,j,k-1) = ZERO
                   elseif (phys_bc(1,1) .eq. INLET) then
                      umac(is,j,k-1) = u(is-1,j,k-1,1)
                   elseif (phys_bc(1,1) .eq. OUTLET) then
                      umac(is,j,k-1) = min(umacr(is,j),ZERO)
                   endif
                else if(i .eq. ie+1) then
                   ! hi side
                   if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                      umac(ie+1,j,k-1) = ZERO
                   elseif (phys_bc(1,2) .eq. INLET) then
                      umac(ie+1,j,k-1) = u(ie+1,j,k-1,1)
                   elseif (phys_bc(1,2) .eq. OUTLET) then
                      umac(ie+1,j,k-1) = max(umacl(ie+1,j),ZERO)
                   endif
                endif
             enddo
          enddo

          !******************************************************************
          ! 12.Compute vmac  (is  :ie  ,js  :je+1,k-1)
          !******************************************************************

          do j=js,je+1
             do i=is,ie
                ! extrapolate to edges
                vmacl(i,j) = uly(i,j,kp,2) &
                     - (dt4/hx)*(uimhx(i+1,j-1,kp,1)+uimhx(i,j-1,kp,1))*(vimhxz(i+1,j-1,kp)-vimhxz(i,j-1,kp)) &
                     - (dt4/hz)*(uimhz(i  ,j-1,kc,3)+uimhz(i,j-1,kp,3))*(vimhzx(i  ,j-1,kc)-vimhzx(i,j-1,kp))
                vmacr(i,j) = ury(i,j,kp,2) &
                     - (dt4/hx)*(uimhx(i+1,j  ,kp,1)+uimhx(i,j  ,kp,1))*(vimhxz(i+1,j  ,kp)-vimhxz(i,j  ,kp)) &
                     - (dt4/hz)*(uimhz(i  ,j  ,kc,3)+uimhz(i,j  ,kp,3))*(vimhzx(i  ,j  ,kc)-vimhzx(i,j  ,kp))

                ! if use_minion is true, we have already accounted for source terms
                ! in uly and ury; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   vmacl(i,j) = vmacl(i,j) + dt2*force(i,j-1,k-1,2)
                   vmacr(i,j) = vmacr(i,j) + dt2*force(i,j  ,k-1,2)
                endif

                ! solve Riemann problem
                uavg = HALF*(vmacl(i,j)+vmacr(i,j))
                test = ((vmacl(i,j) .le. ZERO .and. vmacr(i,j) .ge. ZERO) .or. &
                     (abs(vmacl(i,j)+vmacr(i,j)) .lt. eps))
                vmac(i,j,k-1) = merge(vmacl(i,j),vmacr(i,j),uavg .gt. ZERO)
                vmac(i,j,k-1) = merge(ZERO,vmac(i,j,k-1),test)

                ! Apply boundary conditions
                if(j .eq. js) then
                   ! lo side
                   if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                      vmac(i,js,k-1) = ZERO
                   elseif (phys_bc(2,1) .eq. INLET) then
                      vmac(i,js,k-1) = u(i,js-1,k-1,2)
                   elseif (phys_bc(2,1) .eq. OUTLET) then
                      vmac(i,js,k-1) = min(vmacr(i,js),ZERO)
                   endif
                else if(j .eq. je+1) then
                   ! hi side
                   if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                      vmac(i,je+1,k-1) = ZERO
                   elseif (phys_bc(2,2) .eq. INLET) then
                      vmac(i,je+1,k-1) = u(i,je+1,k-1,2)
                   elseif (phys_bc(2,2) .eq. OUTLET) then
                      vmac(i,je+1,k-1) = max(vmacl(i,je+1),ZERO)
                   endif
                endif
             enddo
          enddo

       endif ! end if(k .gt. ks)

       !******************************************************************
       ! 13. Cycle indices
       !******************************************************************

       kc = 3 - kc
       kp = 3 - kp

    enddo

    deallocate(slopex)
    deallocate(slopey)
    deallocate(slopez)

    deallocate(ulx)
    deallocate(urx)
    deallocate(uimhx)
    deallocate(uly)
    deallocate(ury)
    deallocate(uimhy)
    deallocate(ulz)
    deallocate(urz)
    deallocate(uimhz)

    deallocate(ulyz)
    deallocate(uryz)
    deallocate(uimhyz)

    deallocate(ulzy)
    deallocate(urzy)
    deallocate(uimhzy)

    deallocate(vlxz)
    deallocate(vrxz)
    deallocate(vimhxz)

    deallocate(vlzx)
    deallocate(vrzx)
    deallocate(vimhzx)

    deallocate(wlxy)
    deallocate(wrxy)
    deallocate(wimhxy)

    deallocate(wlyx)
    deallocate(wryx)
    deallocate(wimhyx)

    deallocate(umacl)
    deallocate(umacr)
    deallocate(vmacl)
    deallocate(vmacr)
    deallocate(wmacl)
    deallocate(wmacr)

  end subroutine velpred_3d

end module velpred_module
