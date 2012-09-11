module mkflux_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: mkflux

contains

  subroutine mkflux(mla,sold,sedge,flux,umac,force,divu,dx,dt,the_bc_level, &
                    is_vel,is_conservative)

    use ml_restriction_module, only: ml_edge_restriction_c
    use probin_module, only: use_godunov_debug

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    type(multifab) , intent(in   ) :: divu(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical        , intent(in   ) :: is_vel,is_conservative(:)

    ! local
    integer                  :: n,i,dm,ng,comp,ncomp,bccomp,nlevs,gid
    integer                  :: lo(get_dim(sold(1))),hi(get_dim(sold(1)))
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: dp(:,:,:,:)

    nlevs = mla%nlevel
    dm    = mla%dim

    ng = nghost(sold(1))
    ncomp = multifab_ncomp(sold(1))

    if(is_vel) then
       bccomp = 1
    else
       bccomp = dm+1
    endif

    do n=1,nlevs
       do i = 1, nfabs(sold(n))
          gid = global_index(sold(n),i)
          sop    => dataptr(sold(n), i)
          sepx   => dataptr(sedge(n,1), i)
          sepy   => dataptr(sedge(n,2), i)
          fluxpx => dataptr(flux(n,1), i)
          fluxpy => dataptr(flux(n,2), i)
          ump    => dataptr(umac(n,1), i)
          vmp    => dataptr(umac(n,2), i)
          fp     => dataptr(force(n) , i)
          dp     => dataptr(divu(n), i)
          lo = lwb(get_box(sold(n), i))
          hi = upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             if(use_godunov_debug) then
                call mkflux_debug_2d(sop(:,:,1,:),  &
                                     sepx(:,:,1,:), sepy(:,:,1,:), &
                                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                     ump(:,:,1,1), vmp(:,:,1,1), &
                                     fp(:,:,1,:), dp(:,:,1,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(gid,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(gid,:,:,bccomp:bccomp+ncomp-1),&
                                     ng, is_conservative)
             else
                call mkflux_2d(sop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               fp(:,:,1,:), dp(:,:,1,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(gid,:,:), &
                               the_bc_level(n)%adv_bc_level_array(gid,:,:,bccomp:bccomp+ncomp-1),&
                               ng, is_conservative)
             endif
          case (3)
             sepz   => dataptr(sedge(n,3), i)
             fluxpz => dataptr(flux(n,3), i)
             wmp    => dataptr(umac(n,3), i)
             if(use_godunov_debug) then
                call mkflux_debug_3d(sop(:,:,:,:), &
                                     sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                     fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                     ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                     fp(:,:,:,:), dp(:,:,:,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(gid,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(gid,:,:,bccomp:bccomp+ncomp-1),&
                                     ng, is_conservative)
             else
                call mkflux_3d(sop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               fp(:,:,:,:), dp(:,:,:,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(gid,:,:), &
                               the_bc_level(n)%adv_bc_level_array(gid,:,:,bccomp:bccomp+ncomp-1),&
                               ng, is_conservative)
             endif
          end select
       end do
    enddo

    do n = nlevs,2,-1
       do comp = 1, ncomp
          if(is_conservative(comp)) then
             do i = 1, dm
                call ml_edge_restriction_c(flux(n-1,i),comp,flux(n,i),comp, &
                                           mla%mba%rr(n-1,:),i,1)
             enddo
          endif
       enddo

    enddo
    
  end subroutine mkflux


  subroutine mkflux_2d(s,sedgex,sedgey,fluxx,fluxy,umac,vmac,force,divu,lo,dx,dt,is_vel, &
                       phys_bc,adv_bc,ng,is_conservative)

    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_minion

    integer, intent(in) :: lo(:),ng

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
    real(kind=dp_t), intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:)

    real(kind=dp_t),intent(in) :: dt,dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)
    logical        ,intent(in) :: is_vel, is_conservative(:)

    ! Local variables
    real(kind=dp_t), allocatable:: slopex(:,:,:)
    real(kind=dp_t), allocatable:: slopey(:,:,:)

    real(kind=dp_t) hx, hy, dt2, dt4, savg
    real(kind=dp_t) :: abs_eps, eps, umax

    integer :: hi(2)
    integer :: i,j,is,js,ie,je
    integer :: jc,jp
    integer :: comp,ncomp

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

    call slopex_2d(s,slopex,lo,hi,ng,ncomp,adv_bc)
    call slopey_2d(s,slopey,lo,hi,ng,ncomp,adv_bc)

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
    allocate(simhx(lo(1):hi(1)+1,2))

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

    abs_eps = 1.0d-8

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
    !     5. Cycle indices
    !  enddo
    !
    !*************************************
    ! End pseudo code
    !*************************************

    ! loop over components
    do comp=1,ncomp

       jc = 1
       jp = 2

       do j=js-1,je+1

          !******************************************************************
          ! 1. Compute simhx(is:ie+1,j)
          !******************************************************************

          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,jc) = s(i-1,j,comp) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,comp)
             srx(i,jc) = s(i  ,j,comp) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,comp)

             ! add source terms
             if(use_minion) then
                slx(i,jc) = slx(i,jc) + dt2*force(i-1,j,comp)
                srx(i,jc) = srx(i,jc) + dt2*force(i  ,j,comp)
             endif

             ! add divu contribution
             if(use_minion .and. is_conservative(comp)) then
                slx(i,jc) = slx(i,jc) - dt2*s(i-1,j,comp)*divu(i-1,j)
                srx(i,jc) = srx(i,jc) - dt2*s(i  ,j,comp)*divu(i  ,j)
             endif
          end do
          
          ! impose lo side bc's
          if (phys_bc(1,1) .eq. INLET) then
             slx(is,jc) = s(is-1,j,comp)
             srx(is,jc) = s(is-1,j,comp)
          else if (phys_bc(1,1) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slx(is,jc) = ZERO
                srx(is,jc) = ZERO
             else
                slx(is,jc) = srx(is,jc)
             end if
          else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slx(is,jc) = ZERO
                srx(is,jc) = ZERO
             else
                slx(is,jc) = srx(is,jc)
             end if
          else if (phys_bc(1,1) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 1) then
                slx(is,jc) = min(srx(is,jc),ZERO)
                srx(is,jc) = min(srx(is,jc),ZERO)
             else
                slx(is,jc) = srx(is,jc)
             end if
          end if

          ! impose hi side bc's
          if (phys_bc(1,2) .eq. INLET) then
             slx(ie+1,jc) = s(ie+1,j,comp)
             srx(ie+1,jc) = s(ie+1,j,comp)
          else if (phys_bc(1,2) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slx(ie+1,jc) = ZERO
                srx(ie+1,jc) = ZERO
             else
                srx(ie+1,jc) = slx(ie+1,jc)
             end if
          else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slx(ie+1,jc) = ZERO
                srx(ie+1,jc) = ZERO
             else
                srx(ie+1,jc) = slx(ie+1,jc)
             end if
          else if (phys_bc(1,2) .eq. OUTLET) then       
             if (is_vel .and. comp .eq. 1) then
                slx(ie+1,jc) = max(slx(ie+1,jc),ZERO)
                srx(ie+1,jc) = max(slx(ie+1,jc),ZERO)
             else
                srx(ie+1,jc) = slx(ie+1,jc)
             end if
          end if

          do i=is,ie+1
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
                sly(i,jc) = s(i,j-1,comp) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,comp)
                sry(i,jc) = s(i,j  ,comp) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,comp)

                ! add source terms
                if(use_minion) then
                   sly(i,jc) = sly(i,jc) + dt2*force(i,j-1,comp)
                   sry(i,jc) = sry(i,jc) + dt2*force(i,j  ,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   sly(i,jc) = sly(i,jc) - dt2*s(i,j-1,comp)*divu(i,j-1)
                   sry(i,jc) = sry(i,jc) - dt2*s(i,j  ,comp)*divu(i,j  )
                endif
             end do
             
             ! impose lo side bc's
             if(j .eq. js) then
                if (phys_bc(2,1) .eq. INLET) then
                   sly(is-1:ie+1,jc) = s(is-1:ie+1,js-1,comp)
                   sry(is-1:ie+1,jc) = s(is-1:ie+1,js-1,comp)
                else if (phys_bc(2,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      sly(is-1:ie+1,jc) = ZERO
                      sry(is-1:ie+1,jc) = ZERO
                   else
                      sly(is-1:ie+1,jc) = sry(is-1:ie+1,jc)
                   end if
                else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sly(is-1:ie+1,jc) = ZERO
                      sry(is-1:ie+1,jc) = ZERO
                   else
                      sly(is-1:ie+1,jc) = sry(is-1:ie+1,jc)
                   end if
                else if (phys_bc(2,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 2) then
                      sly(is-1:ie+1,jc) = min(sry(is-1:ie+1,jc),ZERO)
                      sry(is-1:ie+1,jc) = min(sry(is-1:ie+1,jc),ZERO)
                   else
                      sly(is-1:ie+1,jc) = sry(is-1:ie+1,js)
                   end if
                end if
             end if

             if(j .eq. je+1) then
                ! impose hi side bc's
                if (phys_bc(2,2) .eq. INLET) then
                   sly(is-1:ie+1,jc) = s(is-1:ie+1,je+1,comp)
                   sry(is-1:ie+1,jc) = s(is-1:ie+1,je+1,comp)
                else if (phys_bc(2,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      sly(is-1:ie+1,jc) = ZERO
                      sry(is-1:ie+1,jc) = ZERO
                   else
                      sry(is-1:ie+1,jc) = sly(is-1:ie+1,jc)
                   end if
                else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sly(is-1:ie+1,jc) = ZERO
                      sry(is-1:ie+1,jc) = ZERO
                   else
                      sry(is-1:ie+1,jc) = sly(is-1:ie+1,jc)
                   end if
                else if (phys_bc(2,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 2) then
                      sly(is-1:ie+1,jc) = max(sly(is-1:ie+1,jc),ZERO)
                      sry(is-1:ie+1,jc) = max(sly(is-1:ie+1,jc),ZERO)
                   else
                      sry(is-1:ie+1,jc) = sly(is-1:ie+1,jc)
                   end if
                end if
             end if

             do i=is-1,ie+1
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
                if(is_conservative(comp)) then
                   sedgely(i) = sly(i,jc) &
                        - (dt2/hx)*(simhx(i+1,jp)*umac(i+1,j-1) - simhx(i,jp)*umac(i,j-1)) &
                        + (dt2/hx)*s(i,j-1,comp)*(umac(i+1,j-1)-umac(i,j-1))
                   sedgery(i) = sry(i,jc) &
                        - (dt2/hx)*(simhx(i+1,jc)*umac(i+1,j  ) - simhx(i,jc)*umac(i,j  )) &
                        + (dt2/hx)*s(i,j  ,comp)*(umac(i+1,j  )-umac(i,j  ))
                else
                   sedgely(i) = sly(i,jc) &
                        - (dt4/hx)*(umac(i+1,j-1)+umac(i,j-1))*(simhx(i+1,jp)-simhx(i,jp))
                   sedgery(i) = sry(i,jc) &
                        - (dt4/hx)*(umac(i+1,j  )+umac(i,j  ))*(simhx(i+1,jc)-simhx(i,jc))
                endif

                ! if use_minion is true, we have already accounted for source terms
                ! in sly and sry; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   sedgely(i) = sedgely(i) + dt2*force(i,j-1,comp)
                   sedgery(i) = sedgery(i) + dt2*force(i,j  ,comp)
                endif

                ! if use_minion .and. is_conservative, we have already accounted for divu
                ! in slx and srx; otherwise, we account for it here
                if(.not.(use_minion .and. is_conservative(comp))) then
                   sedgely(i) = sedgely(i) - dt2*s(i,j-1,comp)*divu(i,j-1)
                   sedgery(i) = sedgery(i) - dt2*s(i,j  ,comp)*divu(i,j  )
                endif

                ! make sedgey by solving Riemann problem
                ! boundary conditions enforced outside of i,j loop
                sedgey(i,j,comp) = merge(sedgely(i),sedgery(i),vmac(i,j) .gt. ZERO)
                savg = HALF*(sedgely(i)+sedgery(i))
                sedgey(i,j,comp) = merge(sedgey(i,j,comp),savg,abs(vmac(i,j)) .gt. eps)
             end do
             
             ! impose lo side bc's
             if(j .eq. js) then
                if (phys_bc(2,1) .eq. INLET) then
                   sedgey(is:ie,js,comp) = s(is:ie,js-1,comp)
                else if (phys_bc(2,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      sedgey(is:ie,js,comp) = ZERO
                   else
                      sedgey(is:ie,js,comp) = sedgery(is:ie)
                   end if
                else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sedgey(is:ie,js,comp) = ZERO
                   else
                      sedgey(is:ie,js,comp) = sedgery(is:ie)
                   end if
                else if (phys_bc(2,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 2) then
                      sedgey(is:ie,js,comp) = min(sedgery(is:ie),ZERO)
                   else
                      sedgey(is:ie,js,comp) = sedgery(is:ie)
                   end if
                end if
             end if

             ! impose hi side bc's
             if(j .eq. je+1) then
                if (phys_bc(2,2) .eq. INLET) then
                   sedgey(is:ie,je+1,comp) = s(is:ie,je+1,comp)
                else if (phys_bc(2,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      sedgey(is:ie,je+1,comp) = ZERO
                   else
                      sedgey(is:ie,je+1,comp) = sedgely(is:ie)
                   end if
                else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sedgey(is:ie,je+1,comp) = ZERO
                   else
                      sedgey(is:ie,je+1,comp) = sedgely(is:ie)
                   end if
                else if (phys_bc(2,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 2) then
                      sedgey(is:ie,je+1,comp) = max(sedgely(is:ie),ZERO)
                   else
                      sedgey(is:ie,je+1,comp) = sedgely(is:ie)
                   end if
                end if
             end if

             do i=is,ie
                ! create fluxes
                if(is_conservative(comp)) then
                   fluxy(i,j,comp) = sedgey(i,j,comp)*vmac(i,j)
                endif
             enddo

          endif ! end if(j .gt. js-1)

          if (j .gt. js) then

             !******************************************************************
             ! 4. Compute sedgex(is:ie+1,j-1)
             !******************************************************************

             do i=is,ie+1
                ! make sedgelx, sedgerx
                if(is_conservative(comp)) then
                   sedgelx(i) = slx(i,jp) &
                        - (dt2/hy)*(simhy(i-1,jc)*vmac(i-1,j) - simhy(i-1,jp)*vmac(i-1,j-1)) &
                        + (dt2/hy)*s(i-1,j-1,comp)*(vmac(i-1,j)-vmac(i-1,j-1))
                   sedgerx(i) = srx(i,jp) &
                        - (dt2/hy)*(simhy(i  ,jc)*vmac(i  ,j) - simhy(i  ,jp)*vmac(i  ,j-1)) &
                        + (dt2/hy)*s(i  ,j-1,comp)*(vmac(i  ,j)-vmac(i  ,j-1))
                else
                   sedgelx(i) = slx(i,jp) &
                        - (dt4/hy)*(vmac(i-1,j)+vmac(i-1,j-1))*(simhy(i-1,jc)-simhy(i-1,jp))
                   sedgerx(i) = srx(i,jp) &
                        - (dt4/hy)*(vmac(i  ,j)+vmac(i  ,j-1))*(simhy(i  ,jc)-simhy(i  ,jp))
                endif

                ! if use_minion is true, we have already accounted for source terms
                ! in slx and srx; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   sedgelx(i) = sedgelx(i) + dt2*force(i-1,j-1,comp)
                   sedgerx(i) = sedgerx(i) + dt2*force(i  ,j-1,comp)
                endif

                ! if use_minion .and. is_conservative, we have already accounted for divu
                ! in slx and srx; otherwise, we account for it here
                if(.not.(use_minion .and. is_conservative(comp))) then
                   sedgelx(i) = sedgelx(i) - dt2*s(i-1,j-1,comp)*divu(i-1,j-1)
                   sedgerx(i) = sedgerx(i) - dt2*s(i  ,j-1,comp)*divu(i  ,j-1)
                endif

                ! make sedgex by solving Riemann problem
                ! boundary conditions enforced outside of i,j loop
                sedgex(i,j-1,comp) = merge(sedgelx(i),sedgerx(i),umac(i,j-1) .gt. ZERO)
                savg = HALF*(sedgelx(i)+sedgerx(i))
                sedgex(i,j-1,comp) = merge(sedgex(i,j-1,comp),savg,abs(umac(i,j-1)) .gt. eps)
             enddo
             
             ! impose lo side bc's
             if (phys_bc(1,1) .eq. INLET) then
                sedgex(is,j-1,comp) = s(is-1,j-1,comp)
             else if (phys_bc(1,1) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(is,j-1,comp) = ZERO
                else
                   sedgex(is,j-1,comp) = sedgerx(is)
                end if
             else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgex(is,j-1,comp) = ZERO
                else
                   sedgex(is,j-1,comp) = sedgerx(is)
                end if
             else if (phys_bc(1,1) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(is,j-1,comp) = min(sedgerx(is),ZERO)
                else
                   sedgex(is,j-1,comp) = sedgerx(is)
                end if
             end if

             ! impose hi side bc's
             if (phys_bc(1,2) .eq. INLET) then
                sedgex(ie+1,j-1,comp) = s(ie+1,j-1,comp)
             else if (phys_bc(1,2) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(ie+1,j-1,comp) = ZERO
                else
                   sedgex(ie+1,j-1,comp) = sedgelx(ie+1)
                end if
             else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgex(ie+1,j-1,comp) = ZERO
                else
                   sedgex(ie+1,j-1,comp) = sedgelx(ie+1)
                end if
             else if (phys_bc(1,2) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(ie+1,j-1,comp) = max(sedgelx(ie+1),ZERO)
                else
                   sedgex(ie+1,j-1,comp) = sedgelx(ie+1)
                end if
             end if

             ! create fluxes
             if(is_conservative(comp)) then
                do i=is,ie+1
                   fluxx(i,j-1,comp) = sedgex(i,j-1,comp)*umac(i,j-1)
                enddo
             endif

          endif ! end if(j .gt. js)

          !******************************************************************
          ! 5. Cycle indices
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

  end subroutine mkflux_2d

  subroutine mkflux_debug_2d(s,sedgex,sedgey,fluxx,fluxy,umac,vmac,force,divu,lo,dx,dt, &
                             is_vel,phys_bc,adv_bc,ng,is_conservative)

    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_minion

    integer, intent(in) :: lo(:),ng

    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::  fluxx(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::  fluxy(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
    real(kind=dp_t), intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:)

    real(kind=dp_t),intent(in) :: dt,dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)
    logical        ,intent(in) :: is_vel, is_conservative(:)

    ! Local variables
    real(kind=dp_t), allocatable:: slopex(:,:,:)
    real(kind=dp_t), allocatable:: slopey(:,:,:)

    real(kind=dp_t) hx, hy, dt2, dt4, savg
    real(kind=dp_t) :: abs_eps, eps, umax

    integer :: hi(2)
    integer :: i,j,is,js,ie,je
    integer :: comp,ncomp

    ! these correspond to s_L^x, etc.
    real(kind=dp_t), allocatable:: slx(:,:),srx(:,:)
    real(kind=dp_t), allocatable:: sly(:,:),sry(:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    real(kind=dp_t), allocatable:: simhx(:,:),simhy(:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    real(kind=dp_t), allocatable:: sedgelx(:,:),sedgerx(:,:)
    real(kind=dp_t), allocatable:: sedgely(:,:),sedgery(:,:)

    ncomp = size(s,dim=3)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,ncomp))

    call slopex_2d(s,slopex,lo,hi,ng,ncomp,adv_bc)
    call slopey_2d(s,slopey,lo,hi,ng,ncomp,adv_bc)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse direction
    allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1))

    abs_eps = 1.0d-8

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

    ! loop over components
    do comp = 1,ncomp

       !******************************************************************
       ! Create s_{\i-\half\e_x}^x, etc.
       !******************************************************************

       ! loop over appropriate x-faces
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j) = s(i-1,j,comp) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,comp)
             srx(i,j) = s(i  ,j,comp) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,comp)

             ! add source terms
             if(use_minion) then
                slx(i,j) = slx(i,j) + dt2*force(i-1,j,comp)
                srx(i,j) = srx(i,j) + dt2*force(i  ,j,comp)
             endif

             ! add divu contribution
             if(use_minion .and. is_conservative(comp)) then
                slx(i,j) = slx(i,j) - dt2*s(i-1,j,comp)*divu(i-1,j)
                srx(i,j) = srx(i,j) - dt2*s(i  ,j,comp)*divu(i  ,j)
             endif
          end do
       end do
       
       ! impose lo side bc's
       if (phys_bc(1,1) .eq. INLET) then
          slx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
          srx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
       else if (phys_bc(1,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slx(is,js-1:je+1) = ZERO
             srx(is,js-1:je+1) = ZERO
          else
             slx(is,js-1:je+1) = srx(is,js-1:je+1)
          end if
       else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slx(is,js-1:je+1) = ZERO
             srx(is,js-1:je+1) = ZERO
          else
             slx(is,js-1:je+1) = srx(is,js-1:je+1)
          end if
       else if (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slx(is,js-1:je+1) = min(srx(is,js-1:je+1),ZERO)
             srx(is,js-1:je+1) = min(srx(is,js-1:je+1),ZERO)
          else
             slx(is,js-1:je+1) = srx(is,js-1:je+1)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(1,2) .eq. INLET) then
          slx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
          srx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
       else if (phys_bc(1,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1) = ZERO
             srx(ie+1,js-1:je+1) = ZERO
          else
             srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
          end if
       else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slx(ie+1,js-1:je+1) = ZERO
             srx(ie+1,js-1:je+1) = ZERO
          else
             srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
          end if
       else if (phys_bc(1,2) .eq. OUTLET) then       
          if (is_vel .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1) = max(slx(ie+1,js-1:je+1),ZERO)
             srx(ie+1,js-1:je+1) = max(slx(ie+1,js-1:je+1),ZERO)
          else
             srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
          end if
       end if

       do j=js-1,je+1
          do i=is,ie+1
             ! make simhx by solving Riemann problem
             simhx(i,j) = merge(slx(i,j),srx(i,j),umac(i,j) .gt. ZERO)
             savg = HALF*(slx(i,j)+srx(i,j))
             simhx(i,j) = merge(simhx(i,j),savg,abs(umac(i,j)) .gt. eps)
          enddo
       enddo

       ! loop over appropriate y-faces
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j) = s(i,j-1,comp) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,comp)
             sry(i,j) = s(i,j  ,comp) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,comp)

             ! add source terms
             if(use_minion) then
                sly(i,j) = sly(i,j) + dt2*force(i,j-1,comp)
                sry(i,j) = sry(i,j) + dt2*force(i,j  ,comp)
             endif

             ! add divu contribution
             if(use_minion .and. is_conservative(comp)) then
                sly(i,j) = sly(i,j) - dt2*s(i,j-1,comp)*divu(i,j-1)
                sry(i,j) = sry(i,j) - dt2*s(i,j  ,comp)*divu(i,j  )
             endif
          end do
       end do
       
       ! impose lo side bc's
       if (phys_bc(2,1) .eq. INLET) then
          sly(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
          sry(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
       else if (phys_bc(2,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,js) = ZERO
             sry(is-1:ie+1,js) = ZERO
          else
             sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
          end if
       else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sly(is-1:ie+1,js) = ZERO
             sry(is-1:ie+1,js) = ZERO
          else
             sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
          end if
       else if (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,js) = min(sry(is-1:ie+1,js),ZERO)
             sry(is-1:ie+1,js) = min(sry(is-1:ie+1,js),ZERO)
          else
             sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(2,2) .eq. INLET) then
          sly(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
          sry(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
       else if (phys_bc(2,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1) = ZERO
             sry(is-1:ie+1,je+1) = ZERO
          else
             sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
          end if
       else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sly(is-1:ie+1,je+1) = ZERO
             sry(is-1:ie+1,je+1) = ZERO
          else
             sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
          end if
       else if (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1) = max(sly(is-1:ie+1,je+1),ZERO)
             sry(is-1:ie+1,je+1) = max(sly(is-1:ie+1,je+1),ZERO)
          else
             sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
          end if
       end if

       do j=js,je+1
          do i=is-1,ie+1
             ! make simhy by solving Riemann problem
             simhy(i,j) = merge(sly(i,j),sry(i,j),vmac(i,j) .gt. ZERO)
             savg = HALF*(sly(i,j)+sry(i,j))
             simhy(i,j) = merge(simhy(i,j),savg,abs(vmac(i,j)) .gt. eps)
          enddo
       enddo

       !******************************************************************
       ! Create sedgelx, etc.
       !******************************************************************

       ! loop over appropriate x-faces
       do j=js,je
          do i=is,ie+1
             ! make sedgelx, sedgerx
             if(is_conservative(comp)) then
                sedgelx(i,j) = slx(i,j) &
                     - (dt2/hy)*(simhy(i-1,j+1)*vmac(i-1,j+1) - simhy(i-1,j)*vmac(i-1,j)) &
                     + (dt2/hy)*s(i-1,j,comp)*(vmac(i-1,j+1)-vmac(i-1,j))
                sedgerx(i,j) = srx(i,j) &
                     - (dt2/hy)*(simhy(i  ,j+1)*vmac(i  ,j+1) - simhy(i  ,j)*vmac(i  ,j)) &
                     + (dt2/hy)*s(i  ,j,comp)*(vmac(i  ,j+1)-vmac(i  ,j))
             else
                sedgelx(i,j) = slx(i,j) &
                     - (dt4/hy)*(vmac(i-1,j+1)+vmac(i-1,j))*(simhy(i-1,j+1)-simhy(i-1,j))
                sedgerx(i,j) = srx(i,j) &
                     - (dt4/hy)*(vmac(i  ,j+1)+vmac(i  ,j))*(simhy(i  ,j+1)-simhy(i  ,j))
             endif

             ! if use_minion is true, we have already accounted for source terms
             ! in slx and srx; otherwise, we need to account for them here.
             if(.not. use_minion) then
                sedgelx(i,j) = sedgelx(i,j) + dt2*force(i-1,j,comp)
                sedgerx(i,j) = sedgerx(i,j) + dt2*force(i  ,j,comp)
             endif

             ! if use_minion .and. is_conservative, we have already accounted for divu
             ! in slx and srx; otherwise, we account for it here
             if(.not.(use_minion .and. is_conservative(comp))) then
                sedgelx(i,j) = sedgelx(i,j) - dt2*s(i-1,j,comp)*divu(i-1,j)
                sedgerx(i,j) = sedgerx(i,j) - dt2*s(i  ,j,comp)*divu(i  ,j)
             endif

             ! make sedgex by solving Riemann problem
             ! boundary conditions enforced outside of i,j loop
             sedgex(i,j,comp) = merge(sedgelx(i,j),sedgerx(i,j),umac(i,j) .gt. ZERO)
             savg = HALF*(sedgelx(i,j)+sedgerx(i,j))
             sedgex(i,j,comp) = merge(sedgex(i,j,comp),savg,abs(umac(i,j)) .gt. eps)
          enddo
       enddo
       
       ! impose lo side bc's
       if (phys_bc(1,1) .eq. INLET) then
          sedgex(is,js:je,comp) = s(is-1,js:je,comp)
       else if (phys_bc(1,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(is,js:je,comp) = ZERO
          else
             sedgex(is,js:je,comp) = sedgerx(is,js:je)
          end if
       else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sedgex(is,js:je,comp) = ZERO
          else
             sedgex(is,js:je,comp) = sedgerx(is,js:je)
          end if
       else if (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(is,js:je,comp) = min(sedgerx(is,js:je),ZERO)
          else
             sedgex(is,js:je,comp) = sedgerx(is,js:je)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(1,2) .eq. INLET) then
          sedgex(ie+1,js:je,comp) = s(ie+1,js:je,comp)
       else if (phys_bc(1,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(ie+1,js:je,comp) = ZERO
          else
             sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
          end if
       else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sedgex(ie+1,js:je,comp) = ZERO
          else
             sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
          end if
       else if (phys_bc(1,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(ie+1,js:je,comp) = max(sedgelx(ie+1,js:je),ZERO)
          else
             sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
          end if
       end if

       ! create fluxes
       do j=js,je
          do i=is,ie+1
             if(is_conservative(comp)) then
                fluxx(i,j,comp) = sedgex(i,j,comp)*umac(i,j)
             endif
          enddo
       enddo

       ! loop over appropriate y-faces
       do j=js,je+1
          do i=is,ie
             ! make sedgely, sedgery
             if(is_conservative(comp)) then
                sedgely(i,j) = sly(i,j) &
                     - (dt2/hx)*(simhx(i+1,j-1)*umac(i+1,j-1) - simhx(i,j-1)*umac(i,j-1)) &
                     + (dt2/hx)*s(i,j-1,comp)*(umac(i+1,j-1)-umac(i,j-1))
                sedgery(i,j) = sry(i,j) &
                     - (dt2/hx)*(simhx(i+1,j  )*umac(i+1,j  ) - simhx(i,j  )*umac(i,j  )) &
                     + (dt2/hx)*s(i,j  ,comp)*(umac(i+1,j  )-umac(i,j  ))
             else
                sedgely(i,j) = sly(i,j) &
                     - (dt4/hx)*(umac(i+1,j-1)+umac(i,j-1))*(simhx(i+1,j-1)-simhx(i,j-1))
                sedgery(i,j) = sry(i,j) &
                     - (dt4/hx)*(umac(i+1,j  )+umac(i,j  ))*(simhx(i+1,j  )-simhx(i,j  ))
             endif

             ! if use_minion is true, we have already accounted for source terms
             ! in sly and sry; otherwise, we need to account for them here.
             if(.not. use_minion) then
                sedgely(i,j) = sedgely(i,j) + dt2*force(i,j-1,comp)
                sedgery(i,j) = sedgery(i,j) + dt2*force(i,j  ,comp)
             endif

             ! if use_minion .and. is_conservative, we have already accounted for divu
             ! in sly and sry; otherwise, we account for it here
             if(.not.(use_minion .and. is_conservative(comp))) then
                sedgely(i,j) = sedgely(i,j) - dt2*s(i,j-1,comp)*divu(i,j-1)
                sedgery(i,j) = sedgery(i,j) - dt2*s(i,j  ,comp)*divu(i,j  )
             endif

             ! make sedgey by solving Riemann problem
             ! boundary conditions enforced outside of i,j loop
             sedgey(i,j,comp) = merge(sedgely(i,j),sedgery(i,j),vmac(i,j) .gt. ZERO)
             savg = HALF*(sedgely(i,j)+sedgery(i,j))
             sedgey(i,j,comp) = merge(sedgey(i,j,comp),savg,abs(vmac(i,j)) .gt. eps)
          enddo
       enddo
       
       ! impose lo side bc's
       if (phys_bc(2,1) .eq. INLET) then
          sedgey(is:ie,js,comp) = s(is:ie,js-1,comp)
       else if (phys_bc(2,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(is:ie,js,comp) = ZERO
          else
             sedgey(is:ie,js,comp) = sedgery(is:ie,js)
          end if
       else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sedgey(is:ie,js,comp) = ZERO
          else
             sedgey(is:ie,js,comp) = sedgery(is:ie,js)
          end if
       else if (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(is:ie,js,comp) = min(sedgery(is:ie,js),ZERO)
          else
             sedgey(is:ie,js,comp) = sedgery(is:ie,js)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(2,2) .eq. INLET) then
          sedgey(is:ie,je+1,comp) = s(is:ie,je+1,comp)
       else if (phys_bc(2,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(is:ie,je+1,comp) = ZERO
          else
             sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
          end if
       else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sedgey(is:ie,je+1,comp) = ZERO
          else
             sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
          end if
       else if (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(is:ie,je+1,comp) = max(sedgely(is:ie,je+1),ZERO)
          else
             sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
          end if
       end if

       ! create fluxes
       do j=js,je+1
          do i=is,ie
             if(is_conservative(comp)) then
                fluxy(i,j,comp) = sedgey(i,j,comp)*vmac(i,j)
             endif
          enddo
       enddo

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

  end subroutine mkflux_debug_2d

  subroutine mkflux_3d(s,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz,umac,vmac,wmac,force, &
                       divu,lo,dx,dt,is_vel,phys_bc,adv_bc,ng,is_conservative)

    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_minion

    integer, intent(in) :: lo(:),ng

    real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
    real(kind=dp_t),intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)
    real(kind=dp_t),intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)

    real(kind=dp_t),intent(in) :: dt,dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)
    logical        ,intent(in) :: is_vel, is_conservative(:)

    ! Local variables
    real(kind=dp_t), allocatable:: slopex(:,:,:,:)
    real(kind=dp_t), allocatable:: slopey(:,:,:,:)
    real(kind=dp_t), allocatable:: slopez(:,:,:,:)

    real(kind=dp_t) hx, hy, hz, dt2, dt3, dt4, dt6, savg
    real(kind=dp_t) :: abs_eps, eps, umax

    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke
    integer :: kc,kp
    integer :: comp,ncomp

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
       call slopex_2d(s(:,:,k,:),slopex(:,:,k,:),lo,hi,ng,ncomp,adv_bc)
       call slopey_2d(s(:,:,k,:),slopey(:,:,k,:),lo,hi,ng,ncomp,adv_bc)
    end do
    call slopez_3d(s,slopez,lo,hi,ng,ncomp,adv_bc)

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

    abs_eps = 1.0d-8

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
    !  do k=ks-1,ke+1
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
    !     13. Cycle indices
    !  enddo
    !
    !*************************************
    ! End pseudo code
    !*************************************

    ! loop over components
    do comp=1,ncomp

       kc = 1
       kp = 2

       do k=ks-1,ke+1

          !******************************************************************
          ! 1. Compute simhx (is  :ie+1,js-1:je+1,k)
          !******************************************************************

          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,kc) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,comp)
                srx(i,j,kc) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i  ,j,k,comp)

                ! add source terms
                if(use_minion) then
                   slx(i,j,kc) = slx(i,j,kc) + dt2*force(i-1,j,k,comp)
                   srx(i,j,kc) = srx(i,j,kc) + dt2*force(i  ,j,k,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   slx(i,j,kc) = slx(i,j,kc) - dt2*s(i-1,j,k,comp)*divu(i-1,j,k)
                   srx(i,j,kc) = srx(i,j,kc) - dt2*s(i  ,j,k,comp)*divu(i  ,j,k)
                endif
             end do
          end do
          
          ! impose lo side bc's
          if (phys_bc(1,1) .eq. INLET) then
             slx(is,js-1:je+1,kc) = s(is-1,js-1:je+1,k,comp)
             srx(is,js-1:je+1,kc) = s(is-1,js-1:je+1,k,comp)
          else if (phys_bc(1,1) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slx(is,js-1:je+1,kc) = ZERO
                srx(is,js-1:je+1,kc) = ZERO
             else
                slx(is,js-1:je+1,kc) = srx(is,js-1:je+1,kc)
             end if
          else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slx(is,js-1:je+1,kc) = ZERO
                srx(is,js-1:je+1,kc) = ZERO
             else
                slx(is,js-1:je+1,kc) = srx(is,js-1:je+1,kc)
             end if
          else if (phys_bc(1,1) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 1) then
                slx(is,js-1:je+1,kc) = min(srx(is,js-1:je+1,kc),ZERO)
                srx(is,js-1:je+1,kc) = min(srx(is,js-1:je+1,kc),ZERO)
             else
                slx(is,js-1:je+1,kc) = srx(is,js-1:je+1,kc)
             end if
          end if

          ! impose hi side bc's
          if (phys_bc(1,2) .eq. INLET) then
             slx(ie+1,js-1:je+1,kc) = s(ie+1,js-1:je+1,k,comp)
             srx(ie+1,js-1:je+1,kc) = s(ie+1,js-1:je+1,k,comp)
          else if (phys_bc(1,2) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slx(ie+1,js-1:je+1,kc) = ZERO
                srx(ie+1,js-1:je+1,kc) = ZERO
             else
                srx(ie+1,js-1:je+1,kc) = slx(ie+1,js-1:je+1,kc)
             end if
          else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slx(ie+1,js-1:je+1,kc) = ZERO
                srx(ie+1,js-1:je+1,kc) = ZERO
             else
                srx(ie+1,js-1:je+1,kc) = slx(ie+1,js-1:je+1,kc)
             end if
          else if (phys_bc(1,2) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 1) then
                slx(ie+1,js-1:je+1,kc) = max(slx(ie+1,js-1:je+1,kc),ZERO)
                srx(ie+1,js-1:je+1,kc) = max(slx(ie+1,js-1:je+1,kc),ZERO)
             else
                srx(ie+1,js-1:je+1,kc) = slx(ie+1,js-1:je+1,kc)
             end if
          end if

          do j=js-1,je+1
             do i=is,ie+1
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
                sly(i,j,kc) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,comp)
                sry(i,j,kc) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j  ,k,comp)

                ! add source terms
                if(use_minion) then
                   sly(i,j,kc) = sly(i,j,kc) + dt2*force(i,j-1,k,comp)
                   sry(i,j,kc) = sry(i,j,kc) + dt2*force(i,j  ,k,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   sly(i,j,kc) = sly(i,j,kc) - dt2*s(i,j-1,k,comp)*divu(i,j-1,k)
                   sry(i,j,kc) = sry(i,j,kc) - dt2*s(i,j  ,k,comp)*divu(i,j  ,k)
                endif
             end do
          end do
          
          ! impose lo side bc's
          if (phys_bc(2,1) .eq. INLET) then
             sly(is-1:ie+1,js,kc) = s(is-1:ie+1,js-1,k,comp)
             sry(is-1:ie+1,js,kc) = s(is-1:ie+1,js-1,k,comp)
          else if (phys_bc(2,1) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sly(is-1:ie+1,js,kc) = ZERO
                sry(is-1:ie+1,js,kc) = ZERO
             else
                sly(is-1:ie+1,js,kc) = sry(is-1:ie+1,js,kc)
             end if
          else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                sly(is-1:ie+1,js,kc) = ZERO
                sry(is-1:ie+1,js,kc) = ZERO
             else
                sly(is-1:ie+1,js,kc) = sry(is-1:ie+1,js,kc)
             end if
          else if (phys_bc(2,1) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 2) then
                sly(is-1:ie+1,js,kc) = min(sry(is-1:ie+1,js,kc),ZERO)
                sry(is-1:ie+1,js,kc) = min(sry(is-1:ie+1,js,kc),ZERO)
             else
                sly(is-1:ie+1,js,kc) = sry(is-1:ie+1,js,kc)
             end if
          end if

          ! impose hi side bc's
          if (phys_bc(2,2) .eq. INLET) then
             sly(is-1:ie+1,je+1,kc) = s(is-1:ie+1,je+1,k,comp)
             sry(is-1:ie+1,je+1,kc) = s(is-1:ie+1,je+1,k,comp)
          else if (phys_bc(2,2) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sly(is-1:ie+1,je+1,kc) = ZERO
                sry(is-1:ie+1,je+1,kc) = ZERO
             else
                sry(is-1:ie+1,je+1,kc) = sly(is-1:ie+1,je+1,kc)
             end if
          else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                sly(is-1:ie+1,je+1,kc) = ZERO
                sry(is-1:ie+1,je+1,kc) = ZERO
             else
                sry(is-1:ie+1,je+1,kc) = sly(is-1:ie+1,je+1,kc)
             end if
          else if (phys_bc(2,2) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 2) then
                sly(is-1:ie+1,je+1,kc) = max(sly(is-1:ie+1,je+1,kc),ZERO)
                sry(is-1:ie+1,je+1,kc) = max(sly(is-1:ie+1,je+1,kc),ZERO)
             else
                sry(is-1:ie+1,je+1,kc) = sly(is-1:ie+1,je+1,kc)
             end if
          end if

          do j=js,je+1
             do i=is-1,ie+1
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
                if(is_conservative(comp)) then
                   slxy(i,j,kc) = slx(i,j,kc) - (dt3/hy)*(simhy(i-1,j+1,kc)*vmac(i-1,j+1,k) - simhy(i-1,j,kc)*vmac(i-1,j,k))
                   srxy(i,j,kc) = srx(i,j,kc) - (dt3/hy)*(simhy(i  ,j+1,kc)*vmac(i  ,j+1,k) - simhy(i  ,j,kc)*vmac(i  ,j,k))
                else
                   slxy(i,j,kc) = slx(i,j,kc) - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhy(i-1,j+1,kc)-simhy(i-1,j,kc))
                   srxy(i,j,kc) = srx(i,j,kc) - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k))*(simhy(i  ,j+1,kc)-simhy(i  ,j,kc))
                endif
             end do
          end do
          
          ! impose lo side bc's
          if (phys_bc(1,1) .eq. INLET) then
             slxy(is,js:je,kc) = s(is-1,js:je,k,comp)
             srxy(is,js:je,kc) = s(is-1,js:je,k,comp)
          else if (phys_bc(1,1) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slxy(is,js:je,kc) = ZERO
                srxy(is,js:je,kc) = ZERO
             else
                slxy(is,js:je,kc) = srxy(is,js:je,kc)
             end if
          else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slxy(is,js:je,kc) = ZERO
                srxy(is,js:je,kc) = ZERO
             else
                slxy(is,js:je,kc) = srxy(is,js:je,kc)
             end if
          else if (phys_bc(1,1) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 1) then
                slxy(is,js:je,kc) = min(srxy(is,js:je,kc),ZERO)
                srxy(is,js:je,kc) = min(srxy(is,js:je,kc),ZERO)
             else
                slxy(is,js:je,kc) = srxy(is,js:je,kc)
             end if
          end if

          ! impose hi side bc's
          if (phys_bc(1,2) .eq. INLET) then
             slxy(ie+1,js:je,kc) = s(ie+1,js:je,k,comp)
             srxy(ie+1,js:je,kc) = s(ie+1,js:je,k,comp)
          else if (phys_bc(1,2) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                slxy(ie+1,js:je,kc) = ZERO
                srxy(ie+1,js:je,kc) = ZERO
             else
                srxy(ie+1,js:je,kc) = slxy(ie+1,js:je,kc)
             end if
          else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slxy(ie+1,js:je,kc) = ZERO
                srxy(ie+1,js:je,kc) = ZERO
             else
                srxy(ie+1,js:je,kc) = slxy(ie+1,js:je,kc)
             end if
          else if (phys_bc(1,2) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 1) then
                slxy(ie+1,js:je,kc) = max(slxy(ie+1,js:je,kc),ZERO)
                srxy(ie+1,js:je,kc) = max(slxy(ie+1,js:je,kc),ZERO)
             else
                srxy(ie+1,js:je,kc) = slxy(ie+1,js:je,kc)
             end if
          end if

          do j=js,je
             do i=is,ie+1
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
                if(is_conservative(comp)) then
                   slyx(i,j,kc) = sly(i,j,kc) - (dt3/hx)*(simhx(i+1,j-1,kc)*umac(i+1,j-1,k) - simhx(i,j-1,kc)*umac(i,j-1,k))
                   sryx(i,j,kc) = sry(i,j,kc) - (dt3/hx)*(simhx(i+1,j  ,kc)*umac(i+1,j  ,k) - simhx(i,j  ,kc)*umac(i,j  ,k))
                else
                   slyx(i,j,kc) = sly(i,j,kc) - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhx(i+1,j-1,kc)-simhx(i,j-1,kc))
                   sryx(i,j,kc) = sry(i,j,kc) - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k))*(simhx(i+1,j  ,kc)-simhx(i,j  ,kc))
                endif
             end do
          end do
          
          ! impose lo side bc's
          if (phys_bc(2,1) .eq. INLET) then
             slyx(is:ie,js,kc) = s(is:ie,js-1,k,comp)
             sryx(is:ie,js,kc) = s(is:ie,js-1,k,comp)
          else if (phys_bc(2,1) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                slyx(is:ie,js,kc) = ZERO
                sryx(is:ie,js,kc) = ZERO
             else
                slyx(is:ie,js,kc) = sryx(is:ie,js,kc)
             end if
          else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slyx(is:ie,js,kc) = ZERO
                sryx(is:ie,js,kc) = ZERO
             else
                slyx(is:ie,js,kc) = sryx(is:ie,js,kc)
             end if
          else if (phys_bc(2,1) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 2) then
                slyx(is:ie,js,kc) = min(sryx(is:ie,js,kc),ZERO)
                sryx(is:ie,js,kc) = min(sryx(is:ie,js,kc),ZERO)
             else
                slyx(is:ie,js,kc) = sryx(is:ie,js,kc)
             end if
          end if

          ! impose hi side bc's
          if (phys_bc(2,2) .eq. INLET) then
             slyx(is:ie,je+1,kc) = s(is:ie,je+1,k,comp)
             sryx(is:ie,je+1,kc) = s(is:ie,je+1,k,comp)
          else if (phys_bc(2,2) .eq. SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                slyx(is:ie,je+1,kc) = ZERO
                sryx(is:ie,je+1,kc) = ZERO
             else
                sryx(is:ie,je+1,kc) = slyx(is:ie,je+1,kc)
             end if
          else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
             if (is_vel) then
                slyx(is:ie,je+1,kc) = ZERO
                sryx(is:ie,je+1,kc) = ZERO
             else
                sryx(is:ie,je+1,kc) = slyx(is:ie,je+1,kc)
             end if
          else if (phys_bc(2,2) .eq. OUTLET) then
             if (is_vel .and. comp .eq. 2) then
                slyx(is:ie,je+1,kc) = max(slyx(is:ie,je+1,kc),ZERO)
                sryx(is:ie,je+1,kc) = max(slyx(is:ie,je+1,kc),ZERO)
             else
                sryx(is:ie,je+1,kc) = slyx(is:ie,je+1,kc)
             end if
          end if

          do j=js,je+1
             do i=is,ie
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
                   slz(i,j,kc) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,comp)
                   srz(i,j,kc) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k  ,comp)

                   ! add source terms
                   if(use_minion) then
                      slz(i,j,kc) = slz(i,j,kc) + dt2*force(i,j,k-1,comp)
                      srz(i,j,kc) = srz(i,j,kc) + dt2*force(i,j,k  ,comp)
                   endif

                   ! add divu contribution
                   if(use_minion .and. is_conservative(comp)) then
                      slz(i,j,kc) = slz(i,j,kc) - dt2*s(i,j,k-1,comp)*divu(i,j,k-1)
                      srz(i,j,kc) = srz(i,j,kc) - dt2*s(i,j,k  ,comp)*divu(i,j,k  )
                   endif
                end do
             end do
             
             ! impose lo side bc's
             if(k .eq. ks) then
                if (phys_bc(3,1) .eq. INLET) then
                   slz(is-1:ie+1,js-1:je+1,kc) = s(is-1:ie+1,js-1:je+1,ks-1,comp)
                   srz(is-1:ie+1,js-1:je+1,kc) = s(is-1:ie+1,js-1:je+1,ks-1,comp)
                else if (phys_bc(3,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slz(is-1:ie+1,js-1:je+1,kc) = ZERO
                      srz(is-1:ie+1,js-1:je+1,kc) = ZERO
                   else
                      slz(is-1:ie+1,js-1:je+1,kc) = srz(is-1:ie+1,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slz(is-1:ie+1,js-1:je+1,kc) = ZERO
                      srz(is-1:ie+1,js-1:je+1,kc) = ZERO
                   else
                      slz(is-1:ie+1,js-1:je+1,kc) = srz(is-1:ie+1,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slz(is-1:ie+1,js-1:je+1,kc) = min(srz(is-1:ie+1,js-1:je+1,kc),ZERO)
                      srz(is-1:ie+1,js-1:je+1,kc) = min(srz(is-1:ie+1,js-1:je+1,kc),ZERO)
                   else
                      slz(is-1:ie+1,js-1:je+1,kc) = srz(is-1:ie+1,js-1:je+1,kc)
                   end if
                end if
             end if

             ! impose hi side bc's
             if(k .eq. ke+1) then
                if (phys_bc(3,2) .eq. INLET) then
                   slz(is-1:ie+1,js-1:je+1,kc) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
                   srz(is-1:ie+1,js-1:je+1,kc) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
                else if (phys_bc(3,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slz(is-1:ie+1,js-1:je+1,kc) = ZERO
                      srz(is-1:ie+1,js-1:je+1,kc) = ZERO
                   else
                      srz(is-1:ie+1,js-1:je+1,kc) = slz(is-1:ie+1,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slz(is-1:ie+1,js-1:je+1,kc) = ZERO
                      srz(is-1:ie+1,js-1:je+1,kc) = ZERO
                   else
                      srz(is-1:ie+1,js-1:je+1,kc) = slz(is-1:ie+1,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slz(is-1:ie+1,js-1:je+1,kc) = max(slz(is-1:ie+1,js-1:je+1,kc),ZERO)
                      srz(is-1:ie+1,js-1:je+1,kc) = max(slz(is-1:ie+1,js-1:je+1,kc),ZERO)
                   else
                      srz(is-1:ie+1,js-1:je+1,kc) = slz(is-1:ie+1,js-1:je+1,kc)
                   end if
                end if
             end if

             do j=js-1,je+1
                do i=is-1,ie+1
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
                   if(is_conservative(comp)) then
                      sedgelz(i,j) = slz(i,j,kc) &
                           - (dt2/hx)*(simhxy(i+1,j,kp)*umac(i+1,j,k-1) - simhxy(i,j,kp)*umac(i,j,k-1)) &
                           - (dt2/hy)*(simhyx(i,j+1,kp)*vmac(i,j+1,k-1) - simhyx(i,j,kp)*vmac(i,j,k-1)) &
                           + (dt2/hx)*s(i,j,k-1,comp)*(umac(i+1,j  ,k-1)-umac(i,j,k-1)) &
                           + (dt2/hy)*s(i,j,k-1,comp)*(vmac(i  ,j+1,k-1)-vmac(i,j,k-1))
                      sedgerz(i,j) = srz(i,j,kc) &
                           - (dt2/hx)*(simhxy(i+1,j,kc)*umac(i+1,j,k  ) - simhxy(i,j,kc)*umac(i,j,k  )) &
                           - (dt2/hy)*(simhyx(i,j+1,kc)*vmac(i,j+1,k  ) - simhyx(i,j,kc)*vmac(i,j,k  )) &
                           + (dt2/hx)*s(i,j,k  ,comp)*(umac(i+1,j  ,k  )-umac(i,j,k  )) &
                           + (dt2/hy)*s(i,j,k  ,comp)*(vmac(i  ,j+1,k  )-vmac(i,j,k  ))
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
                      sedgelz(i,j) = sedgelz(i,j) + dt2*force(i,j,k-1,comp)
                      sedgerz(i,j) = sedgerz(i,j) + dt2*force(i,j,k  ,comp)
                   endif

                   ! if use_minion .and. is_conservative, we have already accounted for divu
                   ! in slz and srz; otherwise, we account for it here
                   if(.not.(use_minion .and. is_conservative(comp))) then
                      sedgelz(i,j) = sedgelz(i,j) - dt2*s(i,j,k-1,comp)*divu(i,j,k-1)
                      sedgerz(i,j) = sedgerz(i,j) - dt2*s(i,j,k  ,comp)*divu(i,j,k  )
                   endif

                   ! make sedgez by solving Riemann problem
                   ! boundary conditions enforced outside of i,j,k loop
                   sedgez(i,j,k,comp) = merge(sedgelz(i,j),sedgerz(i,j),wmac(i,j,k) .gt. ZERO)
                   savg = HALF*(sedgelz(i,j)+sedgerz(i,j))
                   sedgez(i,j,k,comp) = merge(sedgez(i,j,k,comp),savg,abs(wmac(i,j,k)) .gt. eps)
                enddo
             enddo
             
             ! impose lo side bc's
             if(k .eq. ks) then
                if (phys_bc(3,1) .eq. INLET) then
                   sedgez(is:ie,js:je,ks,comp) = s(is:ie,js:je,ks-1,comp)
                else if (phys_bc(3,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      sedgez(is:ie,js:je,ks,comp) = ZERO
                   else
                      sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je)
                   end if
                else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sedgez(is:ie,js:je,ks,comp) = ZERO
                   else
                      sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je)
                   end if
                else if (phys_bc(3,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      sedgez(is:ie,js:je,ks,comp) = min(sedgerz(is:ie,js:je),ZERO)
                   else
                      sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je)
                   end if
                end if
             end if

             ! impose hi side bc's
             if(k .eq. ke+1) then
                if (phys_bc(3,2) .eq. INLET) then
                   sedgez(is:ie,js:je,ke+1,comp) = s(is:ie,js:je,ke+1,comp)
                else if (phys_bc(3,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      sedgez(is:ie,js:je,ke+1,comp) = ZERO
                   else
                      sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je)
                   end if
                else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      sedgez(is:ie,js:je,ke+1,comp) = ZERO
                   else
                      sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je)
                   end if
                else if (phys_bc(3,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      sedgez(is:ie,js:je,ke+1,comp) = max(sedgelz(is:ie,js:je),ZERO)
                   else
                      sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je)
                   end if
                end if
             end if

             ! create fluxes
             if (is_conservative(comp)) then
                do j=js,je
                   do i=is,ie
                      fluxz(i,j,k,comp) = sedgez(i,j,k,comp)*wmac(i,j,k)
                   enddo
                enddo
             endif

             !******************************************************************
             ! 7. Compute simhzx(is  :ie  ,js-1:je+1,k)
             !******************************************************************

             do j=js-1,je+1
                do i=is,ie
                   ! make slzx, srzx by updating 1D extrapolation
                   if(is_conservative(comp)) then
                      slzx(i,j,kc) = slz(i,j,kc) - (dt3/hx)*(simhx(i+1,j,kp)*umac(i+1,j,k-1) - simhx(i,j,kp)*umac(i,j,k-1))
                      srzx(i,j,kc) = srz(i,j,kc) - (dt3/hx)*(simhx(i+1,j,kc)*umac(i+1,j,k  ) - simhx(i,j,kc)*umac(i,j,k  ))
                   else
                      slzx(i,j,kc) = slz(i,j,kc) - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhx(i+1,j,kp)-simhx(i,j,kp))
                      srzx(i,j,kc) = srz(i,j,kc) - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhx(i+1,j,kc)-simhx(i,j,kc))
                   endif
                end do
             end do
             
             ! impose lo side bc's
             if(k .eq. ks) then
                if (phys_bc(3,1) .eq. INLET) then
                   slzx(is:ie,js-1:je+1,kc) = s(is:ie,js-1:je+1,ks-1,comp)
                   srzx(is:ie,js-1:je+1,kc) = s(is:ie,js-1:je+1,ks-1,comp)
                else if (phys_bc(3,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzx(is:ie,js-1:je+1,kc) = ZERO
                      srzx(is:ie,js-1:je+1,kc) = ZERO
                   else
                      slzx(is:ie,js-1:je+1,kc) = srzx(is:ie,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slzx(is:ie,js-1:je+1,kc) = ZERO
                      srzx(is:ie,js-1:je+1,kc) = ZERO
                   else
                      slzx(is:ie,js-1:je+1,kc) = srzx(is:ie,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slzx(is:ie,js-1:je+1,kc) = min(srzx(is:ie,js-1:je+1,kc),ZERO)
                      srzx(is:ie,js-1:je+1,kc) = min(srzx(is:ie,js-1:je+1,kc),ZERO)
                   else
                      slzx(is:ie,js-1:je+1,kc) = srzx(is:ie,js-1:je+1,kc)
                   end if
                end if
             end if

             ! impose hi side bc's
             if(k .eq. ke+1) then
                if (phys_bc(3,2) .eq. INLET) then
                   slzx(is:ie,js-1:je+1,kc) = s(is:ie,js-1:je+1,ke+1,comp)
                   srzx(is:ie,js-1:je+1,kc) = s(is:ie,js-1:je+1,ke+1,comp)
                else if (phys_bc(3,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzx(is:ie,js-1:je+1,kc) = ZERO
                      srzx(is:ie,js-1:je+1,kc) = ZERO
                   else
                      srzx(is:ie,js-1:je+1,kc) = slzx(is:ie,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slzx(is:ie,js-1:je+1,kc) = ZERO
                      srzx(is:ie,js-1:je+1,kc) = ZERO
                   else
                      srzx(is:ie,js-1:je+1,kc) = slzx(is:ie,js-1:je+1,kc)
                   end if
                else if (phys_bc(3,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slzx(is:ie,js-1:je+1,kc) = max(slzx(is:ie,js-1:je+1,kc),ZERO)
                      srzx(is:ie,js-1:je+1,kc) = max(slzx(is:ie,js-1:je+1,kc),ZERO)
                   else
                      srzx(is:ie,js-1:je+1,kc) = slzx(is:ie,js-1:je+1,kc)
                   end if
                end if
             end if

             do j=js-1,je+1
                do i=is,ie
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
                   if(is_conservative(comp)) then
                      slzy(i,j,kc) = slz(i,j,kc) - (dt3/hy)*(simhy(i,j+1,kp)*vmac(i,j+1,k-1) - simhy(i,j,kp)*vmac(i,j,k-1))
                      srzy(i,j,kc) = srz(i,j,kc) - (dt3/hy)*(simhy(i,j+1,kc)*vmac(i,j+1,k  ) - simhy(i,j,kc)*vmac(i,j,k  ))
                   else
                      slzy(i,j,kc) = slz(i,j,kc) - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhy(i,j+1,kp)-simhy(i,j,kp))
                      srzy(i,j,kc) = srz(i,j,kc) - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhy(i,j+1,kc)-simhy(i,j,kc))
                   endif
                end do
             end do
             
             ! impose lo side bc's
             if(k .eq. ks) then
                if (phys_bc(3,1) .eq. INLET) then
                   slzy(is-1:ie+1,js:je,kc) = s(is-1:ie+1,js:je,ks-1,comp)
                   srzy(is-1:ie+1,js:je,kc) = s(is-1:ie+1,js:je,ks-1,comp)
                else if (phys_bc(3,1) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzy(is-1:ie+1,js:je,kc) = ZERO
                      srzy(is-1:ie+1,js:je,kc) = ZERO
                   else
                      slzy(is-1:ie+1,js:je,kc) = srzy(is-1:ie+1,js:je,kc)
                   end if
                else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slzy(is-1:ie+1,js:je,kc) = ZERO
                      srzy(is-1:ie+1,js:je,kc) = ZERO
                   else
                      slzy(is-1:ie+1,js:je,kc) = srzy(is-1:ie+1,js:je,kc)
                   end if
                else if (phys_bc(3,1) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slzy(is-1:ie+1,js:je,kc) = min(srzy(is-1:ie+1,js:je,kc),ZERO)
                      srzy(is-1:ie+1,js:je,kc) = min(srzy(is-1:ie+1,js:je,kc),ZERO)
                   else
                      slzy(is-1:ie+1,js:je,kc) = srzy(is-1:ie+1,js:je,kc)
                   end if
                end if
             end if

             if(k .eq. ke+1) then
                ! impose hi side bc's
                if (phys_bc(3,2) .eq. INLET) then
                   slzy(is-1:ie+1,js:je,kc) = s(is-1:ie+1,js:je,ke+1,comp)
                   srzy(is-1:ie+1,js:je,kc) = s(is-1:ie+1,js:je,ke+1,comp)
                else if (phys_bc(3,2) .eq. SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzy(is-1:ie+1,js:je,kc) = ZERO
                      srzy(is-1:ie+1,js:je,kc) = ZERO
                   else
                      srzy(is-1:ie+1,js:je,kc) = slzy(is-1:ie+1,js:je,kc)
                   end if
                else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel) then
                      slzy(is-1:ie+1,js:je,kc) = ZERO
                      srzy(is-1:ie+1,js:je,kc) = ZERO
                   else
                      srzy(is-1:ie+1,js:je,kc) = slzy(is-1:ie+1,js:je,kc)
                   end if
                else if (phys_bc(3,2) .eq. OUTLET) then
                   if (is_vel .and. comp .eq. 3) then
                      slzy(is-1:ie+1,js:je,kc) = max(slzy(is-1:ie+1,js:je,kc),ZERO)
                      srzy(is-1:ie+1,js:je,kc) = max(slzy(is-1:ie+1,js:je,kc),ZERO)
                   else
                      srzy(is-1:ie+1,js:je,kc) = slzy(is-1:ie+1,js:je,kc)
                   end if
                end if
             end if

             do j=js,je
                do i=is-1,ie+1
                   ! make simhzy by solving Riemann problem
                   simhzy(i,j,kc) = merge(slzy(i,j,kc),srzy(i,j,kc),wmac(i,j,k) .gt. ZERO)
                   savg = HALF*(slzy(i,j,kc)+srzy(i,j,kc))
                   simhzy(i,j,kc) = merge(simhzy(i,j,kc),savg,abs(wmac(i,j,k)) .gt. eps)
                enddo
             enddo

          endif ! end if (k .gt. ks-1)

          if (k .gt. ks) then

             !******************************************************************
             ! 9. Compute simhxz(is  :ie+1,js-1:je+1,k-1)
             !******************************************************************

             do j=js-1,je+1
                do i=is,ie+1
                   ! make slxz, srxz by updating 1D extrapolation
                   if(is_conservative(comp)) then
                      slxz(i,j,kp) = slx(i,j,kp) - (dt3/hz)*(simhz(i-1,j,kc)*wmac(i-1,j,k) - simhz(i-1,j,kp)*wmac(i-1,j,k-1))
                      srxz(i,j,kp) = srx(i,j,kp) - (dt3/hz)*(simhz(i  ,j,kc)*wmac(i  ,j,k) - simhz(i  ,j,kp)*wmac(i  ,j,k-1))
                   else
                      slxz(i,j,kp) = slx(i,j,kp) - (dt6/hz)*(wmac(i-1,j,k)+wmac(i-1,j,k-1))*(simhz(i-1,j,kc)-simhz(i-1,j,kp))
                      srxz(i,j,kp) = srx(i,j,kp) - (dt6/hz)*(wmac(i  ,j,k)+wmac(i  ,j,k-1))*(simhz(i  ,j,kc)-simhz(i  ,j,kp))
                   endif
                end do
             end do

             ! impose lo side bc's
             if (phys_bc(1,1) .eq. INLET) then
                slxz(is,js-1:je+1,kp) = s(is-1,js-1:je+1,k-1,comp)
                srxz(is,js-1:je+1,kp) = s(is-1,js-1:je+1,k-1,comp)
             else if (phys_bc(1,1) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   slxz(is,js-1:je+1,kp) = ZERO
                   srxz(is,js-1:je+1,kp) = ZERO
                else
                   slxz(is,js-1:je+1,kp) = srxz(is,js-1:je+1,kp)
                end if
             else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   slxz(is,js-1:je+1,kp) = ZERO
                   srxz(is,js-1:je+1,kp) = ZERO
                else
                   slxz(is,js-1:je+1,kp) = srxz(is,js-1:je+1,kp)
                end if
             else if (phys_bc(1,1) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   slxz(is,js-1:je+1,kp) = min(srxz(is,js-1:je+1,kp),ZERO)
                   srxz(is,js-1:je+1,kp) = min(srxz(is,js-1:je+1,kp),ZERO)
                else
                   slxz(is,js-1:je+1,kp) = srxz(is,js-1:je+1,kp)
                end if
             end if

             ! impose hi side bc's
             if (phys_bc(1,2) .eq. INLET) then
                slxz(ie+1,js-1:je+1,kp) = s(ie+1,js-1:je+1,k-1,comp)
                srxz(ie+1,js-1:je+1,kp) = s(ie+1,js-1:je+1,k-1,comp)
             else if (phys_bc(1,2) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   slxz(ie+1,js-1:je+1,kp) = ZERO
                   srxz(ie+1,js-1:je+1,kp) = ZERO
                else
                   srxz(ie+1,js-1:je+1,kp) = slxz(ie+1,js-1:je+1,kp)
                end if
             else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   slxz(ie+1,js-1:je+1,kp) = ZERO
                   srxz(ie+1,js-1:je+1,kp) = ZERO
                else
                   srxz(ie+1,js-1:je+1,kp) = slxz(ie+1,js-1:je+1,kp)
                end if
             else if (phys_bc(1,2) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   slxz(ie+1,js-1:je+1,kp) = max(slxz(ie+1,js-1:je+1,kp),ZERO)
                   srxz(ie+1,js-1:je+1,kp) = max(slxz(ie+1,js-1:je+1,kp),ZERO)
                else
                   srxz(ie+1,js-1:je+1,kp) = slxz(ie+1,js-1:je+1,kp)
                end if
             end if

             do j=js-1,je+1
                do i=is,ie+1
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
                   if(is_conservative(comp)) then
                      slyz(i,j,kp) = sly(i,j,kp) - (dt3/hz)*(simhz(i,j-1,kc)*wmac(i,j-1,k) - simhz(i,j-1,kp)*wmac(i,j-1,k-1))
                      sryz(i,j,kp) = sry(i,j,kp) - (dt3/hz)*(simhz(i,j  ,kc)*wmac(i,j  ,k) - simhz(i,j  ,kp)*wmac(i,j  ,k-1))
                   else
                      slyz(i,j,kp) = sly(i,j,kp) - (dt6/hz)*(wmac(i,j-1,k)+wmac(i,j-1,k-1))*(simhz(i,j-1,kc)-simhz(i,j-1,kp))
                      sryz(i,j,kp) = sry(i,j,kp) - (dt6/hz)*(wmac(i,j  ,k)+wmac(i,j  ,k-1))*(simhz(i,j  ,kc)-simhz(i,j  ,kp))
                   endif
                end do
             end do
             
             ! impose lo side bc's
             if (phys_bc(2,1) .eq. INLET) then
                slyz(is-1:ie+1,js,kp) = s(is-1:ie+1,js-1,k-1,comp)
                sryz(is-1:ie+1,js,kp) = s(is-1:ie+1,js-1,k-1,comp)
             else if (phys_bc(2,1) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   slyz(is-1:ie+1,js,kp) = ZERO
                   sryz(is-1:ie+1,js,kp) = ZERO
                else
                   slyz(is-1:ie+1,js,kp) = sryz(is-1:ie+1,js,kp)
                end if
             else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   slyz(is-1:ie+1,js,kp) = ZERO
                   sryz(is-1:ie+1,js,kp) = ZERO
                else
                   slyz(is-1:ie+1,js,kp) = sryz(is-1:ie+1,js,kp)
                end if
             else if (phys_bc(2,1) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 2) then
                   slyz(is-1:ie+1,js,kp) = min(sryz(is-1:ie+1,js,kp),ZERO)
                   sryz(is-1:ie+1,js,kp) = min(sryz(is-1:ie+1,js,kp),ZERO)
                else
                   slyz(is-1:ie+1,js,kp) = sryz(is-1:ie+1,js,kp)
                end if
             end if

             ! impose hi side bc's
             if (phys_bc(2,2) .eq. INLET) then
                slyz(is-1:ie+1,je+1,kp) = s(is-1:ie+1,je+1,k-1,comp)
                sryz(is-1:ie+1,je+1,kp) = s(is-1:ie+1,je+1,k-1,comp)
             else if (phys_bc(2,2) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   slyz(is-1:ie+1,je+1,kp) = ZERO
                   sryz(is-1:ie+1,je+1,kp) = ZERO
                else
                   sryz(is-1:ie+1,je+1,kp) = slyz(is-1:ie+1,je+1,kp)
                end if
             else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   slyz(is-1:ie+1,je+1,kp) = ZERO
                   sryz(is-1:ie+1,je+1,kp) = ZERO
                else
                   sryz(is-1:ie+1,je+1,kp) = slyz(is-1:ie+1,je+1,kp)
                end if
             else if (phys_bc(2,2) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 2) then
                   slyz(is-1:ie+1,je+1,kp) = max(slyz(is-1:ie+1,je+1,kp),ZERO)
                   sryz(is-1:ie+1,je+1,kp) = max(slyz(is-1:ie+1,je+1,kp),ZERO)
                else
                   sryz(is-1:ie+1,je+1,kp) = slyz(is-1:ie+1,je+1,kp)
                end if
             end if

             do j=js,je+1
                do i=is-1,ie+1
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
                   if(is_conservative(comp)) then
                      sedgelx(i,j) = slx(i,j,kp) &
                           - (dt2/hy)*(simhyz(i-1,j+1,kp)*vmac(i-1,j+1,k-1) - simhyz(i-1,j,kp)*vmac(i-1,j,k-1)) &
                           - (dt2/hz)*(simhzy(i-1,j  ,kc)*wmac(i-1,j  ,k  ) - simhzy(i-1,j,kp)*wmac(i-1,j,k-1)) &
                           + (dt2/hy)*s(i-1,j,k-1,comp)*(vmac(i-1,j+1,k-1)-vmac(i-1,j,k-1)) &
                           + (dt2/hz)*s(i-1,j,k-1,comp)*(wmac(i-1,j  ,k  )-wmac(i-1,j,k-1))
                      sedgerx(i,j) = srx(i,j,kp) &
                           - (dt2/hy)*(simhyz(i  ,j+1,kp)*vmac(i  ,j+1,k-1) - simhyz(i  ,j,kp)*vmac(i  ,j,k-1)) &
                           - (dt2/hz)*(simhzy(i  ,j  ,kc)*wmac(i  ,j  ,k  ) - simhzy(i  ,j,kp)*wmac(i  ,j,k-1)) &
                           + (dt2/hy)*s(i  ,j,k-1,comp)*(vmac(i  ,j+1,k-1)-vmac(i  ,j,k-1)) &
                           + (dt2/hz)*s(i  ,j,k-1,comp)*(wmac(i  ,j  ,k  )-wmac(i  ,j,k-1))
                   else
                      sedgelx(i,j) = slx(i,j,kp) &
                           - (dt4/hy)*(vmac(i-1,j+1,k-1)+vmac(i-1,j,k-1))*(simhyz(i-1,j+1,kp)-simhyz(i-1,j,kp)) &
                           - (dt4/hz)*(wmac(i-1,j  ,k  )+wmac(i-1,j,k-1))*(simhzy(i-1,j  ,kc)-simhzy(i-1,j,kp))
                      sedgerx(i,j) = srx(i,j,kp) &
                           - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i  ,j,k-1))*(simhyz(i  ,j+1,kp)-simhyz(i  ,j,kp)) &
                           - (dt4/hz)*(wmac(i  ,j  ,k  )+wmac(i  ,j,k-1))*(simhzy(i  ,j  ,kc)-simhzy(i  ,j,kp))
                   endif

                   ! if use_minion is true, we have already accounted for source terms
                   ! in slx and srx; otherwise, we need to account for them here.
                   if(.not. use_minion) then
                      sedgelx(i,j) = sedgelx(i,j) + dt2*force(i-1,j,k-1,comp)
                      sedgerx(i,j) = sedgerx(i,j) + dt2*force(i  ,j,k-1,comp)
                   endif

                   ! if use_minion .and. is_conservative, we have already accounted for divu
                   ! in slx and srx; otherwise, we account for it here
                   if(.not.(use_minion .and. is_conservative(comp))) then
                      sedgelx(i,j) = sedgelx(i,j) - dt2*s(i-1,j,k-1,comp)*divu(i-1,j,k-1)
                      sedgerx(i,j) = sedgerx(i,j) - dt2*s(i  ,j,k-1,comp)*divu(i  ,j,k-1)
                   endif

                   ! make sedgex by solving Riemann problem
                   ! boundary conditions enforced outside of i,j,k loop
                   sedgex(i,j,k-1,comp) = merge(sedgelx(i,j),sedgerx(i,j),umac(i,j,k-1) .gt. ZERO)
                   savg = HALF*(sedgelx(i,j)+sedgerx(i,j))
                   sedgex(i,j,k-1,comp) = merge(sedgex(i,j,k-1,comp),savg,abs(umac(i,j,k-1)) .gt. eps)
                end do
             end do
             
             ! impose lo side bc's
             if (phys_bc(1,1) .eq. INLET) then
                sedgex(is,js:je,k-1,comp) = s(is-1,js:je,k-1,comp)
             else if (phys_bc(1,1) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(is,js:je,k-1,comp) = ZERO
                else
                   sedgex(is,js:je,k-1,comp) = sedgerx(is,js:je)
                end if
             else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgex(is,js:je,k-1,comp) = ZERO
                else
                   sedgex(is,js:je,k-1,comp) = sedgerx(is,js:je)
                end if
             else if (phys_bc(1,1) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(is,js:je,k-1,comp) = min(sedgerx(is,js:je),ZERO)
                else
                   sedgex(is,js:je,k-1,comp) = sedgerx(is,js:je)
                end if
             end if

             ! impose hi side bc's
             if (phys_bc(1,2) .eq. INLET) then
                sedgex(ie+1,js:je,k-1,comp) = s(ie+1,js:je,k-1,comp)
             else if (phys_bc(1,2) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(ie+1,js:je,k-1,comp) = ZERO
                else
                   sedgex(ie+1,js:je,k-1,comp) = sedgelx(ie+1,js:je)
                end if
             else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgex(ie+1,js:je,k-1,comp) = ZERO
                else
                   sedgex(ie+1,js:je,k-1,comp) = sedgelx(ie+1,js:je)
                end if
             else if (phys_bc(1,2) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(ie+1,js:je,k-1,comp) = max(sedgelx(ie+1,js:je),ZERO)
                else
                   sedgex(ie+1,js:je,k-1,comp) = sedgelx(ie+1,js:je)
                end if
             end if

             ! create fluxes
             if (is_conservative(comp)) then
                do j=js,je
                   do i=is,ie+1
                      fluxx(i,j,k-1,comp) = sedgex(i,j,k-1,comp)*umac(i,j,k-1)
                   enddo
                enddo
             endif

             !******************************************************************
             ! 12.Compute sedgey(is  :ie  ,js  :je+1,k-1)
             !******************************************************************

             do j=js,je+1
                do i=is,ie
                   ! make sedgely, sedgery
                   if(is_conservative(comp)) then
                      sedgely(i,j) = sly(i,j,kp) &
                           - (dt2/hx)*(simhxz(i+1,j-1,kp)*umac(i+1,j-1,k-1) - simhxz(i,j-1,kp)*umac(i,j-1,k-1)) &
                           - (dt2/hz)*(simhzx(i  ,j-1,kc)*wmac(i  ,j-1,k  ) - simhzx(i,j-1,kp)*wmac(i,j-1,k-1)) &
                           + (dt2/hx)*s(i,j-1,k-1,comp)*(umac(i+1,j-1,k-1)-umac(i,j-1,k-1)) &
                           + (dt2/hz)*s(i,j-1,k-1,comp)*(wmac(i  ,j-1,k  )-wmac(i,j-1,k-1))
                      sedgery(i,j) = sry(i,j,kp) &
                           - (dt2/hx)*(simhxz(i+1,j  ,kp)*umac(i+1,j  ,k-1) - simhxz(i,j  ,kp)*umac(i,j  ,k-1)) &
                           - (dt2/hz)*(simhzx(i  ,j  ,kc)*wmac(i  ,j  ,k  ) - simhzx(i,j  ,kp)*wmac(i,j  ,k-1)) &
                           + (dt2/hx)*s(i,j  ,k-1,comp)*(umac(i+1,j  ,k-1)-umac(i,j  ,k-1)) &
                           + (dt2/hz)*s(i,j  ,k-1,comp)*(wmac(i  ,j  ,k  )-wmac(i,j  ,k-1))
                   else
                      sedgely(i,j) = sly(i,j,kp) &
                           - (dt4/hx)*(umac(i+1,j-1,k-1)+umac(i,j-1,k-1))*(simhxz(i+1,j-1,kp)-simhxz(i,j-1,kp)) &
                           - (dt4/hz)*(wmac(i  ,j-1,k  )+wmac(i,j-1,k-1))*(simhzx(i  ,j-1,kc)-simhzx(i,j-1,kp))
                      sedgery(i,j) = sry(i,j,kp) &
                           - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j  ,k-1))*(simhxz(i+1,j  ,kp)-simhxz(i,j  ,kp)) &
                           - (dt4/hz)*(wmac(i  ,j  ,k  )+wmac(i,j  ,k-1))*(simhzx(i  ,j  ,kc)-simhzx(i,j  ,kp))
                   endif

                   ! if use_minion is true, we have already accounted for source terms
                   ! in sly and sry; otherwise, we need to account for them here.
                   if(.not. use_minion) then
                      sedgely(i,j) = sedgely(i,j) + dt2*force(i,j-1,k-1,comp)
                      sedgery(i,j) = sedgery(i,j) + dt2*force(i,j  ,k-1,comp)
                   endif

                   ! if use_minion .and. is_conservative, we have already accounted for divu
                   ! in sly and sry; otherwise, we account for it here
                   if(.not.(use_minion .and. is_conservative(comp))) then
                      sedgely(i,j) = sedgely(i,j) - dt2*s(i,j-1,k-1,comp)*divu(i,j-1,k-1)
                      sedgery(i,j) = sedgery(i,j) - dt2*s(i,j  ,k-1,comp)*divu(i,j  ,k-1)
                   endif
                   ! make sedgey by solving Riemann problem
                   ! boundary conditions enforced outside of i,j,k loop
                   sedgey(i,j,k-1,comp) = merge(sedgely(i,j),sedgery(i,j),vmac(i,j,k-1) .gt. ZERO)
                   savg = HALF*(sedgely(i,j)+sedgery(i,j))
                   sedgey(i,j,k-1,comp) = merge(sedgey(i,j,k-1,comp),savg,abs(vmac(i,j,k-1)) .gt. eps)
                end do
             end do
             
             ! impose lo side bc's
             if (phys_bc(2,1) .eq. INLET) then
                sedgey(is:ie,js,k-1,comp) = s(is:ie,js-1,k-1,comp)
             else if (phys_bc(2,1) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(is:ie,js,k-1,comp) = ZERO
                else
                   sedgey(is:ie,js,k-1,comp) = sedgery(is:ie,js)
                end if
             else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgey(is:ie,js,k-1,comp) = ZERO
                else
                   sedgey(is:ie,js,k-1,comp) = sedgery(is:ie,js)
                end if
             else if (phys_bc(2,1) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(is:ie,js,k-1,comp) = min(sedgery(is:ie,js),ZERO)
                else
                   sedgey(is:ie,js,k-1,comp) = sedgery(is:ie,js)
                end if
             end if

             ! impose hi side bc's
             if (phys_bc(2,2) .eq. INLET) then
                sedgey(is:ie,je+1,k-1,comp) = s(is:ie,je+1,k-1,comp)
             else if (phys_bc(2,2) .eq. SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(is:ie,je+1,k-1,comp) = ZERO
                else
                   sedgey(is:ie,je+1,k-1,comp) = sedgely(is:ie,je+1)
                end if
             else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
                if (is_vel) then
                   sedgey(is:ie,je+1,k-1,comp) = ZERO
                else
                   sedgey(is:ie,je+1,k-1,comp) = sedgely(is:ie,je+1)
                end if
             else if (phys_bc(2,2) .eq. OUTLET) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(is:ie,je+1,k-1,comp) = max(sedgely(is:ie,je+1),ZERO)
                else
                   sedgey(is:ie,je+1,k-1,comp) = sedgely(is:ie,je+1)
                end if
             end if

             ! create fluxes
             if (is_conservative(comp)) then
                do j=js,je+1
                   do i=is,ie
                      fluxy(i,j,k-1,comp) = sedgey(i,j,k-1,comp)*vmac(i,j,k-1)
                   enddo
                enddo
             endif

          endif ! end if(k .gt. ks)

          !******************************************************************
          ! 13. Cycle indices
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

  end subroutine mkflux_3d

  subroutine mkflux_debug_3d(s,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz,umac,vmac,wmac, &
                             force,divu,lo,dx,dt,is_vel,phys_bc,adv_bc,ng,&
                             is_conservative)
    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_minion

    integer, intent(in) :: lo(:),ng

    real(kind=dp_t),intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:, :)
    real(kind=dp_t),intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(inout) ::  fluxz(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t),intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)
    real(kind=dp_t),intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3) -1:,:)
    real(kind=dp_t),intent(in   ) ::   divu(lo(1)- 1:,lo(2)- 1:,lo(3) -1:)

    real(kind=dp_t),intent(in) :: dt,dx(:)
    integer        ,intent(in) :: phys_bc(:,:)
    integer        ,intent(in) :: adv_bc(:,:,:)
    logical        ,intent(in) :: is_vel, is_conservative(:)

    ! Local variables
    real(kind=dp_t), allocatable:: slopex(:,:,:,:)
    real(kind=dp_t), allocatable:: slopey(:,:,:,:)
    real(kind=dp_t), allocatable:: slopez(:,:,:,:)

    real(kind=dp_t) hx, hy, hz, dt2, dt3, dt4, dt6, savg
    real(kind=dp_t) :: abs_eps, eps, umax

    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke
    integer :: comp,ncomp

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

    ncomp = size(s,dim=4)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,ncomp))

    do k = lo(3)-1,hi(3)+1
       call slopex_2d(s(:,:,k,:),slopex(:,:,k,:),lo,hi,ng,ncomp,adv_bc)
       call slopey_2d(s(:,:,k,:),slopey(:,:,k,:),lo,hi,ng,ncomp,adv_bc)
    end do
    call slopez_3d(s,slopez,lo,hi,ng,ncomp,adv_bc)

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

    abs_eps = 1.0d-8

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

    ! loop over components
    do comp = 1,ncomp

       !******************************************************************
       ! Create s_{\i-\half\e_x}^x, etc.
       !******************************************************************

       ! loop over appropriate x-faces
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,comp)
                srx(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i  ,j,k,comp)

                ! add source terms
                if(use_minion) then
                   slx(i,j,k) = slx(i,j,k) + dt2*force(i-1,j,k,comp)
                   srx(i,j,k) = srx(i,j,k) + dt2*force(i  ,j,k,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   slx(i,j,k) = slx(i,j,k) - dt2*s(i-1,j,k,comp)*divu(i-1,j,k)
                   srx(i,j,k) = srx(i,j,k) - dt2*s(i  ,j,k,comp)*divu(i  ,j,k)
                endif
             end do
          end do
       end do

       ! impose lo side bc's
       if (phys_bc(1,1) .eq. INLET) then
          slx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
          srx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
       else if (phys_bc(1,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slx(is,js-1:je+1,ks-1:ke+1) = ZERO
             srx(is,js-1:je+1,ks-1:ke+1) = ZERO
          else
             slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
          end if
       else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slx(is,js-1:je+1,ks-1:ke+1) = ZERO
             srx(is,js-1:je+1,ks-1:ke+1) = ZERO
          else
             slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
          end if
       else if (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slx(is,js-1:je+1,ks-1:ke+1) = min(srx(is,js-1:je+1,ks-1:ke+1),ZERO)
             srx(is,js-1:je+1,ks-1:ke+1) = min(srx(is,js-1:je+1,ks-1:ke+1),ZERO)
          else
             slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(1,2) .eq. INLET) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
          srx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
       else if (phys_bc(1,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
             srx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
          else
             srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
          end if
       else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
             srx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
          else
             srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
          end if
       else if (phys_bc(1,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slx(ie+1,js-1:je+1,ks-1:ke+1) = max(slx(ie+1,js-1:je+1,ks-1:ke+1),ZERO)
             srx(ie+1,js-1:je+1,ks-1:ke+1) = max(slx(ie+1,js-1:je+1,ks-1:ke+1),ZERO)
          else
             srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
          end if
       end if

       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
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
                sly(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,comp)
                sry(i,j,k) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j  ,k,comp)

                ! add source terms
                if(use_minion) then
                   sly(i,j,k) = sly(i,j,k) + dt2*force(i,j-1,k,comp)
                   sry(i,j,k) = sry(i,j,k) + dt2*force(i,j  ,k,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   sly(i,j,k) = sly(i,j,k) - dt2*s(i,j-1,k,comp)*divu(i,j-1,k)
                   sry(i,j,k) = sry(i,j,k) - dt2*s(i,j  ,k,comp)*divu(i,j  ,k)
                endif
             end do
          end do
       end do
       
       ! impose lo side bc's
       if (phys_bc(2,1) .eq. INLET) then
          sly(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
          sry(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
       else if (phys_bc(2,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,js,ks-1:ke+1) = ZERO
             sry(is-1:ie+1,js,ks-1:ke+1) = ZERO
          else
             sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
          end if
       else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sly(is-1:ie+1,js,ks-1:ke+1) = ZERO
             sry(is-1:ie+1,js,ks-1:ke+1) = ZERO
          else
             sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
          end if
       else if (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,js,ks-1:ke+1) = min(sry(is-1:ie+1,js,ks-1:ke+1),ZERO)
             sry(is-1:ie+1,js,ks-1:ke+1) = min(sry(is-1:ie+1,js,ks-1:ke+1),ZERO)
          else
             sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(2,2) .eq. INLET) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
          sry(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
       else if (phys_bc(2,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
             sry(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
          else
             sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
          end if
       else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             sly(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
             sry(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
          else
             sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
          end if
       else if (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             sly(is-1:ie+1,je+1,ks-1:ke+1) = max(sly(is-1:ie+1,je+1,ks-1:ke+1),ZERO)
             sry(is-1:ie+1,je+1,ks-1:ke+1) = max(sly(is-1:ie+1,je+1,ks-1:ke+1),ZERO)
          else
             sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
          end if
       end if

       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
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
                slz(i,j,k) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,comp)
                srz(i,j,k) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k  ,comp)

                ! add source terms
                if(use_minion) then
                   slz(i,j,k) = slz(i,j,k) + dt2*force(i,j,k-1,comp)
                   srz(i,j,k) = srz(i,j,k) + dt2*force(i,j,k  ,comp)
                endif

                ! add divu contribution
                if(use_minion .and. is_conservative(comp)) then
                   slz(i,j,k) = slz(i,j,k) - dt2*s(i,j,k-1,comp)*divu(i,j,k-1)
                   srz(i,j,k) = srz(i,j,k) - dt2*s(i,j,k  ,comp)*divu(i,j,k  )
                endif
             enddo
          enddo
       enddo
       
       ! impose lo side bc's
       if (phys_bc(3,1) .eq. INLET) then
          slz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
          srz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
       else if (phys_bc(3,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slz(is-1:ie+1,js-1:je+1,ks) = ZERO
             srz(is-1:ie+1,js-1:je+1,ks) = ZERO
          else
             slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
          end if
       else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slz(is-1:ie+1,js-1:je+1,ks) = ZERO
             srz(is-1:ie+1,js-1:je+1,ks) = ZERO
          else
             slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
          end if
       else if (phys_bc(3,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slz(is-1:ie+1,js-1:je+1,ks) = min(srz(is-1:ie+1,js-1:je+1,ks),ZERO)
             srz(is-1:ie+1,js-1:je+1,ks) = min(srz(is-1:ie+1,js-1:je+1,ks),ZERO)
          else
             slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(3,2) .eq. INLET) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
          srz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
       else if (phys_bc(3,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
             srz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
          else
             srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
          end if
       else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
             srz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
          else
             srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
          end if
       else if (phys_bc(3,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slz(is-1:ie+1,js-1:je+1,ke+1) = max(slz(is-1:ie+1,js-1:je+1,ke+1),ZERO)
             srz(is-1:ie+1,js-1:je+1,ke+1) = max(slz(is-1:ie+1,js-1:je+1,ke+1),ZERO)
          else
             srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
          end if
       end if

       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
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
                if(is_conservative(comp)) then
                   slxy(i,j,k) = slx(i,j,k) &
                        - (dt3/hy)*(simhy(i-1,j+1,k)*vmac(i-1,j+1,k) - simhy(i-1,j,k)*vmac(i-1,j,k))
                   srxy(i,j,k) = srx(i,j,k) &
                        - (dt3/hy)*(simhy(i  ,j+1,k)*vmac(i  ,j+1,k) - simhy(i  ,j,k)*vmac(i  ,j,k))
                else
                   slxy(i,j,k) = slx(i,j,k) &
                        - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                   srxy(i,j,k) = srx(i,j,k) &
                        - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k))*(simhy(i  ,j+1,k)-simhy(i  ,j,k))
                endif
             enddo
          enddo
       enddo
       
       ! impose lo side bc's
       if (phys_bc(1,1) .eq. INLET) then
          slxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
          srxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
       else if (phys_bc(1,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slxy(is,js:je,ks-1:ke+1) = ZERO
             srxy(is,js:je,ks-1:ke+1) = ZERO
          else
             slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
          end if
       else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slxy(is,js:je,ks-1:ke+1) = ZERO
             srxy(is,js:je,ks-1:ke+1) = ZERO
          else
             slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
          end if
       else if (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slxy(is,js:je,ks-1:ke+1) = min(srxy(is,js:je,ks-1:ke+1),ZERO)
             srxy(is,js:je,ks-1:ke+1) = min(srxy(is,js:je,ks-1:ke+1),ZERO)
          else
             slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(1,2) .eq. INLET) then
          slxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
          srxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
       else if (phys_bc(1,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slxy(ie+1,js:je,ks-1:ke+1) = ZERO
             srxy(ie+1,js:je,ks-1:ke+1) = ZERO
          else
             srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
          end if
       else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slxy(ie+1,js:je,ks-1:ke+1) = ZERO
             srxy(ie+1,js:je,ks-1:ke+1) = ZERO
          else
             srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
          end if
       else if (phys_bc(1,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slxy(ie+1,js:je,ks-1:ke+1) = max(slxy(ie+1,js:je,ks-1:ke+1),ZERO)
             srxy(ie+1,js:je,ks-1:ke+1) = max(slxy(ie+1,js:je,ks-1:ke+1),ZERO)
          else
             srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
          end if
       end if

       do k=ks-1,ke+1
          do j=js,je
             do i=is,ie+1
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
                if(is_conservative(comp)) then
                   slxz(i,j,k) = slx(i,j,k) &
                        - (dt3/hz)*(simhz(i-1,j,k+1)*wmac(i-1,j,k+1) - simhz(i-1,j,k)*wmac(i-1,j,k))
                   srxz(i,j,k) = srx(i,j,k) &
                        - (dt3/hz)*(simhz(i  ,j,k+1)*wmac(i  ,j,k+1) - simhz(i  ,j,k)*wmac(i  ,j,k))
                else
                   slxz(i,j,k) = slx(i,j,k) &
                        - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k))*(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                   srxz(i,j,k) = srx(i,j,k) &
                        - (dt6/hz)*(wmac(i  ,j,k+1)+wmac(i  ,j,k))*(simhz(i  ,j,k+1)-simhz(i  ,j,k))
                endif
             enddo
          enddo
       enddo

       ! impose lo side bc's
       if (phys_bc(1,1) .eq. INLET) then
          slxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
          srxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
       else if (phys_bc(1,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slxz(is,js-1:je+1,ks:ke) = ZERO
             srxz(is,js-1:je+1,ks:ke) = ZERO
          else
             slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
          end if
       else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slxz(is,js-1:je+1,ks:ke) = ZERO
             srxz(is,js-1:je+1,ks:ke) = ZERO
          else
             slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
          end if
       else if (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slxz(is,js-1:je+1,ks:ke) = min(srxz(is,js-1:je+1,ks:ke),ZERO)
             srxz(is,js-1:je+1,ks:ke) = min(srxz(is,js-1:je+1,ks:ke),ZERO)
          else
             slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(1,2) .eq. INLET) then
          slxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
          srxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
       else if (phys_bc(1,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             slxz(ie+1,js-1:je+1,ks:ke) = ZERO
             srxz(ie+1,js-1:je+1,ks:ke) = ZERO
          else
             srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
          end if
       else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slxz(ie+1,js-1:je+1,ks:ke) = ZERO
             srxz(ie+1,js-1:je+1,ks:ke) = ZERO
          else
             srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
          end if
       else if (phys_bc(1,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 1) then
             slxz(ie+1,js-1:je+1,ks:ke) = max(slxz(ie+1,js-1:je+1,ks:ke),ZERO)
             srxz(ie+1,js-1:je+1,ks:ke) = max(slxz(ie+1,js-1:je+1,ks:ke),ZERO)
          else
             srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
          end if
       end if

       do k=ks,ke
          do j=js-1,je+1
             do i=is,ie+1
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
                if(is_conservative(comp)) then
                   slyx(i,j,k) = sly(i,j,k) &
                        - (dt3/hx)*(simhx(i+1,j-1,k)*umac(i+1,j-1,k) - simhx(i,j-1,k)*umac(i,j-1,k))
                   sryx(i,j,k) = sry(i,j,k) &
                        - (dt3/hx)*(simhx(i+1,j  ,k)*umac(i+1,j  ,k) - simhx(i,j  ,k)*umac(i,j  ,k))
                else
                   slyx(i,j,k) = sly(i,j,k) &
                        - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                   sryx(i,j,k) = sry(i,j,k) &
                        - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k))*(simhx(i+1,j  ,k)-simhx(i,j  ,k))
                endif
             enddo
          enddo
       enddo

       ! impose lo side bc's
       if (phys_bc(2,1) .eq. INLET) then
          slyx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
          sryx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
       else if (phys_bc(2,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             slyx(is:ie,js,ks-1:ke+1) = ZERO
             sryx(is:ie,js,ks-1:ke+1) = ZERO
          else
             slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
          end if
       else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slyx(is:ie,js,ks-1:ke+1) = ZERO
             sryx(is:ie,js,ks-1:ke+1) = ZERO
          else
             slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
          end if
       else if (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             slyx(is:ie,js,ks-1:ke+1) = min(sryx(is:ie,js,ks-1:ke+1),ZERO)
             sryx(is:ie,js,ks-1:ke+1) = min(sryx(is:ie,js,ks-1:ke+1),ZERO)
          else
             slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(2,2) .eq. INLET) then
          slyx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
          sryx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
       else if (phys_bc(2,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             slyx(is:ie,je+1,ks-1:ke+1) = ZERO
             sryx(is:ie,je+1,ks-1:ke+1) = ZERO
          else
             sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
          end if
       else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slyx(is:ie,je+1,ks-1:ke+1) = ZERO
             sryx(is:ie,je+1,ks-1:ke+1) = ZERO
          else
             sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
          end if
       else if (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             slyx(is:ie,je+1,ks-1:ke+1) = max(slyx(is:ie,je+1,ks-1:ke+1),ZERO)
             sryx(is:ie,je+1,ks-1:ke+1) = max(slyx(is:ie,je+1,ks-1:ke+1),ZERO)
          else
             sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
          end if
       end if

       do k=ks-1,ke+1
          do j=js,je+1
             do i=is,ie
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
                if(is_conservative(comp)) then
                   slyz(i,j,k) = sly(i,j,k) &
                        - (dt3/hz)*(simhz(i,j-1,k+1)*wmac(i,j-1,k+1) - simhz(i,j-1,k)*wmac(i,j-1,k))
                   sryz(i,j,k) = sry(i,j,k) &
                        - (dt3/hz)*(simhz(i,j  ,k+1)*wmac(i,j  ,k+1) - simhz(i,j  ,k)*wmac(i,j  ,k))
                else
                   slyz(i,j,k) = sly(i,j,k) &
                        - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k))*(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                   sryz(i,j,k) = sry(i,j,k) &
                        - (dt6/hz)*(wmac(i,j  ,k+1)+wmac(i,j  ,k))*(simhz(i,j  ,k+1)-simhz(i,j  ,k))
                endif
             enddo
          enddo
       enddo

       ! impose lo side bc's
       if (phys_bc(2,1) .eq. INLET) then
          slyz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
          sryz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
       else if (phys_bc(2,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             slyz(is-1:ie+1,js,ks:ke) = ZERO
             sryz(is-1:ie+1,js,ks:ke) = ZERO
          else
             slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
          end if
       else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slyz(is-1:ie+1,js,ks:ke) = ZERO
             sryz(is-1:ie+1,js,ks:ke) = ZERO
          else
             slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
          end if
       else if (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             slyz(is-1:ie+1,js,ks:ke) = min(sryz(is-1:ie+1,js,ks:ke),ZERO)
             sryz(is-1:ie+1,js,ks:ke) = min(sryz(is-1:ie+1,js,ks:ke),ZERO)
          else
             slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(2,2) .eq. INLET) then
          slyz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
          sryz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
       else if (phys_bc(2,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             slyz(is-1:ie+1,je+1,ks:ke) = ZERO
             sryz(is-1:ie+1,je+1,ks:ke) = ZERO
          else
             sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
          end if
       else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slyz(is-1:ie+1,je+1,ks:ke) = ZERO
             sryz(is-1:ie+1,je+1,ks:ke) = ZERO
          else
             sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
          end if
       else if (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 2) then
             slyz(is-1:ie+1,je+1,ks:ke) = max(slyz(is-1:ie+1,je+1,ks:ke),ZERO)
             sryz(is-1:ie+1,je+1,ks:ke) = max(slyz(is-1:ie+1,je+1,ks:ke),ZERO)
          else
             sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
          end if
       end if

       do k=ks,ke
          do j=js,je+1
             do i=is-1,ie+1
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
                if(is_conservative(comp)) then
                   slzx(i,j,k) = slz(i,j,k) &
                        - (dt3/hx)*(simhx(i+1,j,k-1)*umac(i+1,j,k-1) - simhx(i,j,k-1)*umac(i,j,k-1))
                   srzx(i,j,k) = srz(i,j,k) &
                        - (dt3/hx)*(simhx(i+1,j,k  )*umac(i+1,j,k  ) - simhx(i,j,k  )*umac(i,j,k  ))
                else
                   slzx(i,j,k) = slz(i,j,k) &
                        - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                   srzx(i,j,k) = srz(i,j,k) &
                        - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhx(i+1,j,k  )-simhx(i,j,k  ))
                endif
             enddo
          enddo
       enddo
       
       ! impose lo side bc's
       if (phys_bc(3,1) .eq. INLET) then
          slzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks,comp)
          srzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks,comp)
       else if (phys_bc(3,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slzx(is:ie,js-1:je+1,ks) = ZERO
             srzx(is:ie,js-1:je+1,ks) = ZERO
          else
             slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
          end if
       else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slzx(is:ie,js-1:je+1,ks) = ZERO
             srzx(is:ie,js-1:je+1,ks) = ZERO
          else
             slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
          end if
       else if (phys_bc(3,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slzx(is:ie,js-1:je+1,ks) = min(srzx(is:ie,js-1:je+1,ks),ZERO)
             srzx(is:ie,js-1:je+1,ks) = min(srzx(is:ie,js-1:je+1,ks),ZERO)
          else
             slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(3,2) .eq. INLET) then
          slzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
          srzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
       else if (phys_bc(3,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slzx(is:ie,js-1:je+1,ke+1) = ZERO
             srzx(is:ie,js-1:je+1,ke+1) = ZERO
          else
             srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
          end if
       else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slzx(is:ie,js-1:je+1,ke+1) = ZERO
             srzx(is:ie,js-1:je+1,ke+1) = ZERO
          else
             srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
          end if
       else if (phys_bc(3,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slzx(is:ie,js-1:je+1,ke+1) = max(slzx(is:ie,js-1:je+1,ke+1),ZERO)
             srzx(is:ie,js-1:je+1,ke+1) = max(slzx(is:ie,js-1:je+1,ke+1),ZERO)
          else
             srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
          end if
       end if

       do k=ks,ke+1
          do j=js-1,je+1
             do i=is,ie
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
                if(is_conservative(comp)) then
                   slzy(i,j,k) = slz(i,j,k) &
                        - (dt3/hy)*(simhy(i,j+1,k-1)*vmac(i,j+1,k-1) - simhy(i,j,k-1)*vmac(i,j,k-1))
                   srzy(i,j,k) = srz(i,j,k) &
                        - (dt3/hy)*(simhy(i,j+1,k  )*vmac(i,j+1,k  ) - simhy(i,j,k  )*vmac(i,j,k  ))
                else
                   slzy(i,j,k) = slz(i,j,k) &
                        - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                   srzy(i,j,k) = srz(i,j,k) &
                        - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhy(i,j+1,k  )-simhy(i,j,k  ))
                endif
             enddo
          enddo
       enddo

       ! impose lo side bc's
       if (phys_bc(3,1) .eq. INLET) then
          slzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks,comp)
          srzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks,comp)
       else if (phys_bc(3,1) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slzy(is-1:ie+1,js:je,ks) = ZERO
             srzy(is-1:ie+1,js:je,ks) = ZERO
          else
             slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
          end if
       else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slzy(is-1:ie+1,js:je,ks) = ZERO
             srzy(is-1:ie+1,js:je,ks) = ZERO
          else
             slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
          end if
       else if (phys_bc(3,1) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slzy(is-1:ie+1,js:je,ks) = min(srzy(is-1:ie+1,js:je,ks),ZERO)
             srzy(is-1:ie+1,js:je,ks) = min(srzy(is-1:ie+1,js:je,ks),ZERO)
          else
             slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
          end if
       end if

       ! impose hi side bc's
       if (phys_bc(3,2) .eq. INLET) then
          slzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
          srzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
       else if (phys_bc(3,2) .eq. SLIP_WALL) then
          if (is_vel .and. comp .eq. 3) then
             slzy(is-1:ie+1,js:je,ke+1) = ZERO
             srzy(is-1:ie+1,js:je,ke+1) = ZERO
          else
             srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
          end if
       else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
          if (is_vel) then
             slzy(is-1:ie+1,js:je,ke+1) = ZERO
             srzy(is-1:ie+1,js:je,ke+1) = ZERO
          else
             srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
          end if
       else if (phys_bc(3,2) .eq. OUTLET) then
          if (is_vel .and. comp .eq. 3) then
             slzy(is-1:ie+1,js:je,ke+1) = max(slzy(is-1:ie+1,js:je,ke+1),ZERO)
             srzy(is-1:ie+1,js:je,ke+1) = max(slzy(is-1:ie+1,js:je,ke+1),ZERO)
          else
             srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
          end if
       end if

       do k=ks,ke+1
          do j=js,je
             do i=is-1,ie+1
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
                if(is_conservative(comp)) then
                   sedgelx(i,j,k) = slx(i,j,k) &
                        - (dt2/hy)*(simhyz(i-1,j+1,k  )*vmac(i-1,j+1,k  ) - simhyz(i-1,j,k)*vmac(i-1,j,k)) &
                        - (dt2/hz)*(simhzy(i-1,j  ,k+1)*wmac(i-1,j  ,k+1) - simhzy(i-1,j,k)*wmac(i-1,j,k)) &
                        + (dt2/hy)*s(i-1,j,k,comp)*(vmac(i-1,j+1,  k)-vmac(i-1,j,k)) &
                        + (dt2/hz)*s(i-1,j,k,comp)*(wmac(i-1,j  ,k+1)-wmac(i-1,j,k))
                   sedgerx(i,j,k) = srx(i,j,k) &
                        - (dt2/hy)*(simhyz(i  ,j+1,k  )*vmac(i  ,j+1,  k) - simhyz(i  ,j,k)*vmac(i  ,j,k)) &
                        - (dt2/hz)*(simhzy(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) - simhzy(i  ,j,k)*wmac(i  ,j,k)) &
                        + (dt2/hy)*s(i  ,j,k,comp)*(vmac(i  ,j+1,k  )-vmac(i  ,j,k)) &
                        + (dt2/hz)*s(i  ,j,k,comp)*(wmac(i  ,j  ,k+1)-wmac(i  ,j,k))
                else
                   sedgelx(i,j,k) = slx(i,j,k) &
                        - (dt4/hy)*(vmac(i-1,j+1,k  )+vmac(i-1,j,k))*(simhyz(i-1,j+1,k  )-simhyz(i-1,j,k)) &
                        - (dt4/hz)*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k))*(simhzy(i-1,j  ,k+1)-simhzy(i-1,j,k))
                   sedgerx(i,j,k) = srx(i,j,k) &
                        - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i  ,j,k))*(simhyz(i  ,j+1,k  )-simhyz(i  ,j,k)) &
                        - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k))*(simhzy(i  ,j  ,k+1)-simhzy(i  ,j,k))
                endif

                ! if use_minion is true, we have already accounted for source terms
                ! in slx and srx; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   sedgelx(i,j,k) = sedgelx(i,j,k) + dt2*force(i-1,j,k,comp)
                   sedgerx(i,j,k) = sedgerx(i,j,k) + dt2*force(i  ,j,k,comp)
                endif

                ! if use_minion .and. is_conservative, we have already accounted for divu
                ! in slx and srx; otherwise, we account for it here
                if(.not.(use_minion .and. is_conservative(comp))) then
                   sedgelx(i,j,k) = sedgelx(i,j,k) - dt2*s(i-1,j,k,comp)*divu(i-1,j,k)
                   sedgerx(i,j,k) = sedgerx(i,j,k) - dt2*s(i  ,j,k,comp)*divu(i  ,j,k)
                endif

                ! make sedgex by solving Riemann problem
                ! boundary conditions enforced outside of i,j,k loop
                sedgex(i,j,k,comp) = merge(sedgelx(i,j,k),sedgerx(i,j,k),umac(i,j,k) .gt. ZERO)
                savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
                sedgex(i,j,k,comp) = merge(sedgex(i,j,k,comp),savg,abs(umac(i,j,k)) .gt. eps)
             enddo
          enddo
       enddo

       ! sedgex boundary conditions
       do k=ks,ke
          do j=js,je
             ! lo side
             if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(is,j,k,comp) = ZERO
                elseif (is_vel .and. comp .ne. 1) then
                   sedgex(is,j,k,comp) = merge(ZERO,sedgerx(is,j,k),phys_bc(1,1).eq.NO_SLIP_WALL)
                else
                   sedgex(is,j,k,comp) = sedgerx(is,j,k)
                endif
             elseif (phys_bc(1,1) .eq. INLET) then
                sedgex(is,j,k,comp) = s(is-1,j,k,comp)
             elseif (phys_bc(1,1) .eq. OUTLET) then
                if (is_vel .and. comp.eq.1) then
                   sedgex(is,j,k,comp) = MIN(sedgerx(is,j,k),ZERO)
                else
                   sedgex(is,j,k,comp) = sedgerx(is,j,k)
                end if
             endif

             ! hi side
             if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   sedgex(ie+1,j,k,comp) = ZERO
                else if (is_vel .and. comp .ne. 1) then
                   sedgex(ie+1,j,k,comp) = merge(ZERO,sedgelx(ie+1,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                else 
                   sedgex(ie+1,j,k,comp) = sedgelx(ie+1,j,k)
                endif
             elseif (phys_bc(1,2) .eq. INLET) then
                sedgex(ie+1,j,k,comp) = s(ie+1,j,k,comp)
             elseif (phys_bc(1,2) .eq. OUTLET) then
                if (is_vel .and. comp.eq.1) then
                   sedgex(ie+1,j,k,comp) = MAX(sedgelx(ie+1,j,k),ZERO)
                else
                   sedgex(ie+1,j,k,comp) = sedgelx(ie+1,j,k)
                end if
             endif
          enddo
       enddo

       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! make sedgelx, sedgerx
                if(is_conservative(comp)) then
                   fluxx(i,j,k,comp) = sedgex(i,j,k,comp)*umac(i,j,k)
                endif
             enddo
          enddo
       enddo

       ! loop over appropriate y-faces
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! make sedgely, sedgery
                if(is_conservative(comp)) then
                   sedgely(i,j,k) = sly(i,j,k) &
                        - (dt2/hx)*(simhxz(i+1,j-1,k  )*umac(i+1,j-1,k  ) - simhxz(i,j-1,k)*umac(i,j-1,k)) &
                        - (dt2/hz)*(simhzx(i  ,j-1,k+1)*wmac(i  ,j-1,k+1) - simhzx(i,j-1,k)*wmac(i,j-1,k)) &
                        + (dt2/hx)*s(i,j-1,k,comp)*(umac(i+1,j-1,k  )-umac(i,j-1,k)) &
                        + (dt2/hz)*s(i,j-1,k,comp)*(wmac(i  ,j-1,k+1)-wmac(i,j-1,k))
                   sedgery(i,j,k) = sry(i,j,k) &
                        - (dt2/hx)*(simhxz(i+1,j  ,k  )*umac(i+1,j  ,k  ) - simhxz(i,j  ,k)*umac(i,j  ,k)) &
                        - (dt2/hz)*(simhzx(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) - simhzx(i,j  ,k)*wmac(i,j  ,k)) &
                        + (dt2/hx)*s(i,j  ,k,comp)*(umac(i+1,j  ,k  )-umac(i,j  ,k)) &
                        + (dt2/hz)*s(i,j  ,k,comp)*(wmac(i  ,j  ,k+1)-wmac(i,j  ,k))
                else
                   sedgely(i,j,k) = sly(i,j,k) &
                        - (dt4/hx)*(umac(i+1,j-1,k  )+umac(i,j-1,k))*(simhxz(i+1,j-1,k  )-simhxz(i,j-1,k)) &
                        - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k))*(simhzx(i  ,j-1,k+1)-simhzx(i,j-1,k))
                   sedgery(i,j,k) = sry(i,j,k) &
                        - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j  ,k))*(simhxz(i+1,j  ,k  )-simhxz(i,j  ,k)) &
                        - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k))*(simhzx(i  ,j  ,k+1)-simhzx(i,j  ,k))
                endif

                ! if use_minion is true, we have already accounted for source terms
                ! in sly and sry; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   sedgely(i,j,k) = sedgely(i,j,k) + dt2*force(i,j-1,k,comp)
                   sedgery(i,j,k) = sedgery(i,j,k) + dt2*force(i,j  ,k,comp)
                endif

                ! if use_minion .and. is_conservative, we have already accounted for divu
                ! in sly and sry; otherwise, we account for it here
                if(.not.(use_minion .and. is_conservative(comp))) then
                   sedgely(i,j,k) = sedgely(i,j,k) - dt2*s(i,j-1,k,comp)*divu(i,j-1,k)
                   sedgery(i,j,k) = sedgery(i,j,k) - dt2*s(i,j  ,k,comp)*divu(i,j  ,k)
                endif

                ! make sedgey by solving Riemann problem
                ! boundary conditions enforced outside of i,j,k loop
                sedgey(i,j,k,comp) = merge(sedgely(i,j,k),sedgery(i,j,k),vmac(i,j,k) .gt. ZERO)
                savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
                sedgey(i,j,k,comp) = merge(sedgey(i,j,k,comp),savg,abs(vmac(i,j,k)) .gt. eps)
             enddo
          enddo
       enddo

       ! sedgey boundary conditions
       do k=ks,ke
          do i=is,ie
             ! lo side
             if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(i,js,k,comp) = ZERO
                elseif (is_vel .and. comp .ne. 2) then
                   sedgey(i,js,k,comp) = merge(ZERO,sedgery(i,js,k),phys_bc(2,1).eq.NO_SLIP_WALL)
                else 
                   sedgey(i,js,k,comp) = sedgery(i,js,k)
                endif
             elseif (phys_bc(2,1) .eq. INLET) then
                sedgey(i,js,k,comp) = s(i,js-1,k,comp)
             elseif (phys_bc(2,1) .eq. OUTLET) then
                if (is_vel .and. comp.eq.2) then
                   sedgey(i,js,k,comp) = MIN(sedgery(i,js,k),ZERO)
                else
                   sedgey(i,js,k,comp) = sedgery(i,js,k)
                end if
             endif

             ! hi side
             if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   sedgey(i,je+1,k,comp) = ZERO
                elseif (is_vel .and. comp .ne. 2) then
                   sedgey(i,je+1,k,comp) = merge(ZERO,sedgely(i,je+1,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                else 
                   sedgey(i,je+1,k,comp) = sedgely(i,je+1,k)
                endif
             elseif (phys_bc(2,2) .eq. INLET) then
                sedgey(i,je+1,k,comp) = s(i,je+1,k,comp)
             elseif (phys_bc(2,2) .eq. OUTLET) then
                if (is_vel .and. comp.eq.2) then
                   sedgey(i,je+1,k,comp) = MAX(sedgely(i,je+1,k),ZERO)
                else
                   sedgey(i,je+1,k,comp) = sedgely(i,je+1,k)
                end if
             endif
          enddo
       enddo

       ! create fluxes
       if (is_conservative(comp)) then
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                   fluxy(i,j,k,comp) = sedgey(i,j,k,comp)*vmac(i,j,k)
             enddo
          enddo
       enddo
       endif

       ! loop over appropriate z-faces
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! make sedgelz, sedgerz
                if(is_conservative(comp)) then
                   sedgelz(i,j,k) = slz(i,j,k) &
                        - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) - simhxy(i,j,k-1)*umac(i,j,k-1)) &
                        - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
                        + (dt2/hx)*s(i,j,k-1,comp)*(umac(i+1,j  ,k-1)-umac(i,j,k-1)) &
                        + (dt2/hy)*s(i,j,k-1,comp)*(vmac(i  ,j+1,k-1)-vmac(i,j,k-1))
                   sedgerz(i,j,k) = srz(i,j,k) &
                        - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) - simhxy(i,j,k  )*umac(i,j,k  )) &
                        - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) - simhyx(i,j,k  )*vmac(i,j,k  )) &
                        + (dt2/hx)*s(i,j,k  ,comp)*(umac(i+1,j  ,k  )-umac(i,j,k  )) &
                        + (dt2/hy)*s(i,j,k  ,comp)*(vmac(i  ,j+1,k  )-vmac(i,j,k  ))
                else
                   sedgelz(i,j,k) = slz(i,j,k) &
                        - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1))*(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
                        - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1))*(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1))
                   sedgerz(i,j,k) = srz(i,j,k) &
                        - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  ))*(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
                        - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  ))*(simhyx(i  ,j+1,k  )-simhyx(i,j,k  ))
                endif

                ! if use_minion is true, we have already accounted for source terms
                ! in slz and srz; otherwise, we need to account for them here.
                if(.not. use_minion) then
                   sedgelz(i,j,k) = sedgelz(i,j,k) + dt2*force(i,j,k-1,comp)
                   sedgerz(i,j,k) = sedgerz(i,j,k) + dt2*force(i,j,k  ,comp)
                endif

                ! if use_minion .and. is_conservative, we have already accounted for divu
                ! in slz and srz; otherwise, we account for it here
                if(.not.(use_minion .and. is_conservative(comp))) then
                   sedgelz(i,j,k) = sedgelz(i,j,k) - dt2*s(i,j,k-1,comp)*divu(i,j,k-1)
                   sedgerz(i,j,k) = sedgerz(i,j,k) - dt2*s(i,j,k  ,comp)*divu(i,j,k  )
                endif

                ! make sedgez by solving Riemann problem
                ! boundary conditions enforced outside of i,j,k loop
                sedgez(i,j,k,comp) = merge(sedgelz(i,j,k),sedgerz(i,j,k),wmac(i,j,k) .gt. ZERO)
                savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
                sedgez(i,j,k,comp) = merge(sedgez(i,j,k,comp),savg,abs(wmac(i,j,k)) .gt. eps)
             enddo
          enddo
       enddo

       ! sedgez boundary conditions
       do j=js,je
          do i=is,ie
             ! lo side
             if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 3) then
                   sedgez(i,j,ks,comp) = ZERO
                elseif (is_vel .and. comp .ne. 3) then
                   sedgez(i,j,ks,comp) = merge(ZERO,sedgerz(i,j,ks),phys_bc(3,1).eq.NO_SLIP_WALL)
                else 
                   sedgez(i,j,ks,comp) = sedgerz(i,j,ks)
                endif
             elseif (phys_bc(3,1) .eq. INLET) then
                sedgez(i,j,ks,comp) = s(i,j,ks-1,comp)
             elseif (phys_bc(3,1) .eq. OUTLET) then
                if (is_vel .and. comp.eq.3) then
                   sedgez(i,j,ks,comp) = MIN(sedgerz(i,j,ks),ZERO)
                else
                   sedgez(i,j,ks,comp) = sedgerz(i,j,ks)
                end if
             endif

             ! hi side
             if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 3) then
                   sedgez(i,j,ke+1,comp) = ZERO
                elseif (is_vel .and. comp .ne. 3) then
                   sedgez(i,j,ke+1,comp) = merge(ZERO,sedgelz(i,j,ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
                else 
                   sedgez(i,j,ke+1,comp) = sedgelz(i,j,ke+1)
                endif
             elseif (phys_bc(3,2) .eq. INLET) then
                sedgez(i,j,ke+1,comp) = s(i,j,ke+1,comp)
             elseif (phys_bc(3,2) .eq. OUTLET) then
                if (is_vel .and. comp.eq.3) then
                   sedgez(i,j,ke+1,comp) = MAX(sedgelz(i,j,ke+1),ZERO)
                else
                   sedgez(i,j,ke+1,comp) = sedgelz(i,j,ke+1)
                end if
             endif
          enddo
       enddo

       ! create fluxes
       if (is_conservative(comp)) then
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                   fluxz(i,j,k,comp) = sedgez(i,j,k,comp)*wmac(i,j,k)
                enddo
             enddo
          enddo
       endif

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

  end subroutine mkflux_debug_3d

end module mkflux_module
