module scalar_advance_module

  use bl_types
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkforce_module
  use update_module
  use setbc_module
  use layout_module

  implicit none

contains

  subroutine scalar_advance(nlevs,mla,uold,sold,snew,laps,rhohalf,umac,sedge,flux, &
                            ext_scal_force,dx,time,dt,the_bc_level,diff_coef,verbose, &
                            use_godunov_debug,use_minion)

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: laps(:)
    type(multifab) , intent(inout) :: rhohalf(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(multifab) , intent(inout) :: ext_scal_force(:)

    real(kind=dp_t), intent(inout) :: dx(:,:),time,dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: diff_coef
    integer        , intent(in   ) :: verbose
    logical        , intent(in)    :: use_godunov_debug
    logical        , intent(in)    :: use_minion

    ! local variables
    type(multifab), allocatable :: scal_force(:), divu(:)
    logical, allocatable        :: is_conservative(:)

    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer ::  ep(:,:,:,:)
    real(kind=dp_t), pointer ::  fp(:,:,:,:)

    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: lapsp(:,:,:,:)
    real(kind=dp_t), pointer ::  rp(:,:,:,:) 
    real(kind=dp_t), pointer ::  dp(:,:,:,:) 
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)

    integer         :: nscal,comp,dm,d,n
    integer         :: lo(uold(1)%dim),hi(uold(1)%dim)
    integer         :: i,ng_cell,ng_rho
    logical         :: is_vel,make_divu
    real(kind=dp_t) :: diff_fac
    real(kind=dp_t) :: half_dt
    real(kind=dp_t) :: smin, smax

    nscal = ncomp(sold(1))

    allocate(scal_force(nlevs),divu(nlevs))
    allocate(is_conservative(nscal))
    is_conservative(1) = .true.
    is_conservative(2) = .false.
    half_dt = HALF * dt

    ng_cell = uold(1)%ng
    ng_rho  = rhohalf(1)%ng
    dm      = uold(1)%dim

    is_vel  = .false.

    do n = 1, nlevs
       call multifab_build(scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build(divu(n),scal_force(n)%la,1,1)

       call setval(scal_force(n),0.0_dp_t,all=.true.)
       call setval(divu(n),0.0_dp_t,all=.true.)
    enddo

    !***********************************
    ! Create scalar force at time n.
    !***********************************

    diff_fac = ONE
    call mkscalforce(nlevs,scal_force,ext_scal_force,sold,laps,dx,diff_coef,diff_fac)

    do n = 1, nlevs

       !***********************************
       ! Create the edge states of scalar using the MAC velocity 
       !***********************************
       
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle
          sop    => dataptr(sold(n), i)
          uop    => dataptr(uold(n), i)
          sepx   => dataptr(sedge(n,1), i)
          sepy   => dataptr(sedge(n,2), i)
          fluxpx => dataptr(flux(n,1), i)
          fluxpy => dataptr(flux(n,2), i)
          ump    => dataptr(umac(n,1), i)
          vmp    => dataptr(umac(n,2), i)
          fp     => dataptr(scal_force(n) , i)
          dp     => dataptr(divu(n), i)
          lo = lwb(get_box(uold(n), i))
          hi = upb(get_box(uold(n), i))
          select case (dm)
          case (2)
             if(use_godunov_debug) then
                call mkflux_debug_2d(sop(:,:,1,:), uop(:,:,1,:), &
                                     sepx(:,:,1,:), sepy(:,:,1,:), &
                                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                                     ump(:,:,1,1), vmp(:,:,1,1), &
                                     fp(:,:,1,:), dp(:,:,1,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                     ng_cell, use_minion, is_conservative)
             else
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               fp(:,:,1,:), dp(:,:,1,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                               ng_cell, use_minion, is_conservative)
             endif
          case (3)
             sepz   => dataptr(sedge(n,3), i)
             fluxpz => dataptr(flux(n,3), i)
             wmp  => dataptr(umac(n,3), i)
             if(use_godunov_debug) then
                call mkflux_debug_3d(sop(:,:,:,:), uop(:,:,:,:), &
                                     sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                     fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                                     ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                     fp(:,:,:,:), dp(:,:,:,1), &
                                     lo, dx(n,:), dt, is_vel, &
                                     the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                                     ng_cell, use_minion, is_conservative)
             else
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               fp(:,:,:,:), dp(:,:,:,1), &
                               lo, dx(n,:), dt, is_vel, &
                               the_bc_level(n)%phys_bc_level_array(i,:,:), &
                               the_bc_level(n)%adv_bc_level_array(i,:,:,dm+1:dm+nscal), &
                               ng_cell, use_minion, is_conservative)
             endif
          end select
       end do

    enddo ! do n = 1, nlevs

    ! synchronize fluxes
    do n = nlevs, 2, -1
       do comp = 1, nscal
          if(is_conservative(comp)) then
             do d = 1, dm
                call ml_edge_restriction_c(flux(n-1,d),comp,flux(n,d),comp, &
                                           mla%mba%rr(n-1,:),d,1)
             enddo
          endif
       enddo
    enddo

    !***********************************
    ! Create scalar force at time n+1/2.
    !***********************************
    
    diff_fac = HALF
    call mkscalforce(nlevs,scal_force,ext_scal_force,sold,laps,dx,diff_coef,diff_fac)

    !***********************************
    ! Update the scalars with conservative or convective differencing.
    !***********************************
    call update(nlevs,sold,umac,sedge,flux,scal_force,snew,rhohalf,dx,dt,is_vel, &
                is_conservative,the_bc_level,mla)

    if (verbose .ge. 1) then
       do n = 1, nlevs
          do comp = 1, nscal
             smin = multifab_min_c(snew(n),comp) 
             smax = multifab_max_c(snew(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2000) n,smin,smax
             else if (comp .eq. 2) then
                if (parallel_IOProcessor()) write(6,2001) n,smin,smax
             end if
          end do
       end do
    end if

    deallocate(is_conservative)

    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(divu(n))
    enddo

2000 format('... level ', i2,' new min/max : density           ',e17.10,2x,e17.10)
2001 format('... level ', i2,' new min/max :  tracer           ',e17.10,2x,e17.10)

  end subroutine scalar_advance

  subroutine modify_force(a,targ,b,src_b,c,src_c,nc,mult,all)
    integer       , intent(in   )           :: targ, src_b, src_c
    integer       , intent(in   ), optional :: nc
    logical       , intent(in   ), optional :: all
    type(multifab), intent(inout)           :: a
    type(multifab), intent(in   )           :: b
    type(multifab), intent(in   )           :: c
    real(dp_t)    , intent(in   )           :: mult

    ! local
    real(dp_t), pointer :: ap(:,:,:,:)
    real(dp_t), pointer :: bp(:,:,:,:)
    real(dp_t), pointer :: cp(:,:,:,:)
    integer             :: i
    logical             :: lall

    lall = .false.; if ( present(all) ) lall = all
    !$OMP PARALLEL DO PRIVATE(i,ap,bp)
    do i = 1, a%nboxes
       if ( multifab_remote(a,i) ) cycle
       if ( lall ) then
          ap => dataptr(a, i, targ, nc)
          bp => dataptr(b, i, src_b, nc)
          cp => dataptr(c, i, src_c, nc)
       else
          ap => dataptr(a, i, get_ibox(a, i), targ, nc)
          bp => dataptr(b, i, get_ibox(b, i), src_b, nc)
          cp => dataptr(c, i, get_ibox(c, i), src_c, nc)
       end if
       ap = ap + mult*bp*cp
    end do
    !$OMP END PARALLEL DO
  end subroutine modify_force

end module scalar_advance_module
