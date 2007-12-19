module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(nlevs,mla,uold,sold,snew,laps,rhohalf,umac,sedge,sflux, &
                            ext_scal_force,dx,dt,the_bc_level,diff_coef,verbose, &
                            use_godunov_debug,use_minion)

    use mkflux_module
    use mkforce_module
    use update_module
    use bl_constants_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: laps(:)
    type(multifab) , intent(inout) :: rhohalf(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: ext_scal_force(:)
    real(kind=dp_t), intent(inout) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: diff_coef
    integer        , intent(in   ) :: verbose
    logical        , intent(in)    :: use_godunov_debug
    logical        , intent(in)    :: use_minion

    ! local variables
    type(multifab), allocatable :: scal_force(:), divu(:)
    logical       , allocatable :: is_conservative(:)

    integer         :: nscal,comp,d,n
    integer         :: lo(uold(1)%dim),hi(uold(1)%dim)
    logical         :: is_vel
    real(kind=dp_t) :: diff_fac
    real(kind=dp_t) :: smin,smax

    nscal = ncomp(sold(1))

    allocate(scal_force(nlevs),divu(nlevs))
    allocate(is_conservative(nscal))
    is_conservative(1) = .true.
    is_conservative(2) = .false.

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

    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************

    call mkflux(nlevs,sold,uold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_level,mla,is_vel,use_minion,is_conservative,use_godunov_debug)

    !***********************************
    ! Create scalar force at time n+1/2.
    !***********************************
    
    diff_fac = HALF
    call mkscalforce(nlevs,scal_force,ext_scal_force,sold,laps,dx,diff_coef,diff_fac)

    !***********************************
    ! Update the scalars with conservative or convective differencing.
    !***********************************

    call update(nlevs,sold,umac,sedge,sflux,scal_force,snew,rhohalf,dx,dt,is_vel, &
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
