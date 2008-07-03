module scalar_advance_sdc_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use explicit_diffusive_module
  use viscous_module
  use macproject_module
  use rxns_integrator
  use probin_module, only : nscal, diff_coef, diffusion_type, stencil_order, &
                            verbose, mg_verbose, cg_verbose

  implicit none

  private

  public :: scalar_advance_sdc, AD_provisional, sdc_interpolant

contains

  ! number of terms at various times we need to save
  integer :: n_adv = 2
  integer :: n_diff = 2
  common /sdc_params/ n_adv, n_diff


  subroutine sdc_interpolant(soln,adv,diff,t,dt)
    real(dp_t)    , intent(in   )  :: soln(:)
    real(dp_t)    , intent(in   )  :: adv(:,:)
    real(dp_t)    , intent(in   )  :: diff(:,:)
    real(dp_t)    , intent(in   )  :: t,dt

    ! local
    integer             :: i

    do i = 1, nscal-1
       soln(i) = adv(0,i) + diff(0,i) + t*(diff(1,i)-diff(0,i))/dt
    end do

  end subroutine sdc_interpolant

  subroutine scalar_advance_sdc(mla,uold,sold,snew,umac, &
                            ext_scal_force,dx,dt,the_bc_tower)

    use mkflux_module
    use mkforce_module
    use update_module
    use bl_constants_module
    use reactions

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    type(multifab), allocatable :: scal_force(:), divu(:)
    type(multifab), allocatable :: sflux(:,:), sedge(:,:)
    type(multifab), allocatable :: source(:)
    logical       , allocatable :: is_conservative(:)
    logical       , allocatable :: umac_nodal_flag(:)
    ! advection term:  adv_s(time,level) = - div(su) 
    ! these values are used to interpolate int(A+D+R)
    type(multifab) , allocatable :: adv_s(:,:)
    ! diffusion at t = n:  D_s(time,level) = laplacian s
    type(multifab) , allocatable :: D_s(:,:)
    ! an intermediate diffusion term, estimate at n+1 after 
    ! solving adv + diff, not used in the interpolation
    ! given 1st index to match data tyoes with D_s for passing
    ! to integrator
    type(multifab) , allocatable :: D_dfsn(1,:)

    integer :: n_adv = 2
    integer :: n_diff = 2
    common /sdc_params/ n_adv, n_diff


    integer         :: i,j,n, nspecies
    integer         :: comp,bc_comp,nlevs,dme
    logical         :: is_vel
    real(kind=dp_t) :: diff_fac,visc_mu
    real(kind=dp_t) :: smin,smax

    nspecies = nscal -1  ! scal = density + species  
    nlevs = mla%nlevel
    dm    = mla%dim

    allocate(umac_nodal_flag(mla%dim))
    allocate(scal_force(nlevs),divu(nlevs),source(nlevs))
    allocate(sflux(nlevs,dm),sedge(nlevs,dm))
    allocate(is_conservative(nscal))
    allocate(D_s(0:n_diff-1,nlevs),adv_s(0:n_adv-1,nlevs))
    allocate(D_dfsn(1,nlevs))


    is_conservative(:) = .true.
! for using concentrations rather than mass fractions:
!    do comp = 2, nscal
!       is_conservative(2) = .false.
!    end do

    is_vel  = .false.

    do n = 1, nlevs
       call multifab_build( scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build( divu(n),mla%la(n),    1,1)
! have no ghost cells--may need to change
       call multifab_build( source(n),mla%la(n),nspecies,0)
       do j = 0, n_diff-1
          call multifab_build( D_s(j,n),mla%la(n),nspecies,0)
          call setval(D_s(j,n),zero,all=.true.)
       enddo
       call multifab_build( D_dfsn(n),mla%la(n),nspecies,0)
       do j = 0, n_adv-1
          call multifab_build( adv_s(j,n),mla%la(n),nspecies,0)
          call setval(adv_s(j,n),zero,all=.true.)
       enddo

       do i = 1,dm
         umac_nodal_flag(:) = .false.
         umac_nodal_flag(i) = .true.
         call multifab_build(sflux(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call multifab_build(sedge(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call setval(sflux(n,i),ZERO,all=.true.)
         call setval(sedge(n,i),ZERO,all=.true.)
       end do

       call setval(scal_force(n),0.0_dp_t,all=.true.)
       call setval(divu(n),0.0_dp_t,all=.true.)
       call setval(source(n),0.0_dp_t,all=.true.)
    enddo
   


   !--Compute provisional solution at t=n+1 using forward/backward euler------------

    !*********************************************
    ! Compute D(s) = laplacian(s) at time n
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(0,:),sold,comp,bc_comp,dx,the_bc_tower)
       end do
    endif

    !************************************************************
    ! Create scalar force at time n.
    ! Provided as a source term to calculate advection term adv_s
    !************************************************************

    ! Create source/forcing term scal_force
!    diff_fac = ONE
!    call mkscalforce(nlevs,scal_force,ext_scal_force,D_n,diff_fac)

    ! no source term
    do n = 1,nlevs
       call setval(scal_force(n),zero)
    end do

    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************

    call mkflux(mla,sold,uold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_tower%bc_tower_array,is_vel,is_conservative)

    !*************************************************************
    ! Create scalar force at time n+1/2
    ! Using implicit euler, is equal to the external forcing terms
    !*************************************************************
    do n = 1, nlevs
       call multifab_copy_c(scal_force(n),nscal,ext_scal_force(n),&
                            nscal,ext_scal_force(n)%ng)
    enddo

    ! another option to do the same thing
    !diff_fac = zero     !in mkscalforce: diff_fac*source
    !call mkscalforce(nlevs,scal_force,ext_scal_force,source,diff_fac) 

    !*****************************************************************
    ! Update the scalars with conservative or convective differencing.
    !*****************************************************************

    call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s(0,:),dx,dt,is_vel, &
                is_conservative,the_bc_tower%bc_tower_array)

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

    !*********************************************
    ! Compute a provisional D(s) at time n+1
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_dfsn,snew,comp,bc_comp,dx,the_bc_tower)
       end do
    endif

    !***********************
    ! Integrate reactions
    !***********************
    
    do n = 1, nlevs
       call multifab_copy_c(snew(n),nscal,sold(n),nscal,sold(n)%ng)
    enddo

    call react(mla,the_bc_tower,snew,dt,provisional,is_conserv,adv_s,D_dfsn)

    !*********************************************
    ! Compute D(s) at time n+1
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(1),snew,comp,bc_comp,dx,the_bc_tower)
       end do
    endif
    
    !*****************************
    ! Compute aprrox to int(A+D+R)
    !*****************************
    


    deallocate(is_conservative)

    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(source(n))
       call multifab_destroy(D_dfsn(1,n))
       do i = 1,dm
         call multifab_destroy(sflux(n,i))
         call multifab_destroy(sedge(n,i))
       end do
       do i = 1,n_adv
         call multifab_destroy(adv_s(i,n))
       end do
       do i = 1,n_diff
         call multifab_destroy(D_s(i,n))
       end do
    enddo

    deallocate(umac_nodal_flag)
    deallocate(scal_force,divu,sflux,sedge)

    if (diff_coef > ZERO) then
       do  comp = 2, nscal
           bc_comp = dm+comp
           visc_mu = dt*diff_coef
           call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp)
        end do
    end if

2000 format('... level ', i2,' new min/max : density           ',e17.10,2x,e17.10)
2001 format('... level ', i2,' new min/max :  tracer           ',e17.10,2x,e17.10)

  end subroutine scalar_advance

  ! this is used in the last step of calculating the provisional solution
  subroutine AD_provisional(soln,adv,diff,t)
    real(dp_t)    , intent(in   )  :: soln(:)
    real(dp_t)    , intent(in   )  :: adv(:,:)
    real(dp_t)    , intent(in   )  :: diff(:,:)
    real(dp_t)    , intent(in   )  :: t

    ! localpwd

    integer             :: i

    do i = 1, nscal-1
       soln(i) = adv(0,i) + diff(0,i)
    end do

  end subroutine AD_provisional


end module scalar_advance_module
