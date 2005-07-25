subroutine varden()

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use estdt_module
  use vort_module
  use pre_advance_module
  use velocity_advance_module
  use scalar_advance_module
  use macproject_module
  use hgproject_module
  use ml_restriction_module
  use bc_module
  use define_bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use plotfile_module
  use checkpoint_module
  use restart_module

  implicit none

  integer        , parameter :: NPARTICLES_MAX = 40000

  integer    :: narg, farg
  integer    :: max_step,init_iter
  integer    :: plot_int, chk_int
  integer    :: verbose, mg_verbose
  integer    :: dm
  real(dp_t) :: cflfac,init_shrink
  real(dp_t) :: visc_coef
  real(dp_t) :: diff_coef
  real(dp_t) :: stop_time
  real(dp_t) :: time,dt,half_dt,dtold,dt_hold,dt_temp
  real(dp_t) :: visc_mu, pressure_inflow_val, grav
  integer    :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
  integer    :: k,istep,ng_cell,ng_grow
  integer    :: i, d, n, nlevs, n_plot_comps, n_chk_comps, nscal
  integer    :: last_plt_written, last_chk_written
  integer    :: init_step, edge_based
  integer    :: iunit
  integer    :: comp,bc_comp
  logical    :: pmask_x,pmask_y,pmask_z

  real(dp_t) :: prob_hi_x,prob_hi_y,prob_hi_z

  integer     , allocatable :: domain_phys_bc(:,:)
  logical     , allocatable :: pmask(:)
  real(dp_t)  , pointer     :: dx(:,:)
  real(dp_t)  , allocatable :: prob_hi(:)
  type(ml_layout)           :: mla
  type(box)   , allocatable :: domain_boxes(:)

  ! Cell-based quantities
  type(multifab), pointer ::       uold(:) => Null()
  type(multifab), pointer ::       sold(:) => Null()
  type(multifab), pointer ::         gp(:) => Null()
  type(multifab), pointer ::       unew(:) => Null()
  type(multifab), pointer ::       snew(:) => Null()
  type(multifab), pointer ::    rhohalf(:) => Null()
  type(multifab), pointer ::          p(:) => Null()
  type(multifab), pointer ::       vort(:) => Null()
  type(multifab), pointer ::      force(:) => Null()
  type(multifab), pointer :: scal_force(:) => Null()
  type(multifab), pointer ::   plotdata(:) => Null()
  type(multifab), pointer ::    chkdata(:) => Null()

  ! Edge-based quantities
  type(multifab), pointer ::     umac(:,:) => Null()
  type(multifab), pointer ::   utrans(:,:) => Null()
  type(multifab), pointer ::    uedge(:,:) => Null()
  type(multifab), pointer ::    sedge(:,:) => Null()

  real(kind=dp_t), pointer :: uop(:,:,:,:)
  real(kind=dp_t), pointer :: sop(:,:,:,:)
  integer,allocatable      :: lo(:),hi(:)

  character(len=128) :: fname
  character(len=128) :: probin_env
  character(len=128) :: test_set
  character(len=7) :: sd_name
  character(len=20), allocatable :: plot_names(:)
  integer :: un, ierr
  integer :: restart
  integer :: do_initial_projection
  logical :: lexist
  logical :: need_inputs
  logical :: nodal(2)
  logical :: umac_nodal_flag(2)

  type(layout)    :: la
  type(box)       :: fine_domain
  type(ml_boxarray) :: mba

  type(bc_tower) ::  the_bc_tower

  type(bc_level) ::  bc

  namelist /probin/ stop_time
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ max_step
  namelist /probin/ plot_int
  namelist /probin/ chk_int
  namelist /probin/ init_iter
  namelist /probin/ cflfac
  namelist /probin/ init_shrink
  namelist /probin/ visc_coef
  namelist /probin/ diff_coef
  namelist /probin/ test_set
  namelist /probin/ restart
  namelist /probin/ do_initial_projection
  namelist /probin/ bcx_lo
  namelist /probin/ bcx_hi
  namelist /probin/ bcy_lo
  namelist /probin/ bcy_hi
  namelist /probin/ bcz_lo
  namelist /probin/ bcz_hi
  namelist /probin/ pmask_x
  namelist /probin/ pmask_y
  namelist /probin/ pmask_z
  namelist /probin/ verbose
  namelist /probin/ mg_verbose

  nodal = .true.

  ng_cell = 3
  ng_grow = 1

  narg = command_argument_count()

! Defaults
  max_step  = 1
  init_iter = 4
  plot_int  = 0
  chk_int  = 0

  init_shrink = 1.0

  dm = 2
  nscal = 2

! grav = -9.8
  grav = 0.0

  allocate(pmask(dm))

  do_initial_projection  = 1

  need_inputs = .true.
  test_set = ''
  restart  = -1
  last_plt_written = -1
  last_chk_written = -1

  bcx_lo = SLIP_WALL
  bcy_lo = SLIP_WALL
  bcz_lo = SLIP_WALL
  bcx_hi = SLIP_WALL
  bcy_hi = SLIP_WALL
  bcz_hi = SLIP_WALL

  pmask = .FALSE.

  call get_environment_variable('PROBIN', probin_env, status = ierr)
  if ( need_inputs .AND. ierr == 0 ) then
     un = unit_new()
     open(unit=un, file = probin_env, status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .false.
  end if

  farg = 1
  if ( need_inputs .AND. narg >= 1 ) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs = .false.
     end if
  end if

  inquire(file = 'inputs_varden', exist = lexist)
  if ( need_inputs .AND. lexist ) then
     un = unit_new()
     open(unit=un, file = 'inputs_varden', status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .false.
  end if

  pmask(1) = pmask_x
  if (dm > 1) pmask(2) = pmask_y
  if (dm > 2) pmask(3) = pmask_z

  if ( .true. ) then
     do while ( farg <= narg )
        call get_command_argument(farg, value = fname)
        select case (fname)

        case ('--prob_hi_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_x

        case ('--prob_hi_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_y

        case ('--prob_hi_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_z

        case ('--cfl')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) cflfac

        case ('--init_shrink')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) init_shrink

        case ('--visc_coef')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) visc_coef

        case ('--diff_coef')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) diff_coef

        case ('--stop_time')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) stop_time

        case ('--max_step')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) max_step

        case ('--plot_int')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) plot_int

        case ('--chk_int')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) chk_int

        case ('--init_iter')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) init_iter

        case ('--bcx_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_lo
        case ('--bcy_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_lo
        case ('--bcz_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_lo
        case ('--bcx_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_hi
        case ('--bcy_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_hi
        case ('--bcz_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_hi

        case ('--pmask_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask(1)
        case ('--pmask_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask(2)
        case ('--pmask_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask(3)

        case ('--verbose')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) verbose

        case ('--mg_verbose')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) mg_verbose

        case ('--test_set')
           farg = farg + 1
           call get_command_argument(farg, value = test_set)

        case ('--restart')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) restart

        case ('--do_initial_projection')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) do_initial_projection

        case ('--')
           farg = farg + 1
           exit

        case default
           if ( .not. parallel_q() ) then
              write(*,*) 'UNKNOWN option = ', fname
              call bl_error("MAIN")
           end if
        end select

        farg = farg + 1
     end do
  end if

  allocate(prob_hi(dm))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set up plot_names for writing plot files.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(plot_names(2*dm+nscal+1))

  plot_names(1) = "x_vel"
  plot_names(2) = "y_vel"
  if (dm > 2) plot_names(3) = "z_vel"
  plot_names(dm+1) = "density"
  if (nscal > 1) plot_names(dm+2) = "tracer"
  plot_names(dm+nscal+1) = "vort"
  plot_names(dm+nscal+2) = "gpx"
  plot_names(dm+nscal+3) = "gpy"
  if (dm > 2) plot_names(dm+nscal+4) = "gpz"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the arrays and read the restart data if restart >= 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call read_a_hgproj_grid(mba, test_set)
  if (mba%dim .ne. dm) then
    print *,'BOXARRAY DIM NOT CONSISTENT WITH SPECIFIED DIM '
    stop
  end if
  nlevs = mba%nlevel
  call ml_layout_build(mla,mba,pmask)

! Define dx at the base level, then at refined levels.
  allocate(dx(nlevs,dm))

  prob_hi(1) = prob_hi_x
  if (dm > 1) prob_hi(2) = prob_hi_y
  if (dm > 2) prob_hi(3) = prob_hi_z

  do i = 1,dm
    dx(1,i) = prob_hi(i) / float(mba%pd(1)%hi(i)-mba%pd(1)%lo(i)+1)
  end do
  do n = 2,nlevs
    dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  allocate(uold(nlevs),sold(nlevs),p(nlevs),gp(nlevs))

  do n = 1,nlevs
     call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
     call multifab_build(   sold(n), mla%la(n),    dm, ng_cell)
     call multifab_build(     gp(n), mla%la(n),    dm,       1)
     call multifab_build(      p(n), mla%la(n),     1,       1, nodal)
     call setval(uold(n),0.0_dp_t, all=.true.)
     call setval(sold(n),0.0_dp_t, all=.true.)
     call setval(  gp(n),0.0_dp_t, all=.true.)
     call setval(   p(n),0.0_dp_t, all=.true.)
  end do

  if (restart >= 0) &
     call fill_restart_data(nscal,ng_cell,restart,uold,sold,gp,p,dx,time,dt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize all remaining arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(unew(nlevs),snew(nlevs))
  allocate(umac(nlevs,dm),utrans(nlevs,dm))
  allocate(uedge(nlevs,dm),sedge(nlevs,dm))
  allocate(rhohalf(nlevs),vort(nlevs))
  allocate(force(nlevs),scal_force(nlevs))

  do n = nlevs,1,-1
     call multifab_build(   unew(n), mla%la(n),    dm, ng_cell)
     call multifab_build(   snew(n), mla%la(n), nscal, ng_cell)

     call multifab_build(scal_force(n), mla%la(n), nscal, 1)
     call multifab_build(     force(n), mla%la(n),    dm, 1)
     call multifab_build(   rhohalf(n), mla%la(n),     1, 1)
     call multifab_build(      vort(n), mla%la(n),     1, 0)

     do i = 1,dm
       umac_nodal_flag = .false.
       umac_nodal_flag(i) = .true.
       call multifab_build(  umac(n,i), mla%la(n),    1, 1, nodal = umac_nodal_flag)
       call multifab_build(utrans(n,i), mla%la(n),    1, 1, nodal = umac_nodal_flag)
       call multifab_build( uedge(n,i), mla%la(n),   dm, 0, nodal = umac_nodal_flag)
       call multifab_build( sedge(n,i), mla%la(n),nscal, 0, nodal = umac_nodal_flag)
     end do

     call setval(  unew(n),ZERO, all=.true.)
     call setval(  snew(n),ZERO, all=.true.)
     do i = 1,dm
       call setval(  umac(n,i),ZERO, all=.true.)
       call setval(utrans(n,i),ZERO, all=.true.)
       call setval( uedge(n,i),ZERO, all=.true.)
       call setval( sedge(n,i),ZERO, all=.true.)
     end do
 
     call setval(      vort(n),ZERO, all=.true.)
     call setval(   rhohalf(n),ONE, all=.true.)
     call setval(     force(n),ZERO, 1,dm-1,all=.true.)
     call setval(     force(n),grav,dm,   1,all=.true.)
     call setval(scal_force(n),ZERO, all=.true.)
  end do

  la = mla%la(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate the arrays for the boundary conditions at the physical boundaries.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(domain_phys_bc(dm,2))

  allocate(domain_boxes(nlevs))
  do n = 1,nlevs
     domain_boxes(n) = layout_get_pd(mla%la(n))
  end do

! Put the bc values from the inputs file into domain_phys_bc
  domain_phys_bc(1,1) = bcx_lo
  domain_phys_bc(1,2) = bcx_hi
  if (dm > 1) then
     domain_phys_bc(2,1) = bcy_lo
     domain_phys_bc(2,2) = bcy_hi
  end if
  if (dm > 2) then
     domain_phys_bc(3,1) = bcz_lo
     domain_phys_bc(3,2) = bcz_hi
  end if

  do i = 1, dm
     if ( pmask(i) ) domain_phys_bc(i,:) = BC_PER
  end do

! Build the arrays for each grid from the domain_bc arrays.
  call bc_tower_build( the_bc_tower,mla,domain_phys_bc,domain_boxes,nscal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then
     time    = ZERO
     dt      = 1.d20
     dt_temp = ONE
     do n = 1,nlevs
        call initdata(uold(n),sold(n),dx(n,:),prob_hi, &
                      the_bc_tower%bc_tower_array(n),nscal)
     end do
     do n = nlevs,2,-1
        call ml_cc_restriction(uold(n-1),uold(n),mba%rr(n-1,:))
        call ml_cc_restriction(sold(n-1),sold(n),mba%rr(n-1,:))
     end do

!    Note that we use rhohalf, filled with 1 at this point, as a temporary
!       in order to do a constant-density initial projection.
     if (do_initial_projection > 0) then
       call hgproject(mla,uold,rhohalf,p,gp,dx,dt_temp, &
                      the_bc_tower,verbose,mg_verbose)
       do n = 1,nlevs
          call setval( p(n)  ,0.0_dp_t, all=.true.)
          call setval(gp(n)  ,0.0_dp_t, all=.true.)
       end do
     end if

     if (bcy_lo == OUTLET) then
        pressure_inflow_val = .16
        print *,'IMPOSING INFLOW PRESSURE '
        call impose_pressure_bcs(p,mla,pressure_inflow_val)
     end if
  endif

  allocate(lo(dm),hi(dm))
  do n = 1,nlevs
     bc = the_bc_tower%bc_tower_array(n)
     do i = 1, uold(n)%nboxes
       if ( multifab_remote(uold(n), i) ) cycle
       uop => dataptr(uold(n), i)
       sop => dataptr(sold(n), i)
       lo =  lwb(get_box(uold(n), i))
       hi =  upb(get_box(uold(n), i))
       select case (dm)
          case (2)
            do d = 1,dm
              call setbc_2d(uop(:,:,1,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,d), &
                            dx(n,:),d)
            end do
            do d = 1,nscal
              call setbc_2d(sop(:,:,1,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,dm+d), &
                            dx(n,:),dm+d)
            end do
         end select
     end do

     call multifab_fill_boundary(uold(n))
     call multifab_fill_boundary(sold(n))

!    This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,all=.true.)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,all=.true.)

  end do

  dtold = dt

  do n = 1,nlevs
     dt_hold = dt
     call estdt(uold(n),sold(n),gp(n),force(n),dx(n,:),cflfac,dtold,dt)
     dt = min(dt_hold,dt)
  end do
  dt = dt * init_shrink

  half_dt = HALF * dt

  allocate(plotdata(nlevs))
  n_plot_comps = 2*dm + nscal + 1
  do n = 1,nlevs
     call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
  end do

  allocate(chkdata(nlevs))
  n_chk_comps = 2*dm + nscal
  do n = 1,nlevs
     call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin the initial iterations to define an initial pressure field.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then
     if (init_iter > 0) then
       if (verbose .eq. 1) print *,'DOING ',init_iter,' INITIAL ITERATIONS ' 
       do istep = 1,init_iter

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(uold(n),uold(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
           call multifab_fill_ghost_cells(sold(n),sold(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,dm+1,nscal)
           call multifab_fill_ghost_cells(gp(n),gp(n-1),fine_domain, &
                                          ng_grow,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        do n = 1,nlevs
           call advance_premac(uold(n),sold(n),&
                               umac(n,:),uedge(n,:),utrans(n,:),& 
                               gp(n),p(n),force(n), &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               visc_coef,&
                               verbose,mg_verbose)
        end do

        call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose)

        do n = 1,nlevs
           call scalar_advance (uold(n),sold(n),snew(n),rhohalf(n),&
                                umac(n,:),sedge(n,:),utrans(n,:), &
                                scal_force(n),&
                                dx(n,:),time,dt, &
                                the_bc_tower%bc_tower_array(n), &
                                diff_coef,verbose,mg_verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(snew(n),snew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,dm+1,nscal)
        end do

        if (diff_coef > ZERO) then
          comp = 2
          bc_comp = dm+comp
          visc_mu = HALF*dt*diff_coef
          call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp,&
                                 mg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(snew(n))
        end do

        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),rhohalf(n),&
                                 umac(n,:),uedge(n,:), &
                                 utrans(n,:),gp(n),p(n), &
                                 force(n), &
                                 dx(n,:),time,dt, &
                                 the_bc_tower%bc_tower_array(n), &
                                 visc_coef,verbose,mg_verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        if (visc_coef > ZERO) then
           visc_mu = HALF*dt*visc_coef
           call visc_solve(mla,unew,rhohalf,dx,visc_mu,the_bc_tower,mg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(unew(n))
        end do

!       Project the new velocity field.
        call hgproject(mla,unew,rhohalf,p,gp,dx,dt, &
                       the_bc_tower,verbose,mg_verbose)

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        if (verbose .eq. 1) then
           do n = 1,nlevs
              print *,'MAX OF UOLD ', norm_inf(uold(n),1,1),' AT LEVEL ',n 
              print *,'MAX OF VOLD ', norm_inf(uold(n),2,1),' AT LEVEL ',n
              print *,'MAX OF UNEW ', norm_inf(unew(n),1,1),' AT LEVEL ',n 
              print *,'MAX OF VNEW ', norm_inf(unew(n),2,1),' AT LEVEL ',n
           end do
           print *,' '
        end if
       end do
     end if

     do n = 1,nlevs
        call make_vorticity(vort(n),uold(n),dx(n,:), &
                            the_bc_tower%bc_tower_array(n))
     end do

!     This writes a plotfile.
      istep = 0
      do n = 1,nlevs
         call multifab_copy_c(plotdata(n),1           ,uold(n),1,dm)
         call multifab_copy_c(plotdata(n),1+dm        ,sold(n),1,nscal)
         call multifab_copy_c(plotdata(n),1+dm+nscal  ,vort(n),1,1)
         call multifab_copy_c(plotdata(n),1+dm+nscal+1,  gp(n),1,dm)
      end do
      write(unit=sd_name,fmt='("plt",i4.4)') istep
      call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                     mba%pd(1), time, dx(1,:))
      last_plt_written = istep

!     This writes a checkpoint file.

      do n = 1,nlevs
         call multifab_copy_c(chkdata(n),1         ,uold(n),1,dm)
         call multifab_copy_c(chkdata(n),1+dm      ,sold(n),1,nscal)
         call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
      end do
      write(unit=sd_name,fmt='("chk",i4.4)') istep

      call checkpoint_write(sd_name, chkdata, p, mba%rr, dx, time, dt)
      last_chk_written = istep
 
  end if

  if (restart < 0) then
     init_step = 1
  else
     init_step = restart+1
  end if

  if (max_step >= init_step) then

     do istep = init_step,max_step

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(uold(n),uold(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
           call multifab_fill_ghost_cells(sold(n),sold(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,dm+1,nscal)
           call multifab_fill_ghost_cells(gp(n),gp(n-1),fine_domain, &
                                          ng_grow,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        do n = 1,nlevs

           call multifab_fill_boundary(uold(n))
           call multifab_fill_boundary(sold(n))
           call multifab_fill_boundary(gp(n))

!          if (verbose .eq. 1) then
!             print *,'MAX OF UOLD ', norm_inf(uold(n),1,1),' AT LEVEL ',n 
!             print *,'MAX OF VOLD ', norm_inf(uold(n),2,1),' AT LEVEL ',n
!          end if
   
           dtold = dt
           if (istep > 1) &
             call estdt(uold(n),sold(n),gp(n),force(n),dx(n,:),cflfac,dtold,dt)

        end do

        half_dt = HALF * dt

        do n = 1,nlevs

           call advance_premac(uold(n),sold(n),&
                               umac(n,:),uedge(n,:), &
                               utrans(n,:),gp(n),p(n), &
                               force(n), &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               visc_coef,verbose,mg_verbose)
        end do

        call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose)

        do n = 1,nlevs
           call scalar_advance (uold(n),sold(n),snew(n),rhohalf(n),&
                                umac(n,:),sedge(n,:),utrans(n,:),&
                                scal_force(n),&
                                dx(n,:),time,dt, &
                                the_bc_tower%bc_tower_array(n), &
                                diff_coef,verbose,mg_verbose)
        end do

        if (diff_coef > ZERO) then
          comp = 2
          bc_comp = dm+comp
          visc_mu = HALF*dt*diff_coef
          call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp, &
                                 mg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(snew(n))
        end do

        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),rhohalf(n),&
                                 umac(n,:),uedge(n,:),utrans(n,:),&
                                 gp(n),p(n),force(n), &
                                 dx(n,:),time,dt, &
                                 the_bc_tower%bc_tower_array(n), &
                                 visc_coef,verbose,mg_verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(uold(n),uold(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        if (visc_coef > ZERO) then
           visc_mu = HALF*dt*visc_coef
           call visc_solve(mla,unew,rhohalf,dx,visc_mu,the_bc_tower,mg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(unew(n))
        end do

!       Project the new velocity field.
        call hgproject(mla,unew,rhohalf,p,gp,dx,dt, &
                       the_bc_tower,verbose,mg_verbose)

        time = time + dt

        if (verbose .eq. 1) then
           do n = 1,nlevs
              print *,'MAX OF UNEW ', norm_inf(unew(n),1,1),' AT LEVEL ',n, &
                      ' AND TIME ',time
              print *,'MAX OF VNEW ', norm_inf(unew(n),2,1),' AT LEVEL ',n, &
                      ' AND TIME ',time
           end do
           print *,' '
        end if

        write(6,1000) istep,time,dt

        do n = 1,nlevs
           call multifab_copy_c(uold(n),1,unew(n),1,dm)
           call multifab_copy_c(sold(n),1,snew(n),1,nscal)
        end do

         if (plot_int > 0 .and. mod(istep,plot_int) .eq. 0) then
            do n = 1,nlevs
               call make_vorticity(vort(n),uold(n),dx(n,:), &
                                   the_bc_tower%bc_tower_array(n))
               call multifab_copy_c(plotdata(n),1           ,unew(n),1,dm)
               call multifab_copy_c(plotdata(n),1+dm        ,snew(n),1,nscal)
               call multifab_copy_c(plotdata(n),1+dm+nscal  ,vort(n),1,1)
               call multifab_copy_c(plotdata(n),1+dm+nscal+1,  gp(n),1,dm)
            end do
            write(unit=sd_name,fmt='("plt",i4.4)') istep
            call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                           mba%pd(1), time, dx(1,:))
            last_plt_written = istep
         end if

         if (chk_int > 0 .and. mod(istep,chk_int) .eq. 0) then
            do n = 1,nlevs
               call multifab_copy_c(chkdata(n),1         ,unew(n),1,dm)
               call multifab_copy_c(chkdata(n),1+dm      ,snew(n),1,nscal)
               call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
            end do
            write(unit=sd_name,fmt='("chk",i4.4)') istep
            call checkpoint_write(sd_name, chkdata, p, mba%rr, dx, time, dt)
            last_chk_written = istep
         end if

     end do

1000   format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)

       if (last_plt_written .ne. max_step) then
!       This writes a plotfile.
        do n = 1,nlevs
          call make_vorticity(vort(n),uold(n),dx(n,:), &
                              the_bc_tower%bc_tower_array(n))
          call multifab_copy_c(plotdata(n),1           ,unew(n),1,dm)
          call multifab_copy_c(plotdata(n),1+dm        ,snew(n),1,nscal)
          call multifab_copy_c(plotdata(n),1+dm+nscal  ,vort(n),1,1)
          call multifab_copy_c(plotdata(n),1+dm+nscal+1,  gp(n),1,dm)
        end do
        write(unit=sd_name,fmt='("plt",i4.4)') max_step
        call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                       mba%pd(1), time, dx(1,:))
       end if

       if (last_chk_written .ne. max_step) then
!       This writes a checkpoint file.
        do n = 1,nlevs
          call multifab_copy_c(chkdata(n),1         ,unew(n),1,dm)
          call multifab_copy_c(chkdata(n),1+dm      ,snew(n),1,nscal)
          call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
        end do
        write(unit=sd_name,fmt='("chk",i4.4)') max_step
        call checkpoint_write(sd_name, chkdata, p, mba%rr, dx, time, dt)
       end if

  end if

  do n = 1,nlevs
     call multifab_destroy(uold(n))
     call multifab_destroy(unew(n))
     call multifab_destroy(sold(n))
     call multifab_destroy(snew(n))
     do i = 1,dm
       call multifab_destroy(umac(n,i))
       call multifab_destroy(utrans(n,i))
       call multifab_destroy(uedge(n,i))
       call multifab_destroy(sedge(n,i))
     end do
     call multifab_destroy( p(n))
     call multifab_destroy(gp(n))
     call multifab_destroy(rhohalf(n))
     call multifab_destroy(vort(n))
     call multifab_destroy(force(n))
     call multifab_destroy(scal_force(n))
     call multifab_destroy(plotdata(n))
     call multifab_destroy(chkdata(n))
  end do

  deallocate(plot_names)
  deallocate(pmask)
  deallocate(prob_hi)
  deallocate(domain_phys_bc)
  deallocate(domain_boxes)
  deallocate(lo,hi)
  deallocate(plotdata)
  deallocate(chkdata)
  deallocate(uold,sold,p,gp)
  deallocate(dx)
  deallocate(unew,snew)
  deallocate(umac,utrans)
  deallocate(uedge,sedge)
  deallocate(rhohalf,vort)
  deallocate(force,scal_force)

  call bc_tower_destroy(the_bc_tower)

end subroutine varden
