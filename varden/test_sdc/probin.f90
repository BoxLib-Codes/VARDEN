! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space

  implicit none

  integer,save    :: narg, farg
  integer,save    :: dim_in, nlevs, max_levs
  integer,save    :: max_step,init_iter
  integer,save    :: plot_int, chk_int, regrid_int
  integer,save    :: verbose, mg_verbose, cg_verbose
  integer,save    :: do_initial_projection
  integer,save    :: restart
  real(dp_t),save :: cflfac,init_shrink,fixed_dt
  real(dp_t),save :: visc_coef
  real(dp_t),save :: diff_coef
  real(dp_t),save :: stop_time
  real(dp_t),save :: grav
  integer,save    :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
  integer,save    :: ng_cell,ng_grow
  integer,save    :: n_cellx,n_celly,n_cellz
  integer,save    :: ref_ratio, n_error_buf
  integer,save    :: max_grid_size
  integer,save    :: slope_order
  integer,save    :: diffusion_type
  integer,save    :: nscal
  integer,save    :: stencil_order
  logical,save    :: pmask_x,pmask_y,pmask_z
  logical,save    :: pmask_xyz(MAX_SPACEDIM)
  logical,save    :: use_godunov_debug
  logical,save    :: use_minion
  real(dp_t),save :: prob_lo_x,prob_lo_y,prob_lo_z
  real(dp_t),save :: prob_hi_x,prob_hi_y,prob_hi_z
  real(dp_t),save :: max_dt_growth
  integer, save   :: boussinesq

  logical,save    :: reactions
  integer, save   :: sdc_iters
  real(dp_t),save :: k_rxn1, k_rxn2
  logical,save    :: mass_fractions
  integer, save   :: nspec
  integer, save   :: n_rxn_steps

  ! This will be allocated and defined below
  logical   , allocatable, save :: nodal(:)
  logical   , allocatable, save :: pmask(:)
  real(dp_t), allocatable, save :: prob_hi(:)
  real(dp_t), allocatable, save :: prob_lo(:)

  integer, parameter :: MAX_ALLOWED_LEVS = 10

  character(len=128), save :: fixed_grids
  character(len=128), save :: grids_file_name

  namelist /probin/ dim_in
  namelist /probin/ stop_time
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_lo_y
  namelist /probin/ prob_lo_z
  namelist /probin/ max_step
  namelist /probin/ plot_int
  namelist /probin/ chk_int
  namelist /probin/ regrid_int
  namelist /probin/ init_iter
  namelist /probin/ cflfac
  namelist /probin/ init_shrink
  namelist /probin/ fixed_dt
  namelist /probin/ visc_coef
  namelist /probin/ diff_coef
  namelist /probin/ fixed_grids
  namelist /probin/ grids_file_name
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
  namelist /probin/ pmask_xyz
  namelist /probin/ verbose
  namelist /probin/ mg_verbose
  namelist /probin/ cg_verbose
  namelist /probin/ grav
  namelist /probin/ use_godunov_debug
  namelist /probin/ use_minion
  namelist /probin/ max_levs
  namelist /probin/ n_cellx
  namelist /probin/ n_celly
  namelist /probin/ n_cellz
  namelist /probin/ ref_ratio
  namelist /probin/ n_error_buf
  namelist /probin/ max_grid_size
  namelist /probin/ boussinesq
  namelist /probin/ diffusion_type
  namelist /probin/ slope_order


  namelist /probin/ reactions
  namelist /probin/ sdc_iters
  namelist /probin/ k_rxn1
  namelist /probin/ k_rxn2
  namelist /probin/ mass_fractions
  namelist /probin/ nscal
  namelist /probin/ nspec
  namelist /probin/ n_rxn_steps

contains

  subroutine probin_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    
    integer    :: narg, farg

    character(len=128) :: fname
    character(len=128) :: probin_env

    logical :: lexist
    logical :: need_inputs

    integer :: un, ierr

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defaults

    dim_in = 2
    nscal = 2
    nspec = 1

    allocate(nodal(dim_in))
    nodal = .true.

    grav = 0.d0
    boussinesq = 0

    max_step  = 1
    stop_time = -1.d0

    ref_ratio = 2
    n_error_buf = -1
    ng_cell = 3
    ng_grow = 1
    nlevs    = -1
    max_levs = -1

    max_grid_size = 256
!    max_grid_size = 2048

    stencil_order = 2

    init_iter = 4
    plot_int  = 0
    chk_int  = 0
    regrid_int = -1
    prob_lo_x = ZERO
    prob_lo_y = ZERO
    prob_lo_z = ZERO

    verbose = 0
    mg_verbose = 0
    cg_verbose = 0

    init_shrink =  1.0
    fixed_dt    = -1.0

    do_initial_projection  = 1

    need_inputs = .true.
    fixed_grids = ''
    grids_file_name  = ''
    restart  = -1
  
    bcx_lo = SLIP_WALL
    bcy_lo = SLIP_WALL
    bcz_lo = SLIP_WALL
    bcx_hi = SLIP_WALL
    bcy_hi = SLIP_WALL
    bcz_hi = SLIP_WALL

    pmask_x = .false.
    pmask_y = .false.
    pmask_z = .false.
  
    ! 1 = Crank-Nicolson, 2 = Backward Euler
    diffusion_type = 1

    max_dt_growth = 1.1D0

    slope_order = 4

    use_godunov_debug = .false.
    use_minion = .false.

    reactions = .false.
    mass_fractions = .true.
    sdc_iters = -1
    n_rxn_steps = 100
    k_rxn1 = zero
    k_rxn2 = zero


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

    pmask_xyz = (/pmask_x, pmask_y, pmask_z/)

    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--dim_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dim_in

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

       case ('--fixed_dt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fixed_dt

       case ('--visc_coef')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_coef

       case ('--grav')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav

       case ('--max_dt_growth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_dt_growth

       case ('--use_godunov_debug')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_godunov_debug

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

       case ('--regrid_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) regrid_int

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
          read(fname, *) pmask_x
       case ('--pmask_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_y
       case ('--pmask_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_z

       case ('--verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) verbose

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--fixed_grids')
          farg = farg + 1
          call get_command_argument(farg, value = fixed_grids)

       case ('--grids_file_name')
          farg = farg + 1
          call get_command_argument(farg, value = grids_file_name)

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! If n_error_buf hasn't been set in the inputs file.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (n_error_buf < 0 .and. fixed_grids == '') then
       if (regrid_int > 0) then
          n_error_buf = regrid_int
       else
          call bl_error('Cant have n_error_buf and regrid_int both unspecified')
       end if
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Don't set regrid_int and fixed_grids
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (fixed_grids /= '') then
      if (regrid_int > 0) &
         call bl_error('Cant have fixed_grids and regrid_int > 0.')
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize prob_lo and prob_hi
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(prob_lo(dim_in))
    prob_lo(1) = prob_lo_x
    if (dim_in > 1) prob_lo(2) = prob_lo_y
    if (dim_in > 2) prob_lo(3) = prob_lo_z

    allocate(prob_hi(dim_in))
    prob_hi(1) = prob_hi_x
    if (dim_in > 1) prob_hi(2) = prob_hi_y
    if (dim_in > 2) prob_hi(3) = prob_hi_z

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialize pmask
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(pmask(dim_in))
    pmask = .FALSE.
    pmask = pmask_xyz(1:dim_in)

  end subroutine probin_init

end module probin_module
