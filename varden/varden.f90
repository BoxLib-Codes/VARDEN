subroutine varden()

  use BoxLib
  use omp_module
  use f2kcli
  use bl_constants_module
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
  use proj_parameters
  use ml_restriction_module
  use multifab_physbc_module
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
  use fillpatch_module
  use regrid_module
  use multifab_fill_ghost_module
  use viscous_module

  implicit none

  integer    :: narg, farg
  integer    :: max_step,init_iter
  integer    :: plot_int, chk_int, regrid_int
  integer    :: verbose, mg_verbose, cg_verbose
  integer    :: dim_in,dm
  real(dp_t) :: cflfac,init_shrink,fixed_dt
  real(dp_t) :: visc_coef
  real(dp_t) :: diff_coef
  real(dp_t) :: stop_time
  real(dp_t) :: time,dt,dtold,dt_lev,dt_temp
  real(dp_t) :: visc_mu, pressure_inflow_val,grav,nrm1,nrm2
  integer    :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
  integer    :: k,istep,ng_cell,ng_grow,n_cellx,n_celly,n_cellz
  integer    :: i, n, nlevs, ref_ratio, n_error_buf
  integer    ::  n_plot_comps, n_chk_comps, nscal
  integer    :: last_plt_written, last_chk_written
  integer    :: init_step
  integer    :: comp,bc_comp
  logical    :: pmask_x,pmask_y,pmask_z
  integer    :: press_comp, vort_comp
  logical    :: use_godunov_debug
  logical    :: use_minion

  real(dp_t) :: prob_hi_x,prob_hi_y,prob_hi_z

  integer     , allocatable :: domain_phys_bc(:,:)
  integer     , allocatable :: n_cell(:)
  logical     , allocatable :: pmask(:)
  logical     , allocatable :: nodal(:)
  real(dp_t)  , allocatable :: dx(:,:)
  real(dp_t)  , allocatable :: prob_hi(:)
  type(ml_layout)           :: mla,mla_new
  type(box)   , allocatable :: domain_box(:)

  ! Cell-based quantities
  type(multifab), allocatable ::     uold(:)
  type(multifab), allocatable ::     sold(:)
  type(multifab), allocatable ::       gp(:)
  type(multifab), allocatable ::     unew(:)
  type(multifab), allocatable ::     snew(:)
  type(multifab), allocatable ::  rhohalf(:)
  type(multifab), allocatable ::        p(:)
  type(multifab), allocatable :: lapu(:)
  type(multifab), allocatable :: ext_vel_force(:)
  type(multifab), allocatable :: laps(:)
  type(multifab), allocatable :: ext_scal_force(:)
  type(multifab), allocatable :: phi(:)
  type(multifab), allocatable :: Lphi(:)
  type(multifab), allocatable :: alpha(:)
  type(multifab), allocatable :: beta(:) 
  type(multifab), allocatable :: plotdata(:)
  type(multifab), pointer     ::  chkdata(:)
  type(multifab), pointer     ::    chk_p(:)
  type(multifab), pointer     :: p_temp(:)

  ! Regridding quantities
  type(multifab), allocatable ::   uold_rg(:)
  type(multifab), allocatable ::   sold_rg(:)
  type(multifab), allocatable ::     gp_rg(:)
  type(multifab), allocatable ::      p_rg(:)

  ! Edge-based quantities
  type(multifab), allocatable ::   umac(:,:)
  type(multifab), allocatable ::  uedge(:,:)
  type(multifab), allocatable ::  sedge(:,:)
  type(multifab), allocatable ::  uflux(:,:)
  type(multifab), allocatable ::  sflux(:,:)

  real(kind=dp_t), pointer :: uop(:,:,:,:)
  real(kind=dp_t), pointer :: sop(:,:,:,:)
  integer,allocatable      :: lo(:),hi(:)

  character(len=128) :: fname
  character(len=128) :: probin_env
  character(len=128) :: fixed_grids
  character(len=7) :: sd_name
  character(len=20), allocatable :: plot_names(:)
  integer :: un, ierr
  integer :: restart
  integer :: do_initial_projection
  integer :: stencil_order
  logical :: lexist
  logical :: need_inputs

  type(ml_layout) :: mla_temp
  type(layout)    :: la
  type(ml_boxarray) :: mba
  type(boxarray)  :: ba
  type(box), allocatable :: bxs(:)
  integer, allocatable  :: rr(:,:)
  integer  :: nl, max_levs
  logical  :: new_grid

  type(bc_tower) ::  the_bc_tower
  type(bc_level) ::  bc

  namelist /probin/ dim_in
  namelist /probin/ stop_time
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
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
  namelist /probin/ cg_verbose
  namelist /probin/ grav
  namelist /probin/ use_godunov_debug
  namelist /probin/ use_minion
  namelist /probin/ nlevs
  namelist /probin/ n_cellx
  namelist /probin/ n_celly
  namelist /probin/ n_cellz
  namelist /probin/ ref_ratio
  namelist /probin/ n_error_buf

  ref_ratio = 2
  n_error_buf = 2
  ng_cell = 3
  ng_grow = 1
  nlevs = 1

  stencil_order = 2

  narg = command_argument_count()

  ! Defaults
  max_step  = 1
  init_iter = 4
  plot_int  = 0
  regrid_int  = 0
  chk_int  = 0

  init_shrink = 1.0
  fixed_dt = -1.0
  nscal = 2

  grav = 0.0d0

  use_godunov_debug = .false.
  use_minion = .false.

  do_initial_projection  = 1

  need_inputs = .true.
  fixed_grids = ''
  restart  = -1
  last_plt_written = -1
  last_chk_written = -1

  bcx_lo = SLIP_WALL
  bcy_lo = SLIP_WALL
  bcz_lo = SLIP_WALL
  bcx_hi = SLIP_WALL
  bcy_hi = SLIP_WALL
  bcz_hi = SLIP_WALL

  stop_time = -1.d0

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

  call read_command_line()

  dm = dim_in
  press_comp = dm + nscal + 1

  allocate(pmask(dm))
  pmask = .FALSE.

  pmask(1) = pmask_x
  if (dm > 1) pmask(2) = pmask_y
  if (dm > 2) pmask(3) = pmask_z

  allocate(nodal(dm))
  nodal = .true.

  allocate(rr(nlevs-1,dm))
  do n = 1, nlevs-1
     rr(n,:) = ref_ratio
  enddo

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

  ! Restart from a checkpoint file
  if (restart >= 0) then
     n_chk_comps = 2*dm + nscal
     allocate(chkdata(nlevs),chk_p(nlevs))
     do n = 1,nlevs
        call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
        call multifab_build(  chk_p(n), mla%la(n), n_chk_comps, 0)
     end do
     call fill_restart_data(restart,mba,chkdata,chk_p,time,dt)
     call ml_layout_build(mla,mba,pmask)
     nlevs = mba%nlevel
     allocate(uold(nlevs),sold(nlevs),gp(nlevs),p(nlevs))
     do n = 1,nlevs
        call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
        call multifab_build(   sold(n), mla%la(n), nscal, ng_cell)
        call multifab_build(     gp(n), mla%la(n),    dm,  1)
        call multifab_build(      p(n), mla%la(n),     1,  1, nodal)
     end do
     do n = 1,nlevs
        call multifab_copy_c(uold(n),1,chkdata(n),1         ,dm)
        call multifab_copy_c(sold(n),1,chkdata(n),1+dm      ,nscal)
        call multifab_copy_c(  gp(n),1,chkdata(n),1+dm+nscal,dm)
        call multifab_copy_c(   p(n),1,  chk_p(n),1         ,1)
        call multifab_destroy(chkdata(n))
        call multifab_destroy(chk_p(n))
     end do
     deallocate(chkdata,chk_p)
! Read the grid info from the file indicated, fixed grid option
  else if (fixed_grids /= '') then
     call read_a_hgproj_grid(mba, fixed_grids)
     call ml_layout_build(mla,mba,pmask)
     nlevs = mla%nlevel
     allocate(uold(nlevs),sold(nlevs),p(nlevs),gp(nlevs))
     do n = 1,nlevs
        call make_new_state(mla%la(n),uold(n),sold(n),gp(n),p(n))
     enddo

! Adaptive gridding
  else 
     max_levs = nlevs

     ! set up hi & lo to carry indexing info
     allocate(lo(dm),hi(dm))
     allocate(n_cell(dm))
     n_cell(1) = n_cellx
     lo(:) = 0
     hi(1) = n_cellx-1
     if (dm > 1) then   
        n_cell(2) = n_celly
        hi(2) = n_celly - 1        
        if (dm > 2)  then
           n_cell(3) = n_cellz
           hi(3) = n_cellz -1
        endif
     endif

     ! make a multi-level boxarray with nlevs, make one single box over 
     ! entire domain at level 1, make that into a boxarray -> 
     ! multilevel boxarray -> ml layout.  use layout to initialize level 1 of
     ! multifabs
     call ml_boxarray_build_n(mba,max_levs,dm)
     do n = 1, max_levs-1
        mba%rr(n,:) = rr(n,:)
     enddo

     allocate(bxs(max_levs))
     allocate(uold(max_levs),sold(max_levs),p(max_levs),gp(max_levs))

       ! Build the level 1 boxarray
     call box_build_2(bxs(1),lo,hi)
     call boxarray_build_bx(mba%bas(1),bxs(1))

! Ask Ann: ok to build to max_levs or not?
     do n = 2, max_levs
        call box_build_2(bxs(n),lo,lo)
        call boxarray_build_bx(mba%bas(n),bxs(n))
     enddo

     ! build pd(:)
     mba%pd(1) = bxs(1)
     do n = 2, nlevs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     enddo

     ! Build the level 1 layout, has correct la, bas, pd, pmask
     ! higher levels are empty
     call ml_layout_build(mla,mba,pmask)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate state and temp variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nlevs = max_levs here

  allocate(ext_vel_force(nlevs),ext_scal_force(nlevs))
  allocate(lapu(nlevs),laps(nlevs))
  allocate(phi(nlevs),Lphi(nlevs))
  allocate(alpha(nlevs),beta(nlevs))
  allocate(uold_rg(nlevs),sold_rg(nlevs),p_rg(nlevs),gp_rg(nlevs))

  allocate(unew(nlevs),snew(nlevs))
  allocate(umac(nlevs,dm))
  allocate(uedge(nlevs,dm),sedge(nlevs,dm))
  allocate(uflux(nlevs,dm),sflux(nlevs,dm))
  allocate(rhohalf(nlevs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize dx, prob_hi, lo, hi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Define dx at the base level, then at refined levels.  nlevs = max_levs
  allocate(dx(nlevs,dm))

  allocate(prob_hi(dm))
  prob_hi(1) = prob_hi_x
  if (dm > 1) prob_hi(2) = prob_hi_y
  if (dm > 2) prob_hi(3) = prob_hi_z

  do i = 1,dm
     dx(1,i) = prob_hi(i) / float(mba%pd(1)%hi(i)-mba%pd(1)%lo(i)+1)
  end do
  do n = 2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate the arrays for the boundary conditions at the physical boundaries.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(domain_phys_bc(dm,2))

  allocate(domain_box(nlevs))
  do n = 1,nlevs
     domain_box(n) = layout_get_pd(mla%la(n))
  end do

  ! Put the bc values from the inputs file into domain_phys_bc
  domain_phys_bc(1,1) = bcx_lo
  domain_phys_bc(1,2) = bcx_hi
  if (pmask(1)) then
     if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
          call bl_error('MUST HAVE BCX = -1 if PMASK = T')
  end if
  if (dm > 1) then
     domain_phys_bc(2,1) = bcy_lo
     domain_phys_bc(2,2) = bcy_hi
     if (pmask(2)) then
        if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
             call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
     end if
  end if
  if (dm > 2) then
     domain_phys_bc(3,1) = bcz_lo
     domain_phys_bc(3,2) = bcz_hi
     if (pmask(3)) then
        if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
             call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
     end if
  end if

  do i = 1, dm
     if ( pmask(i) ) domain_phys_bc(i,:) = BC_PER
  end do

  ! Build the arrays for each grid from the domain_bc arrays.
  call bc_tower_build( the_bc_tower,mla,domain_phys_bc,domain_box,nscal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     time    = ZERO
     dt      = 1.d20
     dt_temp = ONE
     if (regrid_int > 0) then
        ! Build the level 1 data
        call make_new_state(mla%la(1),uold(1),sold(1),gp(1),p(1)) 

        ! Initialize the level 1 data
        call init_one_lev(1,uold(1),sold(1),dx(1,:),&
             prob_hi,the_bc_tower%bc_tower_array(1),nscal,mla%la(1))

     else 
        call initdata(nlevs,uold,sold,dx,prob_hi,the_bc_tower%bc_tower_array,nscal,mla)

     end if
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Regrid before starting the calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0 .and. nlevs > 1 .and. regrid_int > 0) then

! these all used the data in the old mla that got destroyed, need to rebuild

     do nl = 1, max_levs-1
        ! could use n_error_buf (if i wanted to set it) instead of regrid_int 
        call make_new_grids(mla,mla_new,sold(nl),dx(nl,1),regrid_int,rr,&
             nl,new_grid)     
        
        if (new_grid) then
            do n = 1,nl
               call destroy(sold(n))
               call destroy(uold(n))
               call destroy(gp(n))
               call destroy(p(n))
            enddo
            
            call destroy(mla)
            call ml_layout_build(mla, mla_new%mba, pmask)
            call destroy(mla_new)

            do n = 1,nl+1
               call make_new_state(mla%la(n),uold(n),sold(n),gp(n),p(n))
            enddo

! destroy bc_tower from level 1 initialization, make new bc_tower from mla_new
! Build the arrays for each grid from the domain_bc arrays.
! without this remake "multigrid solve: failed to converge in max_iters"
            call bc_tower_destroy(the_bc_tower)
            call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_box,nscal)
            
! fills the physical region of each level with problem data (blob now)
            call initdata(nl+1,uold,sold,dx,prob_hi,the_bc_tower%bc_tower_array,&
                 nscal,mla)

! the proper thing to do here would be to avg the fine grid onto the down
! onto the coarse grid.  not worrying about this right now

            nlevs = nl+1
     
         else 
            if (nl .eq. 1) then

               call destroy(sold(n))
               call destroy(uold(n))
               call destroy(gp(n))
               call destroy(p(n))
                           
               call destroy(mla)
               call ml_layout_build(mla, mla_new%mba, pmask)
               call destroy(mla_new)
               
               call delete_state(uold,sold,gp,p)
               call make_new_state(mla%la(1),uold(1),sold(1),gp(1),p(1))
               
               call bc_tower_destroy(the_bc_tower)
               call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_box,nscal)
               
               call initdata(1,uold,sold,dx,prob_hi,the_bc_tower%bc_tower_array,&
                    nscal,mla)
               
! the proper thing to do here would be to avg the fine grid onto the down
! onto the coarse grid.  not worrying about this right now

               nlevs = 1
               goto 2000
            else
               call destroy(mla_new)
               goto 2000
            endif
         endif
      enddo
          
! check for proper nesting
2000 write(*,*)'checking nesting'
      if(.not. ml_boxarray_properly_nested(mla%mba)) then
   write(*,*) 'not properly nested'
         do n = nlevs, 3, -1
            call boxarray_build_copy(ba, mla%mba%bas(n))
            call boxarray_coarsen(ba, mla%mba%rr(n-1,:))
            call boxarray_diff(ba, mla%mba%bas(n-1))
            call boxarray_intersection(ba, mla%mba%pd(n-1))
            if ( .not. empty(ba) ) then
               call boxarray_destroy(ba)
! buffer the cells
               
            end if
         enddo
      endif
      

! not sure we need to fill the bdry and apply phys bc inside the loop
 do n = 1,nlevs
         call multifab_fill_boundary(uold(n))
         call multifab_fill_boundary(sold(n))
         
         bc = the_bc_tower%bc_tower_array(n)
         call multifab_physbc(uold(n),1,1,   dm,   bc)
         call multifab_physbc(sold(n),1,dm+1,nscal,bc)
      end do

      call build(mla_temp,mla%mba, pmask)
      do n = 1,nlevs
         call make_new_state(mla_temp%la(n),uold_rg(n),sold_rg(n),gp_rg(n),p_rg(n))
      enddo
      

      write(*,*)'mla ready to go'
      call print(mla)
   end if

! initialize rhohalf, etc.
   call make_temps(mla)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Do the initial projection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     ! Note that we use rhohalf, filled with 1 at this point, as a temporary
     ! in order to do a constant-density initial projection.
     if (do_initial_projection > 0) then
        call hgproject(initial_projection,mla,uold,uold,rhohalf,p,gp,dx,dt_temp, &
                       the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp)
         do n = 1,nlevs
           call setval( p(n)  ,0.0_dp_t, all=.true.)
           call setval(gp(n)  ,0.0_dp_t, all=.true.)
        end do
     end if

     if (bcy_lo == OUTLET) then
        pressure_inflow_val = .16
        call impose_pressure_bcs(p,mla,pressure_inflow_val)
     end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set bc's...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1,nlevs

     call multifab_fill_boundary(uold(n))
     call multifab_fill_boundary(sold(n))

     bc = the_bc_tower%bc_tower_array(n)
     call multifab_physbc(uold(n),1,1,   dm,   bc)
     call multifab_physbc(sold(n),1,dm+1,nscal,bc)

     !    This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,ng=unew(n)%ng)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,ng=snew(n)%ng)

  end do

  dtold = dt
  dt = 1.d20
  do n = 1,nlevs
     call estdt(n,uold(n),sold(n),gp(n),ext_vel_force(n),dx(n,:), &
                cflfac,dtold,dt_lev,verbose)
     dt = min(dt,dt_lev)
  end do
  if (restart < 0) dt = dt * init_shrink
  if (fixed_dt > 0.d0) dt = fixed_dt
  if (stop_time >= 0.d0) then
     if (time+dt > stop_time) dt = min(dt, stop_time - time)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the initial iterations to define an initial pressure field.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     if (init_iter > 0) call initial_iters()

     istep = 0

     if ( plot_int > 0 ) then
        if ( mod(istep,plot_int) .eq. 0 ) then
           call write_plotfile(istep)
           last_plt_written = istep
        end if
     end if

     if ( chk_int > 0 ) then
        if ( mod(istep,chk_int) .eq. 0 ) then
           call write_checkfile(istep)
           last_chk_written = istep
        end if
     end if

     init_step = 1

  else 

     init_step = restart+1

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the real integration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0)) then

     do istep = init_step, max_step

        if ( verbose > 0 ) then
           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''
        end if

        if (nlevs > 1 .and. regrid_int > 0) then

           if (mod(istep-1,regrid_int) .eq. 0) then
              call delete_temps()

! only give sold_rg data, since this is the only thing needed for regridding
              call multifab_copy_c(sold_rg(1),1,sold(1),1)              

              do nl = 1, (max_levs-1)

! assume same step in all spatial directions
! make level n+1 grid from sold
                 call make_new_grids(mla_temp,mla_new,sold_rg(nl),dx(nl,1),&
                      regrid_int,rr,nl,new_grid)

                 if (new_grid) then
                    call destroy(mla_temp)
                    call ml_layout_build(mla_temp, mla_new%mba, pmask)
                    call destroy(mla_new)

! Build the arrays for each grid from the domain_bc arrays.
                    call bc_tower_destroy(the_bc_tower)
                    call bc_tower_build(the_bc_tower,mla_temp,domain_phys_bc,&
                         domain_box,nscal)

                    call delete_state(uold_rg,sold_rg,gp_rg,p_rg)
! again, really only need sold for this.  not sure it's a good idea to carry
! the extra multifabs through this

                    do n = 1,nl+1
                       call make_new_state(mla_temp%la(n),uold_rg(n),&
                            sold_rg(n),gp_rg(n),p_rg(n))
                    enddo

! really only need to make and fillpatch sold, since that's the condition 
! we use to refine
                    do n = 2, nl+1
                       call fillpatch(uold_rg(n),uold(n-1), &
                            ng_cell,mla_temp%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1), &
                            the_bc_tower%bc_tower_array(n  ), &
                            1,1,1,dm)
                       call fillpatch(sold_rg(n),sold(n-1), &
                            ng_cell,mla_temp%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1), &
                            the_bc_tower%bc_tower_array(n  ), &
                            1,1,dm+1,nscal)
                       call fillpatch(gp_rg(n),gp(n-1), &
                            ng_grow,mla_temp%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1), &
                            the_bc_tower%bc_tower_array(n  ), &
                            1,1,1,dm)
                    end do
                    
                    do n = 1,nl
                       call multifab_copy_c(uold_rg(n),1,uold(n),1,dm   )
                       call multifab_copy_c(sold_rg(n),1,sold(n),1,nscal)
                       call multifab_copy_c(gp_rg(n)  ,1,gp(n),  1,dm   )
                       call multifab_copy_c(p_rg(n)   ,1,p(n),   1,1    )
                    end do
                   
                    if (mla%nlevel .gt. nl) then
                       call multifab_copy_c(uold_rg(nl+1),1,uold(nl+1),1,dm   )
                       call multifab_copy_c(sold_rg(nl+1),1,sold(nl+1),1,nscal)
                       call multifab_copy_c(gp_rg(nl+1)  ,1,gp(nl+1),  1,dm   )
                       call multifab_copy_c(p_rg(nl+1)   ,1,p(nl+1),   1,1    )
                    endif

                    nlevs = nl+1
                 else 
                    if (nl == 1) then
                       call destroy(mla_temp)
                       call ml_layout_build(mla_temp, mla_new%mba, pmask)
                       call destroy(mla_new)

                       call bc_tower_destroy(the_bc_tower)
                       call bc_tower_build(the_bc_tower,mla_temp,domain_phys_bc,&
                            domain_box,nscal)

                       call delete_state(uold_rg,sold_rg,gp_rg,p_rg)
                       call make_new_state(mla_temp%la(1),uold_rg(1),&
                            sold_rg(1),gp_rg(1),p_rg(1))
                    
                       call multifab_copy_c(uold_rg(1),1,uold(1),1,dm   )
                       call multifab_copy_c(sold_rg(1),1,sold(1),1,nscal)
                       call multifab_copy_c(gp_rg(1)  ,1,gp(1),  1,dm   )
                       call multifab_copy_c(p_rg(1)   ,1,p(1),   1,1    )
                 
                       nlevs = 1
                       goto 2001
                    else
                       call destroy(mla_new)
                       goto 2001
                    endif
                 endif
              enddo
             
! check for proper nesting

2001          call destroy(mla)
              call build(mla,mla_temp%mba,pmask)

              call delete_state(uold,sold,gp,p)

              do n = 1,nlevs
                 call make_new_state(mla%la(n),uold(n),sold(n),gp(n),p(n))
              enddo

              do nl = 1,nlevs
                 call multifab_copy_c(uold(nl),1,uold_rg(nl),1,dm   )
                 call multifab_copy_c(sold(nl),1,sold_rg(nl),1,nscal)
                 call multifab_copy_c(gp(nl)  ,1,gp_rg(nl),  1,dm   )
                 call multifab_copy_c(p(nl)   ,1,p_rg(nl),   1,1    )
              end do

              do n = 1,nlevs
                 call multifab_fill_boundary(uold(n))
                 call multifab_fill_boundary(sold(n))
                 
                 bc = the_bc_tower%bc_tower_array(n)
                 call multifab_physbc(uold(n),1,1,   dm,   bc)
                 call multifab_physbc(sold(n),1,dm+1,nscal,bc)
              end do
                   
  write(*,*)'mla ready to go'
  call print(mla)

              call make_temps(mla)
!              call delete_state(uold_rg,sold_rg,gp_rg,p_rg)
           end if !  end if mod(istep-1,regrid_int) .eq. 0)
        end if  ! end if (nlevs > 1 .and. regrid_int > 0)

        do n = 2, nlevs
           call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,1,dm)
           call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,dm+1,nscal)
           call multifab_fill_ghost_cells(gp(n),gp(n-1), &
                                          ng_grow,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,1,dm)
        end do

        do n = 1,nlevs
           call multifab_fill_boundary(uold(n))
           call multifab_fill_boundary(sold(n))
           call multifab_fill_boundary(gp(n))
        end do

        if (istep > 1) then
           dtold = dt
           dt = 1.d20
           do n = 1,nlevs
              call estdt(n,uold(n),sold(n),gp(n),ext_vel_force(n),dx(n,:), &
                   cflfac,dtold,dt_lev,verbose)
              dt = min(dt,dt_lev)
           end do
           if (fixed_dt > 0.d0) dt = fixed_dt
           if (stop_time >= 0.d0) then
              if (time+dt > stop_time) dt = stop_time - time
           end if
        end if

        if ( verbose .ge. 1 ) then
           do n = 1,nlevs
              nrm1 = norm_inf(uold(n),1,1)
              nrm2 = norm_inf(uold(n),2,1)
              if ( parallel_IOProcessor() ) write(6,1001) n,time,nrm1,nrm2
           end do
           if ( parallel_IOProcessor() ) print *,' '
        end if

        ! compute Lapu
        if(visc_coef .gt. ZERO) then
           do comp = 1, dm
              do n = 1, nlevs
                 call multifab_copy_c(phi(n),1,uold(n),comp,1,1)
              enddo
              call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,comp, &
                               stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
              do n = 1, nlevs
                 call multifab_copy_c(Lapu(n),comp,Lphi(n),1)
              enddo
           enddo
        endif

        ! compute Laps for passive scalar only
        if(diff_coef .gt. ZERO) then
           do n = 1, nlevs
              call multifab_copy_c(phi(n),1,sold(n),2,1,1)
           enddo
           call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,dm+2, &
                            stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
           do n = 1, nlevs
              call multifab_copy_c(Laps(n),2,Lphi(n),1)
           enddo
        endif

        call advance_premac(nlevs,uold,sold,lapu,umac,gp,ext_vel_force,dx,dt, &
                            the_bc_tower%bc_tower_array,visc_coef,use_godunov_debug, &
                            use_minion,mla)

        call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose, &
                        press_comp)

        call scalar_advance(nlevs,mla,uold,sold,snew,laps,rhohalf,umac,sedge,sflux, &
                            ext_scal_force,dx,dt,the_bc_tower%bc_tower_array, &
                            diff_coef,verbose,use_godunov_debug,use_minion)
        
        if (diff_coef > ZERO) then
           comp = 2
           bc_comp = dm+comp
           visc_mu = HALF*dt*diff_coef
           call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp, &
                                  mg_verbose,cg_verbose,verbose)
        end if

        call velocity_advance(nlevs,mla,uold,unew,sold,lapu,rhohalf,umac,uedge,uflux,gp,p, &
                              ext_vel_force,dx,dt,the_bc_tower%bc_tower_array, &
                              visc_coef,verbose,use_godunov_debug,use_minion)

        if (visc_coef > ZERO) then
           visc_mu = HALF*dt*visc_coef
           call visc_solve(mla,unew,rhohalf,dx,visc_mu,the_bc_tower,mg_verbose, &
                           cg_verbose,verbose)
        end if

        ! Project the new velocity field.
        call hgproject(regular_timestep,mla,unew,uold,rhohalf,p,gp,dx,dt, &
                       the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp)

        time = time + dt

        if ( verbose .ge. 1 ) then
           do n = 1,nlevs
              nrm1 = norm_inf(unew(n),1,1)
              nrm2 = norm_inf(unew(n),2,1)
              if ( parallel_IOProcessor() ) write(6,1002) n,time,nrm1,nrm2
           end do
           if ( parallel_IOProcessor() ) print *,' '
        end if

        write(6,1000) istep,time,dt

        do n = 1,nlevs
           call multifab_copy_c(uold(n),1,unew(n),1,dm)
           call multifab_copy_c(sold(n),1,snew(n),1,nscal)
        end do

        if ( plot_int > 0 ) then
           if ( mod(istep,plot_int) .eq. 0 ) then
              call write_plotfile(istep)
              last_plt_written = istep
           end if
        end if

        if ( chk_int > 0 ) then
           if ( mod(istep,chk_int) .eq. 0 ) then
              call write_checkfile(istep)
              last_chk_written = istep
           end if
        end if

        call print_and_reset_fab_byte_spread()

        if (stop_time >= 0.d0 .and. time >= stop_time) goto 999

     end do ! end do istep = init_step,max_step

999  continue
     if (istep > max_step) istep = max_step

1000 format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)
1001 format('LEVEL: ',i3,' TIME: ',f14.8,' UOLD/VOLD MAX: ',e15.9,1x,e15.9)
1002 format('LEVEL: ',i3,' TIME: ',f14.8,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)

     if (last_plt_written .ne. istep .and. plot_int > 0) call write_plotfile(istep)
     if (last_chk_written .ne. istep .and. chk_int  > 0) call write_checkfile(istep)

  end if

! since these states are made from _rg's, only need to delete those states
!  call delete_state(uold,sold,gp,p)
  call delete_state(uold_rg,sold_rg,gp_rg,p_rg)
  call delete_state(uold,sold,gp,p)
  call delete_temps()

  call bc_tower_destroy(the_bc_tower)

! mla was made from mla_new, as long as mla_new gets detroyed, mla is taken
! care of
  call destroy(mla)
  call destroy(mla_temp)
  call destroy(mba)
  deallocate(bxs)

 print *, 'MEMORY STATS AT END OF RUN '
 print*, ' '
 call print(multifab_mem_stats(),    "    multifab")
 call print(fab_mem_stats(),         "         fab")
 call print(boxarray_mem_stats(),    "    boxarray")
 call print(layout_mem_stats(),      "      layout")
 call print(boxassoc_mem_stats(),    "    boxassoc")
 call print(fgassoc_mem_stats(),     "     fgassoc")
 call print(syncassoc_mem_stats(),   "   syncassoc")
 call print(copyassoc_mem_stats(),   "   copyassoc")
 call print(fluxassoc_mem_stats(),   "   fluxassoc")
 print*, ''


contains

  subroutine make_new_state(la_loc,uold_loc,sold_loc,gp_loc,p_loc)

    type(layout),intent(in   ) :: la_loc
    type(multifab ),intent(inout) :: uold_loc,sold_loc,gp_loc,p_loc


    call multifab_build(   uold_loc, la_loc,    dm, ng_cell)
    call multifab_build(   sold_loc, la_loc, nscal, ng_cell)
    call multifab_build(     gp_loc, la_loc,    dm, ng_grow)
    call multifab_build(      p_loc, la_loc,     1, ng_grow, nodal)

    call setval(  uold_loc,ZERO, all=.true.)
    call setval(  sold_loc,ZERO, all=.true.)
    call setval(    gp_loc,ZERO, all=.true.)
    call setval(     p_loc,ZERO, all=.true.)

  end subroutine make_new_state

  subroutine make_temps(mla_loc)

    type(ml_layout),intent(in   ) :: mla_loc

    ! Local variables
    logical, allocatable :: umac_nodal_flag(:)

    allocate(umac_nodal_flag(mla_loc%dim))

    do n = nlevs,1,-1
       call multifab_build(   unew(n), mla_loc%la(n),    dm, ng_cell)
       call multifab_build(   snew(n), mla_loc%la(n), nscal, ng_cell)
       call multifab_build(rhohalf(n), mla_loc%la(n),     1, 1)
       call multifab_build(ext_vel_force(n),  mla_loc%la(n),    dm, 1)
       call multifab_build(ext_scal_force(n), mla_loc%la(n), nscal, 1)
       call multifab_build(lapu(n), mla_loc%la(n),    dm, 0)
       call multifab_build(laps(n), mla_loc%la(n), nscal, 0)
       call multifab_build(phi(n), mla_loc%la(n),  1, 1)
       call multifab_build(Lphi(n), mla_loc%la(n),  1, 0)
       call multifab_build(alpha(n), mla_loc%la(n),  1, 1)
       call multifab_build(beta(n) , mla_loc%la(n), dm, 1)

       call setval(   unew(n),ZERO, all=.true.)
       call setval(   snew(n),ZERO, all=.true.)
       call setval(rhohalf(n),ONE, all=.true.)
       call setval(ext_vel_force(n) ,ZERO, 1,dm-1,all=.true.)
       call setval(ext_vel_force(n) ,grav,dm,   1,all=.true.)
       call setval(ext_scal_force(n),ZERO, all=.true.)
       call setval(lapu(n), ZERO, all=.true.)
       call setval(laps(n), ZERO, all=.true.)
       call setval(phi(n),  ZERO, all=.true.)
       call setval(Lphi(n),  ZERO, all=.true.)
       call setval(alpha(n),  ZERO, all=.true.)
       call setval(beta(n) , 1.0d0, all=.true.)

       do i = 1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(i) = .true.
          call multifab_build(  umac(n,i), mla_loc%la(n),    1, 1, nodal = umac_nodal_flag)
          call multifab_build( uedge(n,i), mla_loc%la(n),   dm, 0, nodal = umac_nodal_flag)
          call multifab_build( sedge(n,i), mla_loc%la(n),nscal, 0, nodal = umac_nodal_flag)
          call multifab_build( uflux(n,i), mla_loc%la(n),   dm, 0, nodal = umac_nodal_flag)
          call multifab_build( sflux(n,i), mla_loc%la(n),nscal, 0, nodal = umac_nodal_flag)

          call setval(  umac(n,i),ZERO, all=.true.)
          call setval( uedge(n,i),ZERO, all=.true.)
          call setval( sedge(n,i),ZERO, all=.true.)
          call setval( uflux(n,i),ZERO, all=.true.)
          call setval( sflux(n,i),ZERO, all=.true.)
       end do

    end do

  end subroutine make_temps

  subroutine delete_temps()

    do n = 1,nlevs
       call multifab_destroy(unew(n))
       call multifab_destroy(snew(n))
       call multifab_destroy(ext_vel_force(n))
       call multifab_destroy(ext_scal_force(n))
       call multifab_destroy(lapu(n))
       call multifab_destroy(laps(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(Lphi(n))
       call multifab_destroy(alpha(n))
       call multifab_destroy(beta(n))
       do i = 1,dm
          call multifab_destroy(umac(n,i))
          call multifab_destroy(uedge(n,i))
          call multifab_destroy(sedge(n,i))
          call multifab_destroy(uflux(n,i))
          call multifab_destroy(sflux(n,i))
       end do
       call multifab_destroy(rhohalf(n))
    end do

  end subroutine delete_temps

  subroutine delete_state(u,s,gp,p)

    type(multifab), intent(inout) :: u(:),s(:),p(:),gp(:)

!    do n = 1,size(u)
    do n = 1,nlevs
       call multifab_destroy( u(n))
       call multifab_destroy( s(n))
       call multifab_destroy( p(n))
       call multifab_destroy(gp(n))
    end do

  end subroutine delete_state

  subroutine initial_iters()

    if (parallel_IOProcessor() .and. verbose .ge. 1) &
         print *,'DOING ',init_iter,' INITIAL ITERATIONS ' 

    do istep = 1,init_iter

       do n = 2, nlevs
          call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                         ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm)
          call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                         ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,dm+1,nscal)
          call multifab_fill_ghost_cells(gp(n),gp(n-1), &
                                         ng_grow,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm)
       end do

       ! compute Lapu
       if(visc_coef .gt. ZERO) then
          do comp = 1, dm
             do n = 1, nlevs
                call multifab_copy_c(phi(n),1,uold(n),comp,1,1)
             enddo
             call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,comp, &
                              stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
             do n = 1, nlevs
                call multifab_copy_c(Lapu(n),comp,Lphi(n),1)
             enddo
          enddo
       endif

       ! compute Laps for passive scalar only
       if(diff_coef .gt. ZERO) then
          do n = 1, nlevs
             call multifab_copy_c(phi(n),1,sold(n),2,1,1)
          enddo
          call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,dm+2, &
                           stencil_order,mla%mba%rr,mg_verbose,cg_verbose)
          do n = 1, nlevs
             call multifab_copy_c(Laps(n),2,Lphi(n),1)
          enddo
       endif

       call advance_premac(nlevs,uold,sold,lapu,umac,gp,ext_vel_force,dx,dt, &
                           the_bc_tower%bc_tower_array,visc_coef,use_godunov_debug, &
                           use_minion,mla)

       call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose, &
                       press_comp)

       call scalar_advance(nlevs,mla,uold,sold,snew,laps,rhohalf,umac,sedge,sflux, &
                           ext_scal_force,dx,dt,the_bc_tower%bc_tower_array, &
                           diff_coef,verbose,use_godunov_debug,use_minion)

       if (diff_coef > ZERO) then
          comp = 2
          bc_comp = dm+comp
          visc_mu = HALF*dt*diff_coef
          call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp, &
                                 mg_verbose,cg_verbose,verbose)
       end if

       call velocity_advance(nlevs,mla,uold,unew,sold,lapu,rhohalf,umac,uedge,uflux,gp,p, &
                             ext_vel_force,dx,dt,the_bc_tower%bc_tower_array, &
                             visc_coef,verbose,use_godunov_debug,use_minion)

       if (visc_coef > ZERO) then
          visc_mu = HALF*dt*visc_coef
          call visc_solve(mla,unew,rhohalf,dx,visc_mu,the_bc_tower,mg_verbose, &
                          cg_verbose,verbose)
       end if

       ! Project the new velocity field.
       call hgproject(pressure_iters,mla,unew,uold,rhohalf,p,gp,dx,dt, &
                      the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp)

       if ( verbose .ge. 1 ) then
          do n = 1,nlevs
             nrm1 = norm_inf(unew(n),1,1)
             nrm2 = norm_inf(unew(n),2,1)
             if ( parallel_IOProcessor() ) write(6,1003) n,istep,nrm1,nrm2
          end do
          if ( parallel_IOProcessor() ) print *,' '
       end if

    end do

1003 format('LEVEL: ',i3,' ITER: ',i3,' UNEW/VNEW MAX: ',e15.9,1x,e15.9)

  end subroutine initial_iters

  subroutine write_plotfile(istep_to_write)

    integer, intent(in   ) :: istep_to_write
    integer                :: n,n_plot_comps

    allocate(plotdata(nlevs))
    n_plot_comps = 2*dm + nscal + 1

    do n = 1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy_c(plotdata(n),1           ,uold(n),1,dm)
       call multifab_copy_c(plotdata(n),1+dm        ,sold(n),1,nscal)

       vort_comp = 1+dm+nscal
       call make_vorticity(plotdata(n),vort_comp,uold(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       call multifab_copy_c(plotdata(n),1+dm+nscal+1,  gp(n),1,dm)
    end do
    write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
         mba%pd(1), time, dx(1,:))

    do n = 1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)

  end subroutine write_plotfile

  subroutine write_checkfile(istep_to_write)

    integer, intent(in   ) :: istep_to_write
    type(multifab), allocatable :: p_temp(:)

    allocate(chkdata(nlevs))
    allocate(p_temp(nlevs))
    n_chk_comps = 2*dm + nscal
    do n = 1,nlevs
       call multifab_build(p_temp(n),  mla%la(n), 1, ng_grow, nodal)
       call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
       call multifab_copy_c(chkdata(n),1         ,uold(n),1,dm)
       call multifab_copy_c(chkdata(n),1+dm      ,sold(n),1,nscal)
       call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
       call multifab_copy_c(p_temp(n), 1, p(n),1,1)
    end do
    write(unit=sd_name,fmt='("chk",i4.4)') istep_to_write

! uses size(p) which is max_levs.  need to fix
    call checkpoint_write(sd_name, chkdata, p_temp, mla%mba%rr, dx, time, dt, verbose)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
       call multifab_destroy(p_temp(n))
    end do
    deallocate(chkdata)
    deallocate(p_temp)

  end subroutine write_checkfile

  subroutine read_command_line()

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

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--fixed_grids')
          farg = farg + 1
          call get_command_argument(farg, value = fixed_grids)

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

  end subroutine read_command_line

end subroutine varden
