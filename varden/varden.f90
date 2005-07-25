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
  use fillpatch_module

  implicit none

  integer    :: narg, farg
  integer    :: max_step,init_iter
  integer    :: plot_int, chk_int, regrid_int
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
  real(dp_t)  , allocatable :: dx(:,:)
  real(dp_t)  , allocatable :: prob_hi(:)
  type(ml_layout)           :: mla,mla_new
  type(box)   , allocatable :: domain_boxes(:)

  ! Cell-based quantities
  type(multifab), allocatable ::     uold(:)
  type(multifab), allocatable ::     sold(:)
  type(multifab), allocatable ::       gp(:)
  type(multifab), allocatable ::     unew(:)
  type(multifab), allocatable ::     snew(:)
  type(multifab), allocatable ::  rhohalf(:)
  type(multifab), allocatable ::        p(:)
  type(multifab), allocatable ::     vort(:)
  type(multifab), allocatable ::    force(:)
  type(multifab), allocatable ::   sforce(:)
  type(multifab), allocatable :: plotdata(:)
  type(multifab), allocatable ::  chkdata(:)

  type(multifab), allocatable ::   uold_rg(:)
  type(multifab), allocatable ::   sold_rg(:)
  type(multifab), allocatable ::     gp_rg(:)
  type(multifab), allocatable ::      p_rg(:)
  type(multifab), allocatable ::  force_rg(:)
  type(multifab), allocatable :: sforce_rg(:)

  ! Edge-based quantities
  type(multifab), allocatable ::     umac(:,:)
  type(multifab), allocatable ::   utrans(:,:)
  type(multifab), allocatable ::    uedge(:,:)
  type(multifab), allocatable ::    sedge(:,:)

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

  type(layout)    :: la
  type(box)       :: fine_domain
  type(ml_boxarray) :: mba,mba_new

  type(bc_tower) ::  the_bc_tower

  type(bc_level) ::  bc

  namelist /probin/ stop_time
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ max_step
  namelist /probin/ plot_int
  namelist /probin/ regrid_int
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
  regrid_int  = 0
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

  call read_command_line()

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

  allocate(uold(nlevs),sold(nlevs),p(nlevs),gp(nlevs))
  allocate(force(nlevs),sforce(nlevs))
  allocate(uold_rg(nlevs),sold_rg(nlevs),p_rg(nlevs),gp_rg(nlevs))
  allocate(force_rg(nlevs),sforce_rg(nlevs))

  allocate(unew(nlevs),snew(nlevs))
  allocate(umac(nlevs,dm),utrans(nlevs,dm))
  allocate(uedge(nlevs,dm),sedge(nlevs,dm))
  allocate(rhohalf(nlevs),vort(nlevs))

  allocate(plotdata(nlevs),chkdata(nlevs))

  call make_new_state(mla,uold,sold,gp,p,force,sforce)
  call make_temps(mla)

  if (restart >= 0) &
     call fill_restart_data(nscal,ng_cell,restart,uold,sold,gp,p,time,dt)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin the initial iterations to define an initial pressure field.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     if (init_iter > 0) call initial_iters()

     do n = 1,nlevs
        call make_vorticity(vort(n),uold(n),dx(n,:), &
                            the_bc_tower%bc_tower_array(n))
     end do

     istep = 0

     call write_plotfile(istep)
     last_plt_written = istep

     call write_checkfile(istep)
     last_chk_written = istep

     init_step = 1
 
  else 

     init_step = restart+1

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin the real integration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (max_step >= init_step) then

     do istep = init_step,max_step

        if (regrid_int > 0 .and. mod(istep,regrid_int) .eq. 0) then
          call delete_temps()
          call make_new_grids()

          call make_new_state(mla_new,uold_rg,sold_rg,gp_rg,p_rg,force_rg,sforce_rg)

          do n = 1,nlevs
            call multifab_copy_c(  uold_rg(n),1,  uold(n),1,   dm)
            call multifab_copy_c(  sold_rg(n),1,  sold(n),1,nscal)
            call multifab_copy_c(    gp_rg(n),1,    gp(n),1,   dm)
            call multifab_copy_c(     p_rg(n),1,     p(n),1,    1)
            call multifab_copy_c( force_rg(n),1, force(n),1,   dm)
            call multifab_copy_c(sforce_rg(n),1,sforce(n),1,nscal)
          end do

          do n = 2, nlevs
             fine_domain = layout_get_pd(mla_new%la(n))
             call fillpatch(uold_rg(n),uold_rg(n-1),fine_domain, &
                            ng_cell,mla_new%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                            1,1,dm)
             call fillpatch(sold_rg(n),sold_rg(n-1),fine_domain, &
                            ng_cell,mla_new%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                            1,dm+1,nscal)
             call fillpatch(gp_rg(n),gp_rg(n-1),fine_domain, &
                            ng_grow,mla_new%mba%rr(n-1,:), &
                            the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                            1,1,dm)
          end do

          call delete_state(uold,sold,gp,p,force,sforce)

            uold = uold_rg
            sold = sold_rg
              gp =   gp_rg
               p =    p_rg
           force = force_rg
          sforce = sforce_rg

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
          end do

          call make_temps(mla_new)

          ! Build the arrays for each grid from the domain_bc arrays.
          call bc_tower_build( the_bc_tower,mla_new,domain_phys_bc,domain_boxes,nscal)

          call destroy(mla)
          mla = mla_new

        end if

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

           if (verbose .eq. 1) then
              print *,'MAX OF UOLD ', norm_inf(uold(n),1,1),' AT LEVEL ',n 
              print *,'MAX OF VOLD ', norm_inf(uold(n),2,1),' AT LEVEL ',n
              print *,'MAX OF SOLD ', norm_inf(sold(n),1,1),' AT LEVEL ',n
           end if
   
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
                                sforce(n),&
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
            call write_plotfile(istep)
            last_plt_written = istep
         end if

         if (chk_int > 0 .and. mod(istep,chk_int) .eq. 0) then
            call write_checkfile(istep)
            last_chk_written = istep
         end if

     end do

1000   format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)

       if (last_plt_written .ne. max_step) call write_plotfile(max_step)
       if (last_chk_written .ne. max_step) call write_checkfile(max_step)

  end if

  call delete_state(uold   ,sold   ,gp   ,p   ,force   ,sforce   )

  call delete_temps()

  call bc_tower_destroy(the_bc_tower)

  contains

    subroutine make_new_grids()

       type(boxarray) :: ba_loc
       type(box     ) :: bx_loc
       integer        :: lo_bx(2),hi_bx(2)

       if (istep .eq. 2) then
         test_set = "gr0_2d_2"
       else if (istep .eq. 4) then
         test_set = "gr0_2d_4"
       else if (istep .eq. 6) then
         test_set = "gr0_2d_6"
       end if
       call read_a_hgproj_grid(mba_new, test_set)

!      lo_bx = 16
!      hi_bx = 32
!      bx_loc = make_box(lo_bx,hi_bx)
!      call build(ba_loc,bx_loc)
!      ba_loc = mba_new%bas(2)
!      call ml_layout_rebuild(mla_new,mla,ba_loc)

       call ml_layout_build(mla_new,mba_new)

    end subroutine make_new_grids

    subroutine make_new_state(mla_loc,uold_loc,sold_loc,gp_loc,p_loc,force_loc,sforce_loc)
  
     type(ml_layout),intent(in   ) :: mla_loc
     type(multifab ),intent(inout) :: uold_loc(:),sold_loc(:),gp_loc(:),p_loc(:)
     type(multifab ),intent(inout) :: force_loc(:),sforce_loc(:)

     do n = 1,nlevs
        call multifab_build(   uold_loc(n), mla_loc%la(n),    dm, ng_cell)
        call multifab_build(   sold_loc(n), mla_loc%la(n), nscal, ng_cell)
        call multifab_build(     gp_loc(n), mla_loc%la(n),    dm,       1)
        call multifab_build(      p_loc(n), mla_loc%la(n),     1,       1, nodal)
        call multifab_build( sforce_loc(n), mla_loc%la(n), nscal, 1)
        call multifab_build(  force_loc(n), mla_loc%la(n),    dm, 1)

        call setval(  uold_loc(n),ZERO, all=.true.)
        call setval(  sold_loc(n),ZERO, all=.true.)
        call setval(    gp_loc(n),ZERO, all=.true.)
        call setval(     p_loc(n),ZERO, all=.true.)
        call setval( force_loc(n),ZERO, 1,dm-1,all=.true.)
        call setval( force_loc(n),grav,dm,   1,all=.true.)
        call setval(sforce_loc(n),ZERO, all=.true.)
     end do

    end subroutine make_new_state

    subroutine make_temps(mla_loc)
  
     type(ml_layout),intent(in   ) :: mla_loc
!    type(multifab), intent(inout) :: uold(:),sold(:),gp(:),p(:)

     ! Local variables
     logical :: umac_nodal_flag(2)

     do n = nlevs,1,-1
        call multifab_build(   unew(n), mla_loc%la(n),    dm, ng_cell)
        call multifab_build(   snew(n), mla_loc%la(n), nscal, ng_cell)
        call multifab_build(rhohalf(n), mla_loc%la(n),     1, 1)
        call multifab_build(   vort(n), mla_loc%la(n),     1, 0)

        call setval(  unew(n),ZERO, all=.true.)
        call setval(  snew(n),ZERO, all=.true.)
        call setval(   vort(n),ZERO, all=.true.)
        call setval(rhohalf(n),ONE, all=.true.)
   
        do i = 1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(i) = .true.
          call multifab_build(  umac(n,i), mla_loc%la(n),    1, 1, nodal = umac_nodal_flag)
          call multifab_build(utrans(n,i), mla_loc%la(n),    1, 1, nodal = umac_nodal_flag)
          call multifab_build( uedge(n,i), mla_loc%la(n),   dm, 0, nodal = umac_nodal_flag)
          call multifab_build( sedge(n,i), mla_loc%la(n),nscal, 0, nodal = umac_nodal_flag)

          call setval(  umac(n,i),ZERO, all=.true.)
          call setval(utrans(n,i),ZERO, all=.true.)
          call setval( uedge(n,i),ZERO, all=.true.)
          call setval( sedge(n,i),ZERO, all=.true.)
        end do

     end do

     n_plot_comps = 2*dm + nscal + 1
     do n = 1,nlevs
        call multifab_build(plotdata(n), mla_loc%la(n), n_plot_comps, 0)
     end do

     n_chk_comps = 2*dm + nscal
     do n = 1,nlevs
        call multifab_build(chkdata(n), mla_loc%la(n), n_chk_comps, 0)
     end do

    end subroutine make_temps

    subroutine delete_temps()

      do n = 1,nlevs
         call multifab_destroy(unew(n))
         call multifab_destroy(snew(n))
         do i = 1,dm
           call multifab_destroy(umac(n,i))
           call multifab_destroy(utrans(n,i))
           call multifab_destroy(uedge(n,i))
           call multifab_destroy(sedge(n,i))
         end do
         call multifab_destroy(rhohalf(n))
         call multifab_destroy(vort(n))
         call multifab_destroy(plotdata(n))
         call multifab_destroy(chkdata(n))
      end do

    end subroutine delete_temps

    subroutine delete_state(u,s,gp,p,f,sf)

      type(multifab), intent(inout) :: u(:),s(:),p(:),gp(:),f(:),sf(:)

      do n = 1,size(u)
         call multifab_destroy( u(n))
         call multifab_destroy( s(n))
         call multifab_destroy( p(n))
         call multifab_destroy(gp(n))
         call multifab_destroy( f(n))
         call multifab_destroy(sf(n))
      end do

    end subroutine delete_state

    subroutine initial_iters()

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
                                sforce(n),&
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

    end subroutine initial_iters

    subroutine write_plotfile(istep_to_write)

    integer, intent(in   ) :: istep_to_write

       do n = 1,nlevs
          call make_vorticity(vort(n),uold(n),dx(n,:), &
                              the_bc_tower%bc_tower_array(n))
          call multifab_copy_c(plotdata(n),1           ,unew(n),1,dm)
          call multifab_copy_c(plotdata(n),1+dm        ,snew(n),1,nscal)
          call multifab_copy_c(plotdata(n),1+dm+nscal  ,vort(n),1,1)
          call multifab_copy_c(plotdata(n),1+dm+nscal+1,  gp(n),1,dm)
       end do
       write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
       call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                      mba%pd(1), time, dx(1,:))
    end subroutine write_plotfile
   
    subroutine write_checkfile(istep_to_write)

    integer, intent(in   ) :: istep_to_write

      do n = 1,nlevs
         call multifab_copy_c(chkdata(n),1         ,uold(n),1,dm)
         call multifab_copy_c(chkdata(n),1+dm      ,sold(n),1,nscal)
         call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
      end do
      write(unit=sd_name,fmt='("chk",i4.4)') istep_to_write

      call checkpoint_write(sd_name, chkdata, p, mba%rr, dx, time, dt)

    end subroutine write_checkfile

    subroutine read_command_line()

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

    end subroutine read_command_line

end subroutine varden
