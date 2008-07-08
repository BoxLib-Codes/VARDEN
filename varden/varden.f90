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
  use advance_module

  use probin_module, only : dim_in, max_levs, nlevs, ng_cell, ng_grow, init_iter, max_step, &
                            stop_time, restart, chk_int, plot_int, regrid_int, init_shrink, &
                            fixed_dt, bcx_lo, bcy_lo, bcz_lo, bcx_hi, bcy_hi, bcz_hi, &
                            n_cellx, n_celly, n_cellz, prob_lo_x, prob_lo_y, prob_lo_z, &
                            prob_hi_x, prob_hi_y, prob_hi_z, ref_ratio, pmask_xyz, &
                            fixed_grids, grids_file_name, max_grid_size, &
                            do_initial_projection, grav, probin_init

  implicit none

  integer    :: dm
  real(dp_t) :: time,dt,dtold,dt_lev,dt_temp
  integer    :: istep
  integer    :: i, n
  integer    :: n_chk_comps
  integer    :: last_plt_written, last_chk_written
  integer    :: init_step
  integer    :: press_comp, vort_comp

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
  type(multifab), allocatable :: ext_vel_force(:)
  type(multifab), allocatable :: ext_scal_force(:)
  type(multifab), allocatable :: plotdata(:)
  type(multifab), pointer     ::  chkdata(:)
  type(multifab), pointer     ::    chk_p(:)

  ! Regridding quantities
  type(multifab), allocatable ::   uold_rg(:)
  type(multifab), allocatable ::   sold_rg(:)
  type(multifab), allocatable ::     gp_rg(:)
  type(multifab), allocatable ::      p_rg(:)

  integer,allocatable      :: lo(:),hi(:)

  character(len=7 ) :: sd_name
  character(len=20), allocatable :: plot_names(:)

  type(ml_layout) :: mla_temp
  type(boxarray)  :: ba

  type(bc_tower) ::  the_bc_tower
  type(bc_level) ::  bc

  last_plt_written = -1
  last_chk_written = -1

  call probin_init()

  dm = dim_in
  press_comp = dm + nscal + 1

  allocate(pmask(dm))
  pmask = .FALSE.

  pmask = pmask_xyz(1:dm)

  allocate(nodal(dm))
  nodal = .true.

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
  ! Initialize prob_hi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(prob_hi(dm))
  prob_hi(1) = prob_hi_x
  if (dm > 1) prob_hi(2) = prob_hi_y
  if (dm > 2) prob_hi(3) = prob_hi_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the grids and the data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart >= 0) then

     call initialize_from_restart()

  else if (fixed_grids /= '') then

     call initialize_with_fixed_grids()

  else  ! Adaptive gridding

     call initialize_with_adaptive_grids()
     if ( parallel_IOProcessor() .and. verbose) &
        call print(mla,"MLA OUT OF INITIAL GRIDDING ROUTINE")

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Allocate new-time state and temp variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(unew(nlevs),snew(nlevs))
  allocate(ext_vel_force(nlevs),ext_scal_force(nlevs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initial projection if not restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then

     time    = ZERO
     dt_temp = ONE

     ! Note that we use rhohalf, filled with 1 at this point, as a temporary
     ! in order to do a constant-density initial projection.
     if (do_initial_projection > 0) then
        allocate(rhohalf(nlevs))
        do n = 1,nlevs
          call multifab_build(rhohalf(n), mla%la(n),1, 1)
          call setval(rhohalf(n),ONE, all=.true.)
        end do
        call hgproject(initial_projection,mla,uold,uold,rhohalf,p,gp,dx,dt_temp, &
                       the_bc_tower,press_comp)
         do n = 1,nlevs
           call multifab_destroy(rhohalf(n))
           call setval( p(n)  ,0.0_dp_t, all=.true.)
           call setval(gp(n)  ,0.0_dp_t, all=.true.)
        end do
        deallocate(rhohalf)
     end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write the grids into a "grdlog" file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (grids_file_name /= '') &
     call write_grids(grids_file_name,mla,0)

  call make_temps(mla)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Impose bc's on uold and copy to unew
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute the time step.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dt = 1.d20
  dtold = dt
  do n = 1,nlevs
     call estdt(n,uold(n),sold(n),gp(n),ext_vel_force(n),dx(n,:), &
                dtold,dt_lev)
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
  ! Make temporaries for regridding
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (regrid_int > 0) then
     allocate(uold_rg(nlevs),sold_rg(nlevs),p_rg(nlevs),gp_rg(nlevs))
     call build(mla_temp,mla%mba, pmask)
     do n = 1, mla_temp%nlevel
        call make_new_state(mla_temp%la(n),uold_rg(n),sold_rg(n),gp_rg(n),p_rg(n))
     enddo
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the real integration.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ((max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0)) then

     do istep = init_step, max_step

!       if ( verbose > 0 ) then
!          if ( parallel_IOProcessor() ) then
!             print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
!             print*, ' '
!          end if
!          call print(multifab_mem_stats(),    "    multifab")
!          call print(fab_mem_stats(),         "         fab")
!          call print(boxarray_mem_stats(),    "    boxarray")
!          call print(layout_mem_stats(),      "      layout")
!          call print(boxassoc_mem_stats(),    "    boxassoc")
!          call print(fgassoc_mem_stats(),     "     fgassoc")
!          call print(syncassoc_mem_stats(),   "   syncassoc")
!          call print(copyassoc_mem_stats(),   "   copyassoc")
!          call print(fluxassoc_mem_stats(),   "   fluxassoc")
!          if ( parallel_IOProcessor() ) print*, ''
!       end if

        if (nlevs > 1 .and. regrid_int > 0 .and. &
            (mod(istep-1,regrid_int) .eq. 0) ) then

           call delete_temps()
!          call regrid_each_time()
           call make_temps(mla)

        end if  

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
                   dtold,dt_lev)
              dt = min(dt,dt_lev)
           end do
           if (fixed_dt > 0.d0) dt = fixed_dt
           if (stop_time >= 0.d0) then
              if (time+dt > stop_time) dt = stop_time - time
           end if
        end if

        call advance_timestep(istep,mla,sold,uold,snew,unew,gp,p,ext_vel_force,ext_scal_force,&
                              the_bc_tower,dt,time,dx,press_comp,regular_timestep)

        do n = 1,nlevs
           call multifab_copy_c(uold(n),1,unew(n),1,dm)
           call multifab_copy_c(sold(n),1,snew(n),1,nscal)
        end do

        time = time + dt

!       if ( verbose > 0 ) then
!          if ( parallel_IOProcessor() ) then
!             print *, 'MEMORY STATS AT END OF TIMESTEP ', istep
!             print*, ' '
!          end if
!          call print(multifab_mem_stats(),    "    multifab")
!          call print(fab_mem_stats(),         "         fab")
!          call print(boxarray_mem_stats(),    "    boxarray")
!          call print(layout_mem_stats(),      "      layout")
!          call print(boxassoc_mem_stats(),    "    boxassoc")
!          call print(fgassoc_mem_stats(),     "     fgassoc")
!          call print(syncassoc_mem_stats(),   "   syncassoc")
!          call print(copyassoc_mem_stats(),   "   copyassoc")
!          call print(fluxassoc_mem_stats(),   "   fluxassoc")
!          if ( parallel_IOProcessor() ) print*, ''
!        end if

         if ( parallel_IOProcessor() ) then
            write(6,1000) istep,time,dt
         end if

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

      end do ! istep loop

     if (istep > max_step) istep = max_step

     if (last_plt_written .ne. istep .and. plot_int > 0) call write_plotfile(istep)
     if (last_chk_written .ne. istep .and. chk_int  > 0) call write_checkfile(istep)

  end if
  
  call delete_state(uold,sold,gp,p)
  call delete_temps()
  if (regrid_int > 0 .and. max_levs > 1) &
     call delete_state(uold_rg,sold_rg,gp_rg,p_rg)

  call bc_tower_destroy(the_bc_tower)

  call destroy(mla)
  if (regrid_int > 0 .and. max_levs > 1) &
     call destroy(mla_temp)

  if ( verbose > 0 ) then
     if ( parallel_IOProcessor() ) then
        print *, 'MEMORY STATS AT END OF RUN '
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

1000 format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)

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

    do n = nlevs,1,-1
       call multifab_build(   unew(n), mla_loc%la(n),    dm, ng_cell)
       call multifab_build(   snew(n), mla_loc%la(n), nscal, ng_cell)
       call multifab_build(ext_vel_force(n),  mla_loc%la(n),    dm, 1)
       call multifab_build(ext_scal_force(n), mla_loc%la(n), nscal, 1)

       call setval(   unew(n),ZERO, all=.true.)
       call setval(   snew(n),ZERO, all=.true.)
       call setval(ext_vel_force(n) ,ZERO, 1,dm-1,all=.true.)
       call setval(ext_vel_force(n) ,grav,dm,   1,all=.true.)
       call setval(ext_scal_force(n),ZERO, all=.true.)

    end do

  end subroutine make_temps

  subroutine delete_temps()

    do n = 1,nlevs
       call multifab_destroy(unew(n))
       call multifab_destroy(snew(n))
       call multifab_destroy(ext_vel_force(n))
       call multifab_destroy(ext_scal_force(n))
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

       call advance_timestep(istep,mla,sold,uold,snew,unew,gp,p,ext_vel_force,ext_scal_force,&
                             the_bc_tower,dt,time,dx,press_comp,pressure_iters)

    end do

  end subroutine initial_iters

  subroutine initialize_from_restart()

     type(ml_boxarray)         :: mba

     call fill_restart_data(restart,mba,chkdata,chk_p,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     allocate(uold(nlevs),sold(nlevs),gp(nlevs),p(nlevs))
     do n = 1,nlevs
        call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
        call multifab_build(   sold(n), mla%la(n), nscal, ng_cell)
        call multifab_build(     gp(n), mla%la(n),    dm, ng_grow)
        call multifab_build(      p(n), mla%la(n),     1, ng_grow, nodal)
     end do
     do n = 1,nlevs
        call multifab_copy_c(uold(n),1,chkdata(n),1         ,dm)
        call multifab_copy_c(sold(n),1,chkdata(n),1+dm      ,nscal)
        call multifab_copy_c(  gp(n),1,chkdata(n),1+dm+nscal,dm)
        call multifab_copy_c(   p(n),1,  chk_p(n),1         ,1)
        !
        ! The layouts for chkdata and chk_p are built standalone, level
        ! by level, and need to be destroy()d as such as well.
        !
        call destroy(chkdata(n)%la)
        call destroy(chk_p(n)%la)
        call multifab_destroy(chkdata(n))
        call multifab_destroy(chk_p(n))
     end do
     deallocate(chkdata,chk_p)

     call initialize_dx(mba,nlevs)

     call initialize_bc(nlevs)
     do n = 1,nlevs
        call bc_tower_level_build( the_bc_tower,n,mla%la(n))
     end do

  end subroutine initialize_from_restart

  subroutine initialize_with_fixed_grids()

     type(ml_boxarray)         :: mba

     call read_a_hgproj_grid(mba, fixed_grids)
     call ml_layout_build(mla,mba,pmask)

     ! check for proper nesting
     if (.not. ml_boxarray_properly_nested(mla%mba, ng_cell, pmask)) &
         call bl_error('fixed_grids not properly nested')

     nlevs = mla%nlevel
     allocate(uold(nlevs),sold(nlevs),p(nlevs),gp(nlevs))

     do n = 1,nlevs
        call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
        call multifab_build(   sold(n), mla%la(n), nscal, ng_cell)
        call multifab_build(     gp(n), mla%la(n),    dm, ng_grow)
        call multifab_build(      p(n), mla%la(n),     1, ng_grow, nodal)
     end do

     call initialize_dx(mba,nlevs)

     call initialize_bc(nlevs)
     do n = 1,nlevs
        call bc_tower_level_build( the_bc_tower,n,mla%la(n))
     end do

     call initdata(nlevs,uold,sold,dx,prob_hi,the_bc_tower%bc_tower_array,mla)

  end subroutine initialize_with_fixed_grids

  subroutine initialize_with_adaptive_grids()

     integer                   :: buf_wid
     type(layout), allocatable :: la_array(:)
     type(box)   , allocatable :: bxs(:)
     type(boxarray)            :: ba_new,ba_new_comp,ba_old_comp
     type(boxarray)            :: ba_newest
     type(ml_boxarray)         :: mba
     type(list_box)            :: bl

     logical  :: new_grid
     integer  :: nl, buff

     buff = 2

     buf_wid = regrid_int

     ! set up hi & lo to carry indexing info
     allocate(lo(dm),hi(dm))
     lo(:) = 0
     hi(1) = n_cellx-1
     if (dm > 1) then   
        hi(2) = n_celly - 1        
        if (dm > 2)  then
           hi(3) = n_cellz -1
        endif
     endif

     ! mba is big enough to hold max_levs levels
     call ml_boxarray_build_n(mba,max_levs,dm)
     do n = 1, max_levs-1
        mba%rr(n,:) = ref_ratio
     enddo

     if (max_levs > 1) allocate(la_array(max_levs))
     allocate(bxs(max_levs))
     allocate(uold(max_levs),sold(max_levs),p(max_levs),gp(max_levs))

       ! Build the level 1 boxarray
     call box_build_2(bxs(1),lo,hi)
     call boxarray_build_bx(mba%bas(1),bxs(1))
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! build pd(:)
     mba%pd(1) = bxs(1)
     do n = 2, max_levs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     enddo

     ! Need to build pd before making dx
     call initialize_dx(mba,max_levs)

     ! Initialize bc's.
     call initialize_bc(max_levs)

     if (max_levs > 1) then

        ! Build the level 1 layout.
        call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

        ! Build the level 1 data only.
        call make_new_state(la_array(1),uold(1),sold(1),gp(1),p(1)) 

        ! Define bc_tower at level 1.
        call bc_tower_level_build(the_bc_tower,1,la_array(1))

        ! Initialize the level 1 data only.
        call initdata_on_level(uold(1),sold(1),dx(1,:),prob_hi,the_bc_tower%bc_tower_array(1),la_array(1))

        new_grid = .true.
        nl = 1

        do while ( (nl .lt. max_levs) .and. (new_grid) )

           ! Do we need finer grids?
           call make_new_grids(la_array(nl),la_array(nl+1),sold(nl),dx(nl,1),buf_wid,&
                               ref_ratio,nl,max_grid_size,new_grid)
        
           if (new_grid) then

              call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

              ! Build the level nl+1 data only.
              call make_new_state(la_array(nl+1),uold(nl+1),sold(nl+1),gp(nl+1),p(nl+1)) 

              ! Define bc_tower at level nl+1.
              call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
            
             ! fills the physical region of each level with problem data (blob now)
              call initdata_on_level(uold(nl+1),sold(nl+1),dx(nl+1,:),prob_hi,the_bc_tower%bc_tower_array(nl+1),la_array(nl+1))

              nlevs = nl+1
              nl = nl + 1

           endif ! if (new_grid) 

      enddo          

      do n = 1,nl
         call destroy(sold(n))
         call destroy(uold(n))
         call destroy(gp(n))
         call destroy(p(n))
      end do

      nlevs = nl

      ! check for proper nesting
      if (nlevs .ge. 3) then
        
         nl = nlevs - 1
         new_grid = .true.

         do while ( (nl .ge. 2) .and. (new_grid) )

            if (.not. ml_boxarray_properly_nested(mba, ng_cell, pmask, nl, nl+1)) then

                new_grid = .true.
               
                if ( parallel_IOProcessor() .and. verbose) then
                   print *,' '
                   print *,'LEVEL ',nl+1,' grids are not properly nested '
                   print *,' '
                end if

                ! Buffer returns a boxarray "ba_new" that contains everything at level nl 
                !  that the level nl+1 level will need for proper nesting
                call buffer(nl,la_array(nl+1),ba_new,ref_ratio,ng_cell)
                call boxarray_intersection(ba_new,mba%pd(nl))

                ! Make sure the new grids start on even and end on odd
                call boxarray_coarsen(ba_new, ref_ratio)
                call boxarray_refine (ba_new, ref_ratio)

                ! Merge the new boxarray "ba_new" with the existing box_array 
                ! mba%bas(nl) so that we get the union of points.
                call boxarray_complementIn(ba_old_comp,mba%pd(nl),mba%bas(nl))
                call boxarray_intersection(ba_old_comp,ba_new)
                do i = 1, mba%bas(nl)%nboxes
                   call push_back(bl,  mba%bas(nl)%bxs(i))
                end do
                do i = 1, ba_old_comp%nboxes
                   call push_back(bl, ba_old_comp%bxs(i))
                end do
                call build(ba_newest,bl)
                call destroy(bl)
                call boxarray_simplify(ba_newest)
                call boxarray_maxsize(ba_newest,max_grid_size)

                ! Replace mba%bas(nl) by ba_new
                call destroy(mba%bas(nl))
                call copy(mba%bas(nl),ba_newest)
                call destroy(ba_newest)

                ! Double check we got the proper nesting right
                if (.not. ml_boxarray_properly_nested(mba, ng_cell, pmask, nl, nl+1)) &
                  call bl_error('Still not properly nested, darn it')

                ! Destroy the old layout and build a new one.
                call destroy(la_array(nl))
                call layout_build_ba(la_array(nl),mba%bas(nl),mba%pd(nl),pmask)

            else

                new_grid = .false.

            endif  !if not properly nested

            nl = nl - 1

         enddo ! do while
      end if ! if (nlevs .ge. 3)

   end if ! end if (maxlev > 1)

   call ml_layout_restricted_build(mla,mba,nlevs,pmask)

   nlevs = mla%nlevel

   do n = 1,nlevs
      call make_new_state(mla%la(n),uold(n),sold(n),gp(n),p(n)) 
   end do

   call initdata(nlevs,uold,sold,dx,prob_hi,the_bc_tower%bc_tower_array,mla)

   call destroy(mba)
   deallocate(bxs)

  end subroutine initialize_with_adaptive_grids

  subroutine initialize_bc(num_levs)

     integer, intent(in) :: num_levs
     integer, allocatable :: domain_phys_bc(:,:)

     ! Define the physical boundary conditions on the domain
     allocate(domain_phys_bc(dm,2))
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bcx_lo
     domain_phys_bc(1,2) = bcx_hi
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bcy_lo
        domain_phys_bc(2,2) = bcy_hi
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bcz_lo
        domain_phys_bc(3,2) = bcz_hi
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.
     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

     deallocate(domain_phys_bc)

  end subroutine initialize_bc

  subroutine initialize_dx(mba,num_levs)
  
     type(ml_boxarray), intent(in) :: mba
     integer          , intent(in) :: num_levs

     allocate(dx(num_levs,dm))

     do i = 1,dm
        dx(1,i) = prob_hi(i) / float(mba%pd(1)%hi(i)-mba%pd(1)%lo(i)+1)
     end do
     do n = 2,num_levs
        dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
     end do

  end subroutine initialize_dx

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
                                   mla%mba%pd(1), time, dx(1,:))

    do n = 1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)

  end subroutine write_plotfile

  subroutine write_checkfile(istep_to_write)

    integer, intent(in   ) :: istep_to_write

    allocate(chkdata(nlevs))
    n_chk_comps = 2*dm + nscal
    do n = 1,nlevs
       call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
       call multifab_copy_c(chkdata(n),1         ,uold(n),1,dm)
       call multifab_copy_c(chkdata(n),1+dm      ,sold(n),1,nscal)
       call multifab_copy_c(chkdata(n),1+dm+nscal,  gp(n),1,dm)
    end do
    write(unit=sd_name,fmt='("chk",i4.4)') istep_to_write

    call checkpoint_write(nlevs, sd_name, chkdata, p, mla%mba%rr, dx, time, dt, verbose)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
    end do
    deallocate(chkdata)

  end subroutine write_checkfile

  subroutine write_grids(grids_file_name,mla,nstep)

    character(len=128), intent(inout) :: grids_file_name
    type(ml_layout)   , intent(in   ) :: mla
    integer           , intent(in   ) :: nstep

    integer        :: i,d,un,nb,tp(mla%dim)
    type(boxarray) :: mba
    type(box)      :: bx

    un = 11
    tp = 0
   
    if ( parallel_IOProcessor() ) then
       open(un,file=grids_file_name, position='append')
       write(unit=un, fmt='("At step ",i2,":")') nstep
!      write(unit=un, fmt='(i1," levels ")') mla%mba%nlevel
       write(unit=un, fmt='(i2)') mla%mba%nlevel
       do n = 1, mla%mba%nlevel
          nb = mla%mba%bas(n)%nboxes
          bx = ml_layout_get_pd(mla,n)
          write(unit=un, fmt='("   (")', advance = 'no') 
          write(unit=un, fmt='("(" 3(I0,:,", "))', advance = 'no') bx%lo(1:bx%dim)
          write(unit=un, fmt='(") (", 3(I0,:,", "))', advance = 'no') bx%hi(1:bx%dim)
          write(unit=un, fmt='(") (" 3(I0,:,","))', advance = 'no') tp(1:bx%dim)
          write(unit=un, fmt='("))")', advance = 'no' )
          write(unit=un, fmt='(" ",i4)', advance = 'yes') nb
          do i = 1, nb
             bx = mla%mba%bas(n)%bxs(i)
             tp = 0
             write(unit=un, fmt='("      (")', advance = 'no') 
             write(unit=un, fmt='("(", 3(I0,:,", "))', advance = 'no') bx%lo(1:bx%dim)
             write(unit=un, fmt='(") (" 3(I0,:,", "))', advance = 'no') bx%hi(1:bx%dim)
             write(unit=un, fmt='(") (" 3(I0,:,","))', advance = 'no') tp(1:bx%dim)
             write(unit=un, fmt='("))")', advance = 'no' )
             write(unit=un, fmt='(" ")')
          end do
       end do
       write(unit=un, fmt='(" ")')
       close(un)
    end if

  end subroutine write_grids

end subroutine varden
