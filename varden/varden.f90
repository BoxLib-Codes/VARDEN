subroutine varden()
  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use mboxarray_module
  use layout_module
  use multifab_module
  use init_module
  use estdt_module
  use advance_module
  use bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use plotfile_module

  implicit none

  integer    :: narg, farg
  integer    :: max_step,init_iter
  integer    :: plot_int
  integer    :: dm
  real(dp_t) :: cflfac
  real(dp_t) :: visc_coef
  real(dp_t) :: diff_coef
  real(dp_t) :: stop_time
  real(dp_t) :: time,dt,dtold,dt_temp
  real(dp_t) :: unrm
  real(dp_t) :: dx(3)
  integer    :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
  integer    :: istep
  integer    :: n, nlevs, n_plot_comps, nscal

  integer     , allocatable ::   phys_bc(:,:)
  integer     , allocatable :: norm_vel_bc(:,:)
  integer     , allocatable :: tang_vel_bc(:,:)
  integer     , allocatable :: scal_bc(:,:)
  integer     , allocatable :: press_bc(:,:)
  integer     , allocatable :: rrs(:)
  integer     , allocatable :: ref_ratio(:)
  type(layout), allocatable :: la_tower(:)

  type(multifab) :: umac, uedgex, uedgey, uedgez, utrans
  type(multifab) :: sedgex, sedgey, sedgez
  type(multifab), pointer ::       uold(:) => Null()
  type(multifab), pointer ::       unew(:) => Null()
  type(multifab), pointer ::       sold(:) => Null()
  type(multifab), pointer ::       snew(:) => Null()
  type(multifab), pointer ::    rhohalf(:) => Null()
  type(multifab), pointer ::          p(:) => Null()
  type(multifab), pointer ::         gp(:) => Null()
  type(multifab), pointer ::      force(:) => Null()
  type(multifab), pointer :: scal_force(:) => Null()
  type(multifab), pointer ::   plotdata(:) => Null()

  character(len=128) :: fname
  character(len=128) :: probin_env
  character(len=128) :: test_set
  character(len=7) :: sd_name
  character(len=20), allocatable :: var_names(:)
  integer :: un, ierr
  logical :: lexist
  logical :: need_inputs
  logical :: nodal(3)

  type(layout)    :: la
  type(mboxarray) :: mba

  namelist /probin/ stop_time
  namelist /probin/ max_step
  namelist /probin/ plot_int
  namelist /probin/ init_iter
  namelist /probin/ cflfac
  namelist /probin/ visc_coef
  namelist /probin/ diff_coef
  namelist /probin/ test_set
  namelist /probin/ bcx_lo
  namelist /probin/ bcx_hi
  namelist /probin/ bcy_lo
  namelist /probin/ bcy_hi
  namelist /probin/ bcz_lo
  namelist /probin/ bcz_hi

  narg = command_argument_count()

! Defaults
  max_step  = 1
  init_iter = 4
  plot_int  = 0

  nscal = 2

  need_inputs = .TRUE.
  test_set = ''

  bcx_lo         = WALL
  bcy_lo         = WALL
  bcz_lo         = WALL
  bcx_hi         = WALL
  bcy_hi         = WALL
  bcz_hi         = WALL

  call get_environment_variable('PROBIN', probin_env, status = ierr)
  if ( need_inputs .AND. ierr == 0 ) then
     un = unit_new()
     open(unit=un, file = probin_env, status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .FALSE.
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
        need_inputs = .FALSE.
     end if
  end if

  inquire(file = 'inputs_varden', exist = lexist)
  if ( need_inputs .AND. lexist ) then
     un = unit_new()
     open(unit=un, file = 'inputs_varden', status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .FALSE.
  end if

  if ( .true. ) then
     do while ( farg <= narg )
        call get_command_argument(farg, value = fname)
        select case (fname)

        case ('--cfl')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) cflfac

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

        case ('--test_set')
           farg = farg + 1
           call get_command_argument(farg, value = test_set)

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

  call read_a_hgproj_grid(mba, test_set)
  nlevs = mba%nlevel
  dm = mba%dim

  allocate(var_names(2*dm+nscal))

  var_names(1) = "x_vel"
  var_names(2) = "y_vel"
  if (dm > 2) var_names(3) = "z_vel"
  var_names(dm+1) = "density"
  if (nscal > 1) var_names(dm+2) = "tracer"
  var_names(dm+nscal+1) = "gpx"
  var_names(dm+nscal+2) = "gpy"
  if (dm > 2) var_names(dm+nscal+3) = "gpz"

  allocate(la_tower(nlevs))
  allocate(ref_ratio(dm))
  call build(la_tower(1), mba%bas(1))
  do n = 2, nlevs
     ref_ratio = mba%rr(n-1,:)
     call layout_build_pn(la_tower(n), la_tower(n-1), mba%bas(n), ref_ratio)
  end do

  allocate(uold(nlevs),unew(nlevs))
  allocate(sold(nlevs),snew(nlevs))
  allocate(p(nlevs),gp(nlevs),rhohalf(nlevs))
  allocate(force(nlevs),scal_force(nlevs))

  nodal = .true.

  do n = nlevs,1,-1
     call multifab_build(   uold(n), la_tower(n),    dm, 3)
     call multifab_build(   unew(n), la_tower(n),    dm, 3)
     call multifab_build(   sold(n), la_tower(n), nscal, 3)
     call multifab_build(   snew(n), la_tower(n), nscal, 3)

     call multifab_build(scal_force(n), la_tower(n), nscal, 1)
     call multifab_build(     force(n), la_tower(n),    dm, 1)
     call multifab_build(         p(n), la_tower(n),     1, 1, nodal)
     call multifab_build(        gp(n), la_tower(n),    dm, 1)
     call multifab_build(   rhohalf(n), la_tower(n),     1, 1)

     call setval(uold(n),0.0_dp_t, all=.true.)
     call setval(unew(n),0.0_dp_t, all=.true.)
     call setval(sold(n),0.0_dp_t, all=.true.)
     call setval(snew(n),0.0_dp_t, all=.true.)

     call setval(        gp(n),0.0_dp_t, all=.true.)
     call setval(   rhohalf(n),1.0_dp_t, all=.true.)
     call setval(     force(n),0.0_dp_t, all=.true.)
     call setval(scal_force(n),0.0_dp_t, all=.true.)
  end do

  la = la_tower(1)

  call multifab_build(umac       , la,    dm, 1)
  call multifab_build(uedgex     , la,    dm, 1)
  call multifab_build(uedgey     , la,    dm, 1)
  call multifab_build(uedgez     , la,    dm, 1)
  call multifab_build(sedgex     , la, nscal, 1)
  call multifab_build(sedgey     , la, nscal, 1)
  call multifab_build(sedgez     , la, nscal, 1)
  call multifab_build(utrans     , la,    dm, 2)

  call setval(umac,0.0_dp_t, all=.true.)

  dx = 1.d0

  time = 0.d0
  dt = 0.d0

  allocate(rrs(nlevs-1))
  if (nlevs > 1) rrs = 2

! Put the bc values from the inputs file into phys_bc
  allocate(phys_bc(dm,2))
  phys_bc(1,1) = bcx_lo
  phys_bc(1,2) = bcx_hi
  if (dm > 1) then
     phys_bc(2,1) = bcy_lo
     phys_bc(2,2) = bcy_hi
  end if
  if (dm > 2) then
     phys_bc(3,1) = bcz_lo
     phys_bc(3,2) = bcz_hi
  end if

! Translate the phys_bc values into bc's for each component
  allocate(norm_vel_bc(dm,2))
  allocate(tang_vel_bc(dm,2))
  allocate(scal_bc(dm,2))
  allocate(press_bc(dm,2))
  call define_bcs(phys_bc,norm_vel_bc,tang_vel_bc,scal_bc,press_bc,visc_coef)

  dt_temp = 1.0_dp_t
  do n = 1,nlevs
     call initdata(uold(n),sold(n),dx)
     call multifab_fill_boundary(sold(n))

!    We do this here in order to set any Dirichlet boundary conditions.
     unew(n) = uold(n)
     snew(n) = sold(n)

     unrm = norm_inf(uold(n))
     print *,'MAX OF UOLD BEFORE PROJ ',unrm,' AT LEVEL ',n

!    Note that we use rhohalf, filled with 1 at this point, as a temporary
!       in order to do a constant-density initial projection.
     call hgproject(uold(n),rhohalf(n),p(n),gp(n),dx,dt_temp,phys_bc,press_bc)
     call multifab_fill_boundary(uold(n))
  end do

  do n = 1,nlevs
     call setval(gp(n)  ,0.0_dp_t, all=.true.)
  end do

  do n = 1,nlevs
     print *,'MAX OF UOLD AFTER PROJ ',norm_inf(uold(n)),' AT LEVEL ',n
  end do
  print *,' '

  dtold = dt

  do n = 1,nlevs
     call estdt(uold(n),sold(n),dx,cflfac,dtold,dt)
  end do

  if (init_iter > 0) then
     print *,'DOING ',init_iter,' INITIAL ITERATIONS ' 
     do istep = 1,init_iter
      do n = 1,nlevs
         call advance(uold(n),unew(n),sold(n),snew(n),rhohalf(n),&
                      umac,uedgex,uedgey,uedgez, &
                      sedgex,sedgey,sedgez, &
                      utrans,p(n),gp(n),force(n),scal_force(n),&
                      dx,time,dt,phys_bc,norm_vel_bc,tang_vel_bc, &
                      scal_bc,press_bc,visc_coef,diff_coef)
      end do
     end do

     do n = 1,nlevs
        print *,'MAX OF UOLD ',norm_inf(uold(n)),' AT LEVEL ',n
     end do
  end if
  print *,' '

  allocate(plotdata(nlevs))
  n_plot_comps = 2*dm + nscal
  do n = 1,nlevs
     call multifab_build(plotdata(n), la_tower(n), n_plot_comps, 0)
  end do


! This writes a plotfile.
  istep = 0
  do n = 1,nlevs
     call multifab_copy_c(plotdata(n),1         ,uold(n),1,dm)
     call multifab_copy_c(plotdata(n),1+dm      ,sold(n),1,nscal)
     call multifab_copy_c(plotdata(n),1+dm+nscal,  gp(n),1,dm)
  end do
  write(unit=sd_name,fmt='("plt",i4.4)') istep
  call fabio_ml_multifab_write_d(plotdata, rrs, sd_name, var_names)

  if (max_step > 0) then
     print *,'DOING TIME ADVANCE '
       do istep = 1,max_step
       do n = 1,nlevs
          call multifab_fill_boundary(uold(n))
          call multifab_fill_boundary(sold(n))

          dtold = dt
          call estdt(uold(n),sold(n),dx,cflfac,dtold,dt)

          print *,'MAX OF UOLD ',norm_inf(uold(n)),' AT TIME ',time

          call advance(uold(n),unew(n),sold(n),snew(n),rhohalf(n), &
                       umac,uedgex,uedgey,uedgez, &
                       sedgex,sedgey,sedgez, &
                       utrans,p(n),gp(n),force(n),scal_force(n),&
                       dx,time,dt,phys_bc,norm_vel_bc,tang_vel_bc, &
                       scal_bc,press_bc,visc_coef,diff_coef)

          time = time + dt

          write(6,1000) istep,time,dt
          print *,'MAX OF UNEW ',norm_inf(unew(n)),' AT TIME ',time
          print *,'MAX OF SNEW ',norm_inf(snew(n)),' AT TIME ',time
          print *,' '

         end do

         do n = 1,nlevs
           call multifab_copy_c(uold(n),1,unew(n),1,dm)
           call multifab_copy_c(sold(n),1,snew(n),1,nscal)
         end do

         if (plot_int > 0 .and. mod(istep,plot_int) .eq. 0) then
            do n = 1,nlevs
               call multifab_copy_c(plotdata(n),1         ,unew(n),1,dm)
               call multifab_copy_c(plotdata(n),1+dm      ,snew(n),1,nscal)
               call multifab_copy_c(plotdata(n),1+dm+nscal,  gp(n),1,dm)
            end do
            write(unit=sd_name,fmt='("plt",i4.4)') istep
            call fabio_ml_multifab_write_d(plotdata, rrs, sd_name, var_names)
         end if

       end do
1000   format('STEP = ',i4,1x,' TIME = ',f14.10,1x,'DT = ',f14.9)

!      This writes a plotfile.
       do n = 1,nlevs
          call multifab_copy_c(plotdata(n),1         ,unew(n),1,dm)
          call multifab_copy_c(plotdata(n),1+dm      ,snew(n),1,nscal)
          call multifab_copy_c(plotdata(n),1+dm+nscal,  gp(n),1,dm)
       end do
       write(unit=sd_name,fmt='("plt",i4.4)') max_step
       call fabio_ml_multifab_write_d(plotdata, rrs, sd_name, var_names)
  end if

  do n = 1,nlevs
     call multifab_destroy(uold(n))
     call multifab_destroy(unew(n))
     call multifab_destroy(sold(n))
     call multifab_destroy(snew(n))
     call multifab_destroy(gp(n))
     call multifab_destroy( p(n))
     call multifab_destroy(rhohalf(n))
     call multifab_destroy(force(n))
     call multifab_destroy(scal_force(n))
     call multifab_destroy(plotdata(n))
  end do

  deallocate(norm_vel_bc)
  deallocate(tang_vel_bc)
  deallocate(scal_bc)
  deallocate(press_bc)

end subroutine varden
