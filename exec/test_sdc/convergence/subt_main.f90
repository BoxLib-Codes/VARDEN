program main

  use BoxLib
  use parallel
  use layout_module
  use bl_prof_module
  use difference

  implicit none

  real(dp_t) :: r1, r2
  real(dp_t) :: cfl = 0.5
  integer    :: sdc_loc(1,4) = (/2032,2064,2128,2256/) 
  integer    :: strang_loc(1,4) = (/9032,9064,9128,9256/) 
!  integer    :: sdc_loc(1,1) = (/2032/) 
!  integer    :: strang_loc(1,1) = (/9032/) 
  integer    :: files(1,2) = (/10,11/)        
  integer    :: i

  call boxlib_initialize()

  r1 = parallel_wtime()

  call bl_prof_initialize(on = .true.)

  !call layout_set_verbosity(1)

  open(10, FILE='l1_comparison.dat', STATUS='NEW')
  open(11, FILE='l2_comparison.dat', STATUS='NEW')

  do i = 1, 1
     call conv(sdc_loc(i,:),strang_loc(i,:),files(i,:),cfl)
  end do

  do i = 10,11
     close(i)
  end do

  call bl_prof_glean("bl_prof_res")

  call bl_prof_finalize()

  r2 = parallel_wtime() - r1

  call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())

  if (parallel_IOProcessor()) print*, 'Run Time = ', r1

  call boxlib_finalize()

end program main
