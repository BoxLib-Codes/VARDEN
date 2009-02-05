program main

  use BoxLib
  use parallel
  use layout_module
  use bl_prof_module
  use convergence

  implicit none

  real(dp_t) :: r1, r2
  integer    :: loc(3,4) = (/0032,1032,2032, &
                             0064,1064,2064, &
                             0128,1128,2128, &
                             0256,1256,2256/) 
  integer    :: files(3,2) = (/10,12,14, &
                               11,13,15/)        
  integer    :: i

  call boxlib_initialize()

  r1 = parallel_wtime()

  call bl_prof_initialize(on = .true.)

  !call layout_set_verbosity(1)

  open(10, FILE='l1_sdc0.dat', STATUS='NEW')
  open(11, FILE='l2_sdc0.dat', STATUS='NEW')

  open(12, FILE='l1_sdc1.dat', STATUS='NEW')
  open(13, FILE='l2_sdc1.dat', STATUS='NEW')

  open(14, FILE='l1_sdc2.dat', STATUS='NEW')
  open(15, FILE='l2_sdc2.dat', STATUS='NEW')

  do i = 1, 3
     call conv(loc(i,:),files(i,:))
  end do

  do i = 10,15
     close(i)
  end do

  call bl_prof_glean("bl_prof_res")

  call bl_prof_finalize()

  r2 = parallel_wtime() - r1

  call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())

  if (parallel_IOProcessor()) print*, 'Run Time = ', r1

  call boxlib_finalize()

end program main
