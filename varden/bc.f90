module bc_module

  implicit none

  integer, parameter :: UNDEFINED = -3
  integer, parameter :: INTERIOR  = -2
  integer, parameter :: PERIODIC  =  1
  integer, parameter :: WALL      =  2
  integer, parameter :: INLET     =  3
  integer, parameter :: OUTLET    =  4

end module bc_module
