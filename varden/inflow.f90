module inflow_values

   use bl_types,  only: dp_t
   implicit none

   real(kind=dp_t), parameter ::  VX_INLET = 0.0d0
   real(kind=dp_t), parameter ::  VY_INLET = 0.0d0
   real(kind=dp_t), parameter :: RHO_INLET = 1.0d0

end module inflow_values
