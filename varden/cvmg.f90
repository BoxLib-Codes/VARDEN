module cvmg_module

  use bl_types

  implicit none

contains

      function cvmgp(a,b,c) result(r)
        real(dp_t), intent(in) :: a,b,c
        real(dp_t) :: r
        if (c > 0.0_dp_t) then
          r = a
        else
          r = b
        endif
      end function cvmgp

      function cvmgt(a,b,c) result(r)
        real(dp_t), intent(in) :: a,b
        logical, intent(in) :: c
        real(dp_t) :: r
        if (c) then
          r = a
        else
          r = b
        endif
      end function cvmgt

end module cvmg_module
