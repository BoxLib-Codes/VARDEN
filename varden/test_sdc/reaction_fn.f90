module reaction_fn

  use bl_types
  use bl_constants_module

  implicit none

  public :: f_rxn

contains 

  ! make sure to adjust integral to 0 -> dt    
  ! the reactions, used as is for strang splitting
  subroutine f_rxn(soln,u,t)  
     use probin_module, only : k_rxn1, k_rxn2

      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t

      real(kind=dp_t) :: u2,u3,u4

      ! first component is density 
      u2 = merge(u(2),0.d0,u(2) > 0.d0)
      u3 = merge(u(3),0.d0,u(3) > 0.d0)
      u4 = merge(u(4),0.d0,u(4) > 0.d0)

      soln(1) =    -k_rxn1*u3*u2 + half*k_rxn2*u4
      soln(2) =    -k_rxn1*u3*u2 + half*k_rxn2*u4
      soln(3) = two*k_rxn1*u3*u2 -      k_rxn2*u4

    end subroutine f_rxn

  end module reaction_fn
