module hg_multigrid_module
 
  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
 
  implicit none
 
  private
 
  public :: hg_multigrid

contains

  subroutine hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                          press_comp,stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)

    use nodal_stencil_fill_module , only : stencil_fill_nodal_all_mglevels
    use ml_solve_module           , only : ml_nd_solve    
    use nodal_divu_module         , only : divu, subtract_divu_from_rh
    use mg_module           
    use probin_module             , only : mg_verbose, cg_verbose, hg_bottom_solver, max_mg_bottom_nlevels
    use stencil_types_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp
    integer        , intent(in   ) :: stencil_type
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    type(multifab ), intent(inout), optional :: divu_rhs(:)

    type(multifab), allocatable :: coeffs(:)

    integer :: dm, nlevs
    integer :: max_iter
    integer :: n
    integer :: do_diagnostics

    type(bl_prof_timer), save :: bpt

    call build(bpt,"hg_multigrid")

    !! Defaults:

    dm    = mla%dim
    nlevs = mla%nlevel

    allocate(coeffs(nlevs))

    if ( mg_verbose >= 3 ) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    max_iter = 100

    do n = 1,nlevs
       call multifab_build(coeffs(n), mla%la(n), 1, 1)

       ! Initialize ghost cells of coeffs to zero so we don't get uninitialized errors
       !   when using coeffs in the calculation of Anorm
       call multifab_setval(coeffs(n),0.d0,all=.true.)

       ! coeffs = 1/rho
       call multifab_setval(coeffs(n),1.d0)
       call multifab_div_div_c(coeffs(n),1,rhohalf(n),1,1,0)

       call multifab_fill_boundary(coeffs(n))
    end do

    ! Set rh = -divu_rhs here. 
    ! inside ml_nd_solve we will add div(u) to result in rh = div(unew) -divu_rhs. 
    if (present(divu_rhs)) then
       do n=1,nlevs
          call multifab_copy_c(rh(n),1,divu_rhs(n),1,1,1)
          call multifab_mult_mult_s(rh(n),-ONE)
       end do
    end if

    ! ********************************************************************************
    ! Call the solver 
    ! ********************************************************************************

    call ml_nd_solve(mla,rh,phi,coeffs,dx,the_bc_tower,press_comp, &
                     add_divu=.true.,u=unew, &
                     eps = rel_solver_eps, &
                     abs_eps = abs_solver_eps, &
                     bottom_solver = hg_bottom_solver, &
                     max_bottom_nlevel = max_mg_bottom_nlevels, &
                     do_diagnostics = do_diagnostics, &
                     verbose = mg_verbose, &
                     cg_verbose = cg_verbose, &
                     max_iter = max_iter, &
                     stencil_type = stencil_type)

    ! ********************************************************************************
    ! Clean-up...
    ! ********************************************************************************

    do n=1,nlevs
       call multifab_destroy(coeffs(n))
    end do

    deallocate(coeffs)

    call destroy(bpt)

  end subroutine hg_multigrid

  !   ********************************************************************************** !

end module hg_multigrid_module
