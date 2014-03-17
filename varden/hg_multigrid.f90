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

    use nodal_divu_module         , only : enforce_outflow_on_divu_rhs

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

    type(box     )  :: pd
    type(  layout)  :: la

    type(mg_tower) :: mgt(mla%nlevel)
    type(multifab), allocatable :: coeffs(:)

    real(dp_t) :: bottom_solver_eps

    logical :: nodal(mla%dim)
    integer :: i, dm, nlevs
    integer :: bottom_solver, bottom_max_iter
    integer :: max_iter
    integer :: min_width
    integer :: max_nlevel
    integer :: nu1, nu2, nub, cycle_type, smoother
    integer :: n
    integer :: max_nlevel_in
    integer :: do_diagnostics
    integer, allocatable :: lo_inflow(:),hi_inflow(:)

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

    bottom_solver = 1

    if ( hg_bottom_solver >= 0 ) then
        if (hg_bottom_solver == 4 .and. nboxes(phi(1)%la) == 1) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with only one grid -- '
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else if (hg_bottom_solver == 4 .and. max_mg_bottom_nlevels < 2) then
           if (parallel_IOProcessor()) then
              print *,'Dont use hg_bottom_solver == 4 with max_mg_bottom_nlevels < 2'
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else
           bottom_solver = hg_bottom_solver
        end if
    end if

    do n = 1,nlevs
       call multifab_build(coeffs(n), mla%la(n), 1, 1)
       call mkcoeffs(rhohalf(n),coeffs(n))
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
                     bottom_solver = bottom_solver, &
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

  end subroutine hg_multigrid

  !   ********************************************************************************** !


  subroutine mkcoeffs(rho,coeffs)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,dm,ng_r,ng_c

      dm = get_dim(rho)
    ng_r = nghost(rho)
    ng_c = nghost(coeffs)

    do i = 1, nfabs(rho)
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (1)
          call mkcoeffs_1d(cp(:,1,1,1), ng_c, rp(:,1,1,1), ng_r)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), ng_c, rp(:,:,1,1), ng_r)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), ng_c, rp(:,:,:,1), ng_r)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_1d(coeffs,ng_c,rho,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:)

    integer :: i,nx

    nx = size(coeffs,dim=1) - 2

    do i = 1,nx
       coeffs(i) = ONE / rho(i)
    end do

  end subroutine mkcoeffs_1d

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,ng_c,rho,ng_r)

    use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:)

    integer :: i,j
    integer :: nx,ny

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2

    do j = 1,ny
       do i = 1,nx
          coeffs(i,j) = ONE / rho(i,j)
       end do
    end do

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,ng_c,rho,ng_r)

      use bl_constants_module

    integer                        :: ng_c,ng_r
    real(kind=dp_t), intent(inout) :: coeffs(1-ng_c:,1-ng_c:,1-ng_c:)
    real(kind=dp_t), intent(in   ) ::    rho(1-ng_r:,1-ng_r:,1-ng_r:)

    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2
    nz = size(coeffs,dim=3) - 2

!$omp parallel do private(i,j,k)
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             coeffs(i,j,k) = ONE / rho(i,j,k)
          end do
       end do
    end do
!$omp end parallel do

  end subroutine mkcoeffs_3d

end module hg_multigrid_module
