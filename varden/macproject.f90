module macproject_module

  use bl_types
  use define_bc_module
  use multifab_module
  use boxarray_module
  use stencil_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t

contains 

subroutine macproject(umac,rho,dx,press_bc,domain_press_bc,verbose,mg_verbose)

  type(multifab), intent(inout) :: umac
  type(multifab), intent(in   ) :: rho
  integer       , intent(in   ) :: domain_press_bc(:,:)
  type(bc_level), intent(in   ) :: press_bc
  integer       , intent(in   ) :: verbose,mg_verbose
  real(dp_t) :: dx(:)

! Local  
  type(multifab) :: rh,phi,alpha,beta
  type(layout) :: la
  integer  :: dm,stencil_order

  la = umac%la
  dm = rho%dim

  call multifab_build(   rh, la,  1, 0)
  call multifab_build(  phi, la,  1, 1)
  call multifab_build(alpha, la,  1, 1)
  call multifab_build( beta, la, dm, 1)

  if (verbose .eq. 1) then
     print *,' '
     print *,'... begin mac_projection ... '
  end if

  call divumac(umac,rh,dx,verbose)

  call mk_mac_coeffs(rho,beta)

  call setval(alpha,ZERO,all=.true.)
  call setval(  phi,ZERO,all=.true.)

  stencil_order = 1

  call mac_multigrid(la,rh,phi,alpha,beta,dx,domain_press_bc,stencil_order,mg_verbose)

  call mkumac(umac,phi,beta,dx,press_bc,verbose)

  if (verbose .eq. 1) then
     print *,'...   end mac_projection ... '
     print *,' '
  end if

  call multifab_destroy(rh)
  call multifab_destroy(phi)
  call multifab_destroy(alpha)
  call multifab_destroy(beta)

  contains

    subroutine divumac(umac,rh,dx,verbose)

      type(multifab) , intent(in   ) :: umac
      type(multifab) , intent(inout) :: rh
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: verbose
 
      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      real(kind=dp_t)          :: rhsum,rhmax
      real(kind=dp_t)          :: local_sum,local_max
      integer :: i,dm

      dm = rh%dim

      rhmax = ZERO
      rhsum = ZERO

      do i = 1, umac%nboxes
         if ( multifab_remote(umac, i) ) cycle
         ump => dataptr(umac, i)
         rhp => dataptr(rh  , i)
         select case (dm)
            case (2)
              call divumac_2d(ump(:,:,1,:), rhp(:,:,1,1), dx, local_max, local_sum)
            case (3)
              call divumac_3d(ump(:,:,:,:), rhp(:,:,:,1), dx, local_max, local_sum)
         end select
         rhmax = max(rhmax,local_max)
         rhsum = rhsum + local_sum
      end do
      if (verbose .eq. 1) print *,'MAX/SUM OF RHS ',rhmax,rhsum

    end subroutine divumac

    subroutine divumac_2d(umac,rh,dx,rhmax,rhsum)

      real(kind=dp_t), intent(in   ) :: umac(-2:,-2:,:)
      real(kind=dp_t), intent(inout) ::   rh( 0:, 0:)
      real(kind=dp_t), intent(in   ) ::   dx(:)
      real(kind=dp_t), intent(  out) :: rhsum,rhmax

      integer :: i,j

      rhsum = 0.0_dp_t
      rhmax = 0.0_dp_t
      do j = 0, size(rh,dim=2)-1
      do i = 0, size(rh,dim=1)-1
         rh(i,j) = (umac(i+1,j,1) - umac(i,j,1)) / dx(1) + &
                   (umac(i,j+1,2) - umac(i,j,2)) / dx(2)
         rhsum = rhsum + rh(i,j)
         rhmax = max(rhmax,abs(rh(i,j)))
         rh(i,j) = -rh(i,j)
      end do
      end do

    end subroutine divumac_2d

    subroutine divumac_3d(umac,rh,dx,rhmax,rhsum)

      real(kind=dp_t), intent(in   ) :: umac(-2:,-2:,-2:,:)
      real(kind=dp_t), intent(inout) ::   rh( 0:, 0:, 0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      real(kind=dp_t), intent(  out) :: rhsum,rhmax

      integer :: i,j,k

      rhsum = 0.0_dp_t
      rhmax = 0.0_dp_t

      do k = 0,size(rh,dim=3)-1
      do j = 0,size(rh,dim=2)-1
      do i = 0,size(rh,dim=1)-1
         rh(i,j,k) = (umac(i+1,j,k,1) - umac(i,j,k,1)) / dx(1) + &
                     (umac(i,j+1,k,2) - umac(i,j,k,2)) / dx(2) + &
                     (umac(i,j,k+1,3) - umac(i,j,k,3)) / dx(3)
         rhsum = rhsum + rh(i,j,k)
         rhmax = max(rhmax,abs(rh(i,j,k)))
         rh(i,j,k) = -rh(i,j,k)
      end do
      end do
      end do

    end subroutine divumac_3d

    subroutine mk_mac_coeffs(rho,beta)

      type(multifab) , intent(in   ) :: rho
      type(multifab) , intent(inout) :: beta
 
      real(kind=dp_t), pointer :: bp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer :: i,dm,ng

      dm = rho%dim
      ng = rho%ng

      do i = 1, rho%nboxes
         if ( multifab_remote(rho, i) ) cycle
         rp => dataptr(rho , i)
         bp => dataptr(beta, i)
         select case (dm)
            case (2)
              call mk_mac_coeffs_2d(bp(:,:,1,:), rp(:,:,1,1), ng)
            case (3)
              call mk_mac_coeffs_3d(bp(:,:,:,:), rp(:,:,:,1), ng)
         end select
      end do

    end subroutine mk_mac_coeffs

    subroutine mk_mac_coeffs_2d(beta,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: beta( -1:, -1:,:)
      real(kind=dp_t), intent(inout) ::  rho(-ng:,-ng:)

      integer :: i,j
      integer :: nx,ny
 
      nx = size(beta,dim=1) - 2
      ny = size(beta,dim=2) - 2

      do j = 0,ny-1
      do i = 0,nx
         beta(i,j,1) = TWO / (rho(i,j) + rho(i-1,j))
      end do
      end do

      do j = 0,ny
      do i = 0,nx-1
         beta(i,j,2) = TWO / (rho(i,j) + rho(i,j-1))
      end do
      end do

    end subroutine mk_mac_coeffs_2d

    subroutine mk_mac_coeffs_3d(beta,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: beta( -1:, -1:, -1:,:)
      real(kind=dp_t), intent(inout) ::  rho(-ng:,-ng:,-ng:)

      integer :: i,j,k
      integer :: nx,ny,nz
 
      nx = size(beta,dim=1) - 2
      ny = size(beta,dim=2) - 2
      nz = size(beta,dim=3) - 2

      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx
         beta(i,j,k,1) = TWO / (rho(i,j,k) + rho(i-1,j,k))
      end do
      end do
      end do

      do k = 0,nz-1
      do j = 0,ny
      do i = 0,nx-1
         beta(i,j,k,2) = TWO / (rho(i,j,k) + rho(i,j-1,k))
      end do
      end do
      end do

      do k = 0,nz
      do j = 0,ny-1
      do i = 0,nx-1
         beta(i,j,k,3) = TWO / (rho(i,j,k) + rho(i,j,k-1))
      end do
      end do
      end do

    end subroutine mk_mac_coeffs_3d

    subroutine mkumac(umac,phi,beta,dx,press_bc,verbose)

      type(multifab), intent(inout) :: umac
      type(multifab), intent(in   ) :: phi
      type(multifab), intent(in   ) :: beta
      real(dp_t) :: dx(:)
      type(bc_level), intent(in   ) :: press_bc
      integer       , intent(in   ) :: verbose

      integer :: i,dm
 
      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: php(:,:,:,:) 
      real(kind=dp_t), pointer ::  bp(:,:,:,:) 
      real(kind=dp_t)          :: rhmax,local_max

      dm = umac%dim
      rhmax = ZERO

      do i = 1, umac%nboxes
         if ( multifab_remote(umac, i) ) cycle
         ump => dataptr(umac, i)
         php => dataptr(phi , i)
          bp => dataptr(beta, i)
         select case (dm)
            case (2)
              call mkumac_2d(ump(:,:,1,:), php(:,:,1,1), bp(:,:,1,:), dx, &
                             press_bc%bc_level_array(i,:,:), local_max)
            case (3)
              call mkumac_3d(ump(:,:,:,:), php(:,:,:,1), bp(:,:,:,:), dx, &
                             press_bc%bc_level_array(i,:,:), local_max)
         end select
         rhmax = max(rhmax,local_max)
      end do
      if (verbose .eq. 1) print *,'MAX DIVU AFTER PROJECTION',rhmax

    end subroutine mkumac

    subroutine mkumac_2d(umac,phi,beta,dx,press_bc,rhmax)

      real(kind=dp_t), intent(inout) :: umac(-2:,-2:,:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)
      real(kind=dp_t), intent(  out) :: rhmax

      real(kind=dp_t) :: gpx,gpy
      real(kind=dp_t) :: rh
      integer :: i,j,nx,ny

      nx = size(umac,dim=1) - 4
      ny = size(umac,dim=2) - 4

      if (press_bc(1,1) == BC_NEU) then
         do j = 0,ny-1
            phi(-1,j) = phi(0,j)
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do j = 0,ny-1
            phi(-1,j) = -phi(0,j)
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do j = 0,ny-1
            phi(nx,j) = phi(nx-1,j)
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do j = 0,ny-1
            phi(nx,j) = -phi(nx-1,j)
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do i = 0,nx-1
            phi(i,-1) = phi(i,0)
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do i = 0,nx-1
            phi(i,-1) = -phi(i,0)
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do i = 0,nx-1
            phi(i,ny) = phi(i,ny-1)
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do i = 0,nx-1
            phi(i,ny) = -phi(i,ny-1)
         end do
      end if

      do j = 0,ny-1
      do i = 0,nx
         gpx = (phi(i,j) - phi(i-1,j)) / dx(1)
         umac(i,j,1) = umac(i,j,1) - beta(i,j,1)*gpx
      end do
      end do

      do j = 0,ny
      do i = 0,nx-1
         gpy = (phi(i,j) - phi(i,j-1)) / dx(2)
         umac(i,j,2) = umac(i,j,2) - beta(i,j,2)*gpy
      end do
      end do

!     This is just a test
      rhmax = 0.0_dp_t
      do j = 0,ny-1
      do i = 0,nx-1
         rh = (umac(i+1,j,1) - umac(i,j,1)) / dx(1) + &
              (umac(i,j+1,2) - umac(i,j,2)) / dx(2)
         rhmax = max(abs(rh),rhmax)
      end do
      end do

    end subroutine mkumac_2d

    subroutine mkumac_3d(umac,phi,beta,dx,press_bc,rhmax)

      real(kind=dp_t), intent(inout) :: umac(-2:,-2:,-2:,:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)
      real(kind=dp_t), intent(  out) :: rhmax

      real(kind=dp_t) :: rh
      real(kind=dp_t) :: gpx,gpy,gpz
      integer :: i,j,k,nx,ny,nz

      nx = size(umac,dim=1) - 4
      ny = size(umac,dim=2) - 4
      nz = size(umac,dim=3) - 4

      if (press_bc(1,1) == BC_NEU) then
         do k = 0,nz-1
         do j = 0,ny-1
            phi(-1,j,k) = phi(0,j,k)
         end do
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do k = 0,nz-1
         do j = 0,ny-1
            phi(-1,j,k) = -phi(0,j,k)
         end do
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do k = 0,nz-1
         do j = 0,ny-1
            phi(nx,j,k) = phi(nx-1,j,k)
         end do
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do k = 0,nz-1
         do j = 0,ny-1
            phi(nx,j,k) = -phi(nx-1,j,k)
         end do
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do k = 0,nz-1
         do i = 0,nx-1
            phi(i,-1,k) = phi(i,0,k)
         end do
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do k = 0,nz-1
         do i = 0,nx-1
            phi(i,-1,k) = -phi(i,0,k)
         end do
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do k = 0,nz-1
         do i = 0,nx-1
            phi(i,ny,k) = phi(i,ny-1,k)
         end do
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do k = 0,nz-1
         do i = 0,nx-1
            phi(i,ny,k) = -phi(i,ny-1,k)
         end do
         end do
      end if
      if (press_bc(3,1) == BC_NEU) then
         do j = 0,ny-1
         do i = 0,nx-1
            phi(i,j,-1) = phi(i,j,0)
         end do
         end do
      else if (press_bc(3,1) == BC_DIR) then
         do j = 0,ny-1
         do i = 0,nx-1
            phi(i,j,-1) = -phi(i,j,0)
         end do
         end do
      end if
      if (press_bc(3,2) == BC_NEU) then
         do j = 0,ny-1
         do i = 0,nx-1
            phi(i,j,nz) = phi(i,j,nz-1)
         end do
         end do
      else if (press_bc(3,2) == BC_DIR) then
         do j = 0,ny-1
         do i = 0,nx-1
            phi(i,j,nz) = -phi(i,j,nz-1)
         end do
         end do
      end if

      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx
         gpx = (phi(i,j,k) - phi(i-1,j,k)) / dx(1)
         umac(i,j,k,1) = umac(i,j,k,1) - beta(i,j,k,1)*gpx
      end do
      end do
      end do

      do k = 0,nz-1
      do j = 0,ny
      do i = 0,nx-1
         gpy = (phi(i,j,k) - phi(i,j-1,k)) / dx(2)
         umac(i,j,k,2) = umac(i,j,k,2) - beta(i,j,k,2)*gpy
      end do
      end do
      end do

      do k = 0,nz
      do j = 0,ny-1
      do i = 0,nx-1
         gpz = (phi(i,j,k) - phi(i,j,k-1)) / dx(3)
         umac(i,j,k,3) = umac(i,j,k,3) - beta(i,j,k,3)*gpz
      end do
      end do
      end do

!     This is just a test
      rhmax = 0.0_dp_t
      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx-1
         rh = (umac(i+1,j,k,1) - umac(i,j,k,1)) / dx(1) + &
              (umac(i,j+1,k,2) - umac(i,j,k,2)) / dx(2) + &
              (umac(i,j,k+1,3) - umac(i,j,k,3)) / dx(3)
         rhmax = max(abs(rh),rhmax)
      end do
      end do
      end do

    end subroutine mkumac_3d

end subroutine macproject

subroutine mac_multigrid(la,rh,phi,alpha,beta,dx,bc,stencil_order,mg_verbose)
  use BoxLib
  use omp_module
  use f2kcli
  use stencil_module
  use coeffs_module
  use mg_module
  use list_box_module
  use itsol_module
  use sparse_solve_module
  use bl_mem_stat_module
  use box_util_module
  use bl_IO_module

  type(layout),intent(inout) :: la
  integer     ,intent(in   ) :: stencil_order
  integer     ,intent(in   ) :: mg_verbose

  type(boxarray) pdv
  type(box) pd

  type(multifab), allocatable :: coeffs(:)

  real(dp_t), intent(in) :: dx(:)
  integer, intent(in   ) :: bc(:,:)

  type(multifab), intent(in   ) :: alpha, beta
  type(multifab), intent(inout) :: rh, phi
  type( multifab) :: ss
  type( multifab) :: phi0
  type(imultifab) :: mm
  type(sparse) :: sparse_object
  type(mg_tower) :: mgt
  integer i, dm, ns

  integer :: test
  real(dp_t) :: snrm(2)

  ! MG solver defaults
  integer :: bottom_solver, bottom_max_iter
  real(dp_t) :: bottom_solver_eps
  real(dp_t) :: eps
  integer :: max_iter
  integer :: min_width
  integer :: max_nlevel
  integer :: verbose
  integer :: nu1, nu2, gamma, cycle, solver, smoother
  real(dp_t) :: omega
  logical :: nodal(rh%dim)

  real(dp_t) :: xa(rh%dim), xb(rh%dim)

  !! Defaults:

  test           = 0

  dm             = rh%dim

  max_nlevel        = mgt%max_nlevel
  max_iter          = mgt%max_iter
  eps               = mgt%eps
  solver            = mgt%solver
  smoother          = mgt%smoother
  nu1               = mgt%nu1
  nu2               = mgt%nu2
  gamma             = mgt%gamma
  omega             = mgt%omega
  cycle             = mgt%cycle
  bottom_solver     = mgt%bottom_solver
  bottom_solver_eps = mgt%bottom_solver_eps
  bottom_max_iter   = mgt%bottom_max_iter
  min_width         = mgt%min_width
  verbose           = mgt%verbose

! Note: put this here to minimize asymmetries - ASA
  eps = 1.d-12

  bottom_solver = 0
  
  nodal = .FALSE.

  if ( test /= 0 .AND. max_iter == mgt%max_iter ) &
     max_iter = 1000

  ns = 1 + dm*3

  pd = bbox(get_boxarray(la))

  call mg_tower_build(mgt, la, pd, bc, &
       dh = dx(:), &
       ns = ns, &
       smoother = smoother, &
       nu1 = nu1, &
       nu2 = nu2, &
       gamma = gamma, &
       cycle = cycle, &
       omega = omega, &
       bottom_solver = bottom_solver, &
       bottom_max_iter = bottom_max_iter, &
       bottom_solver_eps = bottom_solver_eps, &
       max_iter = max_iter, &
       max_nlevel = max_nlevel, &
       min_width = min_width, &
       eps = eps, &
       verbose = verbose, &
       nodal = nodal)

  allocate(coeffs(mgt%nlevels))
  call multifab_build(coeffs(mgt%nlevels), la, 1+dm, 1)
  call multifab_copy_c(coeffs(mgt%nlevels),1,alpha,1,1,all=.true.)
  call multifab_copy_c(coeffs(mgt%nlevels),2,beta,1,dm,all=.true.)

! Do multigrid solve here 
  xa = 0.0_dp_t
  xb = 0.0_dp_t
  do i = mgt%nlevels-1, 1, -1
    call multifab_build(coeffs(i), mgt%ss(i)%la, 1+dm, 1)
    call setval(coeffs(i), 0.0_dp_t, 1, dm+1, all=.true.)
    call coarsen_coeffs(coeffs(i+1),coeffs(i))
  end do
  do i = mgt%nlevels, 1, -1
     pdv = layout_boxarray(mgt%ss(i)%la)
     call stencil_fill_cc(mgt%ss(i), coeffs(i), mgt%dh(:,i), pdv, &
          mgt%mm(i), xa, xb, pd, stencil_order, bc)
  end do

  if ( bottom_solver == 3 ) then
     call sparse_build(mgt%sparse_object, mgt%ss(1), mgt%mm(1), &
          mgt%ss(1)%la, stencil_order, mgt%verbose)
  end if

! Put the problem in residual-correction form.
  i = mgt%nlevels
  call multifab_build(phi0,la,1,1)
  call multifab_copy_c(phi0,1,phi,1,1,all=.true.)
  call mg_defect(mgt%ss(i),mgt%dd(i),rh,phi0,mgt%mm(i))
  call multifab_copy_c(rh,1,mgt%dd(i),1,1,all=.true.)

  call setval(phi,ZERO,all=.true.)
  call mg_tower_solve(mgt, phi, rh)
  call saxpy(phi,1.d0,phi0)

  do i = mgt%nlevels-1, 1, -1
    call multifab_destroy(coeffs(i))
  end do

  call multifab_fill_boundary(phi)

  call multifab_destroy(phi0)

  call multifab_destroy(coeffs(mgt%nlevels))
  deallocate(coeffs)

  snrm(1) = norm_l2(phi)
  snrm(2) = norm_inf(phi)
  if ( mg_verbose > 0 .and. parallel_IOProcessor() ) &
       print *, 'solution norm = ', snrm(1), "/", snrm(2)

  if ( test == 3 ) then
     call sparse_destroy(sparse_object)
  end if
  if ( test > 0 ) then
     call destroy(ss)
     call destroy(mm)
  end if
  call mg_tower_destroy(mgt)

end subroutine mac_multigrid

end module macproject_module
