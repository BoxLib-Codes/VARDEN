module hgproject_module

  use bl_types
  use bc_module
  use define_bc_module
  use multifab_module
  use boxarray_module
  use stencil_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO   = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE    = 1.0_dp_t
  real (kind = dp_t), private, parameter :: TWO    = 2.0_dp_t
  real (kind = dp_t), private, parameter :: HALF   = 0.5_dp_t
  real (kind = dp_t), private, parameter :: FOURTH = 0.25_dp_t
  real (kind = dp_t), private, parameter :: EIGHTH = 0.125_dp_t

contains 

subroutine hgproject(unew,rhohalf,p,gp,dx,dt,phys_bc,press_bc,mg_verbose)

  type(multifab), intent(inout) :: unew
  type(multifab), intent(inout) :: rhohalf
  type(multifab), intent(inout) :: gp
  type(multifab), intent(inout) :: p
  type(bc_level), intent(in   ) :: phys_bc
  integer       , intent(in   ) :: press_bc(:,:)
  real(dp_t) :: dx(:),dt
  integer       , intent(in   ) :: mg_verbose

! Local  
  type(multifab) :: rh,phi,gphi
  type(layout) :: la
  logical   , allocatable :: nodal(:)
  integer  :: dm,ng,i
  type(box) :: bx

  la = unew%la
  dm = unew%dim
  ng = unew%ng
  allocate(nodal(dm))
  nodal = .True.

  call multifab_build(  rh, la, 1, 1, nodal)
  call multifab_build( phi, la, 1, 1, nodal)
  call multifab_build(gphi, la, dm, 0) 

  print *,' '
  print *,'... begin hg_projection ... '

  call divu(unew,rhohalf,gp,rh,dx,dt,phys_bc)

  call hg_multigrid(la,rh,rhohalf,phi,dx,press_bc,mg_verbose)

  call mkgp(gphi,phi,dx)
  call mkunew(unew,gp,gphi,rhohalf,ng,dt)

  call multifab_copy_c(p,1,phi,1,1,all=.true.)

  print *,'...   end hg_projection ... '
  print *,' '

  call multifab_destroy(rh)
  call multifab_destroy(phi)
  call multifab_destroy(gphi)

  contains

    subroutine divu(unew,rhohalf,gp,rh,dx,dt,phys_bc)

      type(multifab) , intent(in   ) :: unew
      type(multifab) , intent(in   ) :: gp
      type(multifab) , intent(in   ) :: rhohalf
      type(multifab) , intent(inout) :: rh
      real(kind=dp_t), intent(in   ) :: dx(:),dt
      type(bc_level) , intent(in   ) :: phys_bc
 
      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 

      integer :: i,dm,ng
      real(kind=dp_t)          :: rhmax,rhsum
      real(kind=dp_t)          :: local_max,local_sum

      dm = rh%dim
      ng = unew%ng

      rhmax = ZERO
      rhsum = ZERO

      do i = 1, unew%nboxes
         if ( multifab_remote(unew, i) ) cycle
         ump => dataptr(unew   , i)
         gpp => dataptr(gp     , i)
          rp => dataptr(rhohalf, i)
         rhp => dataptr(rh     , i)
         select case (dm)
            case (2)
              call divu_2d(ump(:,:,1,:), rp(:,:,1,1), gpp(:,:,1,:), &
                           rhp(:,:,1,1), dx, dt, phys_bc%bc_level_array(i,:,:), ng, &
                           local_max, local_sum)
            case (3)
              call divu_3d(ump(:,:,:,:), rp(:,:,:,1), gpp(:,:,:,:), &
                           rhp(:,:,:,1), dx, dt, phys_bc%bc_level_array(i,:,:), ng, &
                           local_max, local_sum)
         end select
         rhmax = max(rhmax,local_max)
         rhsum = rhsum + local_sum
      end do
      print *,'HGPROJ RHS MAX / SUM',rhmax,rhsum

    end subroutine divu

    subroutine divu_2d(u,rhohalf,gp,rh,dx,dt,phys_bc,ng,rhmax,rhsum)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(in   ) :: gp(-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:),dt
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(  out) :: rhmax,rhsum

      integer :: i,j,nx,ny
      nx = size(rh,dim=1) - ng
      ny = size(rh,dim=2) - ng

      do j = -1,ny
      do i = -1,nx
         u(i,j,1) = u(i,j,1)/dt + gp(i,j,1)/rhohalf(i,j)
         u(i,j,2) = u(i,j,2)/dt + gp(i,j,2)/rhohalf(i,j)
      end do
      end do

      if (phys_bc(1,1) == WALL) u(-1,:,:) = ZERO
      if (phys_bc(1,2) == WALL) u(nx,:,:) = ZERO
      if (phys_bc(2,1) == WALL) u(:,-1,:) = ZERO
      if (phys_bc(2,2) == WALL) u(:,ny,:) = ZERO

      do j = 0,ny
      do i = 0,nx
         rh(i,j) = HALF * (u(i  ,j,1) + u(i  ,j-1,1) &
                          -u(i-1,j,1) - u(i-1,j-1,1)) / dx(1) + &
                   HALF * (u(i,j  ,2) + u(i-1,j  ,2) &
                          -u(i,j-1,2) - u(i-1,j-1,2)) / dx(2)
      end do
      end do

      if (phys_bc(1,1) == OUTLET) rh( 0,:) = ZERO
      if (phys_bc(1,2) == OUTLET) rh(nx,:) = ZERO
      if (phys_bc(2,1) == OUTLET) rh(:, 0) = ZERO
      if (phys_bc(2,2) == OUTLET) rh(:,ny) = ZERO

      rhmax = ZERO
      rhsum = ZERO
      do j = 0,ny
      do i = 0,nx
         rhmax = max(rhmax,abs(rh(i,j)))
         rhsum = rhsum + rh(i,j)
      end do
      end do

      if ( (phys_bc(1,1) == WALL .or. phys_bc(1,1) == INLET) .and. &
           (phys_bc(1,2) == WALL .or. phys_bc(1,2) == INLET) .and. &
           (phys_bc(2,1) == WALL .or. phys_bc(2,1) == INLET) .and. &
           (phys_bc(2,2) == WALL .or. phys_bc(2,2) == INLET) ) then
         do j = 0,ny
         do i = 0,nx
            rh(i,j) = rh(i,j) - rhsum / dble(nx+1) / dble(ny+1)
         end do
         end do
      end if

      if (phys_bc(1,1) == WALL .or. phys_bc(1,1) == INLET) rh( 0,:) = TWO*rh( 0,:)
      if (phys_bc(1,2) == WALL .or. phys_bc(1,2) == INLET) rh(nx,:) = TWO*rh(nx,:)
      if (phys_bc(2,1) == WALL .or. phys_bc(2,1) == INLET) rh(:, 0) = TWO*rh(:, 0)
      if (phys_bc(2,2) == WALL .or. phys_bc(2,2) == INLET) rh(:,ny) = TWO*rh(:,ny)

    end subroutine divu_2d

    subroutine divu_3d(u,rhohalf,gp,rh,dx,dt,phys_bc,ng,rhmax,rhsum)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(in   ) :: gp(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:),dt
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(  out) :: rhmax,rhsum

      integer :: i,j,k,nx,ny,nz

      nx = size(rh,dim=1) - ng
      ny = size(rh,dim=2) - ng
      nz = size(rh,dim=3) - ng

      if (phys_bc(1,1) == WALL) u(-1,:,:,:) = ZERO
      if (phys_bc(1,2) == WALL) u(nx,:,:,:) = ZERO
      if (phys_bc(2,1) == WALL) u(:,-1,:,:) = ZERO
      if (phys_bc(2,2) == WALL) u(:,ny,:,:) = ZERO
      if (phys_bc(3,1) == WALL) u(:,:,-1,:) = ZERO
      if (phys_bc(3,2) == WALL) u(:,:,nz,:) = ZERO

      do k = -1,nz
      do j = -1,ny
      do i = -1,nx
         u(i,j,k,1) = u(i,j,k,1)/dt + gp(i,j,k,1)/rhohalf(i,j,k)
         u(i,j,k,2) = u(i,j,k,2)/dt + gp(i,j,k,2)/rhohalf(i,j,k)
         u(i,j,k,3) = u(i,j,k,3)/dt + gp(i,j,k,3)/rhohalf(i,j,k)
      end do
      end do
      end do

      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         rh(i,j,k) = (u(i  ,j,k  ,1) + u(i  ,j-1,k  ,1) &
                     +u(i  ,j,k-1,1) + u(i  ,j-1,k-1,1) &
                     -u(i-1,j,k  ,1) - u(i-1,j-1,k  ,1) &
                     -u(i-1,j,k-1,1) - u(i-1,j-1,k-1,1)) / dx(1) + &
                     (u(i,j  ,k  ,2) + u(i-1,j  ,k  ,2) &
                     +u(i,j  ,k-1,2) + u(i-1,j  ,k-1,2) &
                     -u(i,j-1,k  ,2) - u(i-1,j-1,k  ,2) &
                     -u(i,j-1,k-1,2) - u(i-1,j-1,k-1,2)) / dx(2) + &
                     (u(i,j  ,k  ,3) + u(i-1,j  ,k  ,3) &
                     +u(i,j-1,k  ,3) + u(i-1,j-1,k  ,3) &
                     -u(i,j  ,k-1,3) - u(i-1,j  ,k-1,3) &
                     -u(i,j-1,k-1,3) - u(i-1,j-1,k-1,3)) / dx(3)
         rh(i,j,k) = FOURTH*rh(i,j,k)
      end do
      end do
      end do

      if (phys_bc(1,1) == OUTLET) rh( 0,:,:) = ZERO
      if (phys_bc(1,2) == OUTLET) rh(nx,:,:) = ZERO
      if (phys_bc(2,1) == OUTLET) rh(:, 0,:) = ZERO
      if (phys_bc(2,2) == OUTLET) rh(:,ny,:) = ZERO
      if (phys_bc(3,1) == OUTLET) rh(:,:, 0) = ZERO
      if (phys_bc(3,2) == OUTLET) rh(:,:,nz) = ZERO

      rhsum = ZERO
      rhmax = ZERO
      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         rhsum = rhsum + rh(i,j,k)
         rhmax = max(rhmax,abs(rh(i,j,k)))
      end do
      end do
      end do

      if (phys_bc(1,1) == WALL .or. phys_bc(1,1) == INLET) rh( 0,:,:) = TWO*rh( 0,:,:)
      if (phys_bc(1,2) == WALL .or. phys_bc(1,2) == INLET) rh(nx,:,:) = TWO*rh(nx,:,:)
      if (phys_bc(2,1) == WALL .or. phys_bc(2,1) == INLET) rh(:, 0,:) = TWO*rh(:, 0,:)
      if (phys_bc(2,2) == WALL .or. phys_bc(2,2) == INLET) rh(:,ny,:) = TWO*rh(:,ny,:)
      if (phys_bc(3,1) == WALL .or. phys_bc(3,1) == INLET) rh(:,:, 0) = TWO*rh(:,:, 0)
      if (phys_bc(3,2) == WALL .or. phys_bc(3,2) == INLET) rh(:,:,nz) = TWO*rh(:,:,nz)

    end subroutine divu_3d

    subroutine mkgp(gphi,phi,dx)

      type(multifab), intent(inout) :: gphi
      type(multifab), intent(in   ) :: phi
      real(dp_t) :: dx(:)

      integer :: i,dm
 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      dm = phi%dim

      do i = 1, phi%nboxes
         if ( multifab_remote(phi, i) ) cycle
         gph => dataptr(gphi, i)
         pp  => dataptr(phi , i)
         select case (dm)
            case (2)
              call mkgp_2d(gph(:,:,1,:), pp(:,:,1,1), dx)
            case (3)
              call mkgp_3d(gph(:,:,:,:), pp(:,:,:,1), dx)
         end select
      end do

      call multifab_fill_boundary(gphi)

    end subroutine mkgp

    subroutine mkgp_2d(gp,phi,dx)

      real(kind=dp_t), intent(inout) ::  gp(0:,0:,:)
      real(kind=dp_t), intent(inout) :: phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,nx,ny

      nx = size(gp,dim=1)
      ny = size(gp,dim=2)

      do j = 0,ny-1
      do i = 0,nx-1
          gp(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                            phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
          gp(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                            phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
      end do
      end do

    end subroutine mkgp_2d

    subroutine mkgp_3d(gp,phi,dx)

      real(kind=dp_t), intent(inout) ::  gp(0:,0:,0:,1:)
      real(kind=dp_t), intent(inout) :: phi(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k,nx,ny,nz

      nx = size(gp,dim=1)
      ny = size(gp,dim=2)
      nz = size(gp,dim=3)

      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx-1
         gp(i,j,k,1) = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                              +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                              -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                              -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
         gp(i,j,k,2) = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                              +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                              -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                              -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
         gp(i,j,k,3) = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                              +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                              -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                              -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)
      end do
      end do
      end do

    end subroutine mkgp_3d

    subroutine mkunew(unew,gp,gphi,rhohalf,ng,dt)

      type(multifab) , intent(inout) :: unew
      type(multifab) , intent(inout) :: gp
      type(multifab) , intent(in   ) :: gphi
      type(multifab) , intent(in   ) :: rhohalf
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(in   ) :: dt 
      integer :: i,dm
 
      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 

      dm = unew%dim

      do i = 1, unew%nboxes
         if ( multifab_remote(unew, i) ) cycle
         upn => dataptr(unew, i)
         gpp => dataptr(gp  , i)
         gph => dataptr(gphi, i)
         rp  => dataptr(rhohalf, i)
         select case (dm)
            case (2)
              call mkunew_2d(upn(:,:,1,:), gpp(:,:,1,:), gph(:,:,1,:),rp(:,:,1,1),ng,dt)
            case (3)
              call mkunew_3d(upn(:,:,:,:), gpp(:,:,:,:), gph(:,:,:,:),rp(:,:,:,1),ng,dt)
         end select
      end do
      call multifab_fill_boundary(unew)
      call multifab_fill_boundary(gp)

    end subroutine mkunew

    subroutine mkunew_2d(unew,gp,gphi,rhohalf,ng,dt)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) :: unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp(-1:,-1:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(0:,0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer :: i,j,nx,ny

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1

!     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,1) = unew(0:nx,0:ny,1) - gphi(0:nx,0:ny,1)/rhohalf(0:nx,0:ny) 
      unew(0:nx,0:ny,2) = unew(0:nx,0:ny,2) - gphi(0:nx,0:ny,2)/rhohalf(0:nx,0:ny) 

!     Multiply by dt.
      unew(-1:nx,-1:ny,:) = dt * unew(-1:nx,-1:ny,:)
 
!     Replace gradient of pressure by gradient of phi.
      gp(0:nx,0:ny,:) = gphi(0:nx,0:ny,:)

    end subroutine mkunew_2d

    subroutine mkunew_3d(unew,gp,gphi,rhohalf,ng,dt)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) :: unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::   gp(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: gphi( 0:, 0:, 0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer :: nx,ny,nz

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1
      nz = size(gphi,dim=3)-1

!     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,0:nz,1) = &
          unew(0:nx,0:ny,0:nz,1) - gphi(0:nx,0:ny,0:nz,1)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,2) = &
          unew(0:nx,0:ny,0:nz,2) - gphi(0:nx,0:ny,0:nz,2)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,3) = &
          unew(0:nx,0:ny,0:nz,3) - gphi(0:nx,0:ny,0:nz,3)/rhohalf(0:nx,0:ny,0:nz) 

!     Multiply by dt.
      unew(-1:nx,-1:ny,-1:nz,:) = dt * unew(-1:nx,-1:ny,-1:nz,:)

!     Replace (gradient of) pressure by (gradient of) phi.
      gp(0:nx,0:ny,0:nz,:) = gphi(0:nx,0:ny,0:nz,:)

    end subroutine mkunew_3d

end subroutine hgproject

subroutine hg_multigrid(la,rh,rhohalf,phi,dx,press_bc,mg_verbose)
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
  use bl_timer_module
  use box_util_module
  use bl_IO_module

  type(layout),intent(inout) :: la

  type(box     ) pd
  type(boxarray) pdv

  type(multifab), allocatable :: coeffs(:)

  real(dp_t), intent(in) :: dx(:)
  integer, intent(in   ) :: press_bc(:,:)
  integer, intent(in   ) :: mg_verbose

  type(multifab), intent(inout) :: phi, rh
  type(multifab), intent(inout) :: rhohalf
  type( multifab) :: ss
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
  integer :: ng, nc

  type(timer) :: tm(2)

  integer :: stencil_order
  logical :: nodal(rh%dim)

  real(dp_t) :: xa(rh%dim), xb(rh%dim)

  !! Defaults:


  test           = 0

  dm                = rh%dim

  ng                = mgt%ng
  nc                = mgt%nc
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

  bottom_solver = 2
  
  nodal = .TRUE.

  if ( test /= 0 .AND. max_iter == mgt%max_iter ) then
     max_iter = 1000
  end if

  call setval(phi, 0.0_dp_t, all=.true.)

  ns = 3**dm
 
  pd = bbox(get_boxarray(la))

  eps = 1.d-13

  call mg_tower_build(mgt, la, pd, press_bc, &
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
       verbose = mg_verbose, &
       nodal = nodal)

! Note the minus sign here is because we solve (-del dot b grad) phi = rhs,
!   NOT (del dot b grad) phi = rhs
  allocate(coeffs(mgt%nlevels))
  call multifab_build(coeffs(mgt%nlevels), la, 1, 1)
  call setval(coeffs(mgt%nlevels), 0.0_dp_t, 1, all=.true.)

  call mkcoeffs(rhohalf,coeffs(mgt%nlevels))
  call multifab_fill_boundary(coeffs(mgt%nlevels))

! Do multigrid solve here 
  xa = 0.0_dp_t
  xb = 0.0_dp_t
  call timer_start(tm(1))
  do i = mgt%nlevels-1, 1, -1
    call multifab_build(coeffs(i), mgt%ss(i)%la, 1, 1)
    call setval(coeffs(i), 0.0_dp_t, 1, all=.true.)
    call coarsen_coeffs(coeffs(i+1),coeffs(i))
    call multifab_fill_boundary(coeffs(i))
  end do

  do i = mgt%nlevels, 1, -1
     pdv = layout_boxarray(mgt%ss(i)%la)
     call stencil_fill_nodal(mgt%ss(i), coeffs(i), mgt%dh(:,i), &
          mgt%dh(:,mgt%nlevels), &
          mgt%mm(i), mgt%face_type, pd, pdv)
     pd  = coarsen(pd,2)
  end do

  if ( bottom_solver == 3 ) then
     call sparse_build(mgt%sparse_object, mgt%ss(1), mgt%mm(1), &
          mgt%ss(1)%la, stencil_order, mgt%verbose)
  end if
  call timer_stop(tm(1))
  call timer_start(tm(2))
  call mg_tower_solve(mgt, phi, rh)
  call timer_stop(tm(2))
  do i = mgt%nlevels-1, 1, -1
    call multifab_destroy(coeffs(i))
  end do

  call multifab_fill_boundary(phi)

  call multifab_destroy(coeffs(mgt%nlevels))
  deallocate(coeffs)

  snrm(1) = norm_l2(phi)
  snrm(2) = norm_inf(phi)
  if ( parallel_IOProcessor() ) &
       print *, 'solution norm = ', snrm(1), "/", snrm(2)

  if ( test == 3 ) then
     call sparse_destroy(sparse_object)
  end if
  if ( test > 0 ) then
     call destroy(ss)
     call destroy(mm)
  end if
  call mg_tower_destroy(mgt)

end subroutine hg_multigrid

    subroutine mkcoeffs(rho,coeffs)

      type(multifab) , intent(in   ) :: rho
      type(multifab) , intent(inout) :: coeffs

      real(kind=dp_t), pointer :: cp(:,:,:,:)
      real(kind=dp_t), pointer :: rp(:,:,:,:)
      integer :: i,dm,ng

      dm = rho%dim
      ng = rho%ng

      do i = 1, rho%nboxes
         if ( multifab_remote(rho, i) ) cycle
         rp => dataptr(rho   , i)
         cp => dataptr(coeffs, i)
         select case (dm)
            case (2)
              call mkcoeffs_2d(cp(:,:,1,1), rp(:,:,1,1), ng)
            case (3)
              call mkcoeffs_3d(cp(:,:,:,1), rp(:,:,:,1), ng)
         end select
      end do

    end subroutine mkcoeffs

    subroutine mkcoeffs_2d(coeffs,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:)

      integer :: i,j
      integer :: nx,ny

      nx = size(coeffs,dim=1) - 2
      ny = size(coeffs,dim=2) - 2

      do j = 1,ny
      do i = 1,nx
         coeffs(i,j) = -ONE / rho(i,j)
      end do
      end do

    end subroutine mkcoeffs_2d

    subroutine mkcoeffs_3d(coeffs,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:, 0:)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:,1-ng:)

      integer :: i,j,k
      integer :: nx,ny,nz

      nx = size(coeffs,dim=1) - 2
      ny = size(coeffs,dim=2) - 2
      nz = size(coeffs,dim=3) - 2

      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
         coeffs(i,j,k) = -ONE / rho(i,j,k)
      end do
      end do
      end do

    end subroutine mkcoeffs_3d

end module hgproject_module
