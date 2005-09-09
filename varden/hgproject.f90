module hgproject_module

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use boxarray_module
  use init_module
  use stencil_module
  use ml_solve_module
  use ml_restriction_module
  use fabio_module

  implicit none

contains 

subroutine hgproject(mla,unew,rhohalf,p,gp,dx,dt,the_bc_tower, &
                     verbose,mg_verbose,press_comp)

  type(ml_layout), intent(inout) :: mla
  type(multifab ), intent(inout) :: unew(:)
  type(multifab ), intent(inout) :: rhohalf(:)
  type(multifab ), intent(inout) :: gp(:)
  type(multifab ), intent(inout) :: p(:)
  type(bc_tower ), intent(in   ) :: the_bc_tower
  integer        , intent(in   ) :: verbose,mg_verbose
  real(dp_t)     , intent(in   ) :: dx(:,:),dt
  integer        , intent(in   ) :: press_comp

! Local  
  type(multifab), allocatable :: phi(:),gphi(:)
  type(box)      :: pd,bx
  type(layout  ) :: la
  integer        :: n,nlevs,dm,ng,i,RHOCOMP
  integer        :: ng_for_fill
  logical, allocatable :: nodal(:)
  real(dp_t)     :: nrm

  nlevs = mla%nlevel
  dm    = mla%dim
  ng = unew(nlevs)%ng

  allocate(phi(nlevs), gphi(nlevs))
  allocate(nodal(dm))
  nodal = .true.

  do n = 1, nlevs
     call multifab_build( phi(n), mla%la(n), 1, 1, nodal)
     call multifab_build(gphi(n), mla%la(n), dm, 0) 
     call multifab_copy(phi(n),p(n))
     call multifab_mult_mult_s(phi(n),dt,all=.true.)
  end do

  do n = 1, nlevs
  end do

  if (verbose .eq. 1) then
     print *,' '
     nrm = ZERO
     do n = 1, nlevs
        nrm = max(nrm,norm_inf(unew(n)))
     end do
     print *,'... hgproject: max of u before projection ',nrm
  end if

  call create_uvec_for_projection(nlevs,unew,rhohalf,gp,dt,the_bc_tower)

  do n = 1, nlevs
     call setval(phi(n),ZERO,all=.true.)
  end do

  call hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower, &
                    verbose,mg_verbose,press_comp)

  do n = 1,nlevs
     call mkgp(gphi(n),phi(n),dx(n,:))
     call mk_p(   p(n),phi(n),dt)
     call mkunew(unew(n),gp(n),gphi(n),rhohalf(n),ng,dt)
  end do

  do n = nlevs,2,-1
     call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:))
     call ml_cc_restriction(  gp(n-1),  gp(n),mla%mba%rr(n-1,:))
  end do

  if (verbose .eq. 1) then
     nrm = ZERO
     do n = 1, nlevs
        nrm = max(nrm,norm_inf(unew(n)))
     end do
     print *,'... hgproject: max of u after projection ',nrm
     print *,' '
  end if

  do n = 1,nlevs
     call multifab_destroy(phi(n))
     call multifab_destroy(gphi(n))
  end do

  deallocate(phi)
  deallocate(gphi)

  contains

    subroutine create_uvec_for_projection(nlevs,unew,rhohalf,gp,dt,the_bc_tower)

      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(inout) :: gp(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_tower) , intent(in   ) :: the_bc_tower
 
      type(bc_level) :: bc

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 

      integer :: i,n,dm,ng

      dm = unew(nlevs)%dim
      ng = unew(nlevs)%ng

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n)   , i)
            gpp => dataptr(gp(n)     , i)
             rp => dataptr(rhohalf(n), i)
            select case (dm)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), rp(:,:,1,1), gpp(:,:,1,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), rp(:,:,:,1), gpp(:,:,:,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do

    end subroutine create_uvec_for_projection

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

    subroutine mk_p(p,phi,dt)

      type(multifab), intent(inout) :: p
      type(multifab), intent(in   ) :: phi
      real(kind=dp_t),intent(in   ) :: dt
 
      real(kind=dp_t), pointer :: ph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 
      integer                  :: i
      real(kind=dp_t)          :: dt_inv

      dt_inv = ONE/dt

      do i = 1, phi%nboxes
         if ( multifab_remote(phi, i) ) cycle
         ph => dataptr(phi, i)
         pp => dataptr(p  , i)
         pp = dt_inv * ph
      end do

      call multifab_fill_boundary(p)

    end subroutine mk_p

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

    subroutine create_uvec_2d(u,rhohalf,gp,dt,phys_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::       u(-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: i,j,nx,ny
      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,:) = ZERO

      u(-1:nx,-1:ny,1) = u(-1:nx,-1:ny,1) + dt*gp(-1:nx,-1:ny,1)/rhohalf(-1:nx,-1:ny)
      u(-1:nx,-1:ny,2) = u(-1:nx,-1:ny,2) + dt*gp(-1:nx,-1:ny,2)/rhohalf(-1:nx,-1:ny)

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) u(-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) u(nx,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) u(:,-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) u(:,ny,:) = ZERO

    end subroutine create_uvec_2d

!   ********************************************************************************************* !

    subroutine create_uvec_3d(u,rhohalf,gp,dt,phys_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::       u(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)

      integer :: i,j,k,nx,ny,nz

      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2
      nz = size(gp,dim=3) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,-1:nz,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,-1:nz,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gp(-1:nx,-1:ny,-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gp(-1:nx,-1:ny,nz,:) = ZERO

      u(-1:nx,-1:ny,-1:nz,1) = u(-1:nx,-1:ny,-1:nz,1) + &
                           dt*gp(-1:nx,-1:ny,-1:nz,1)/rhohalf(-1:nx,-1:ny,-1:nz)
      u(-1:nx,-1:ny,-1:nz,2) = u(-1:nx,-1:ny,-1:nz,2) + &
                           dt*gp(-1:nx,-1:ny,-1:nz,2)/rhohalf(-1:nx,-1:ny,-1:nz)
      u(-1:nx,-1:ny,-1:nz,3) = u(-1:nx,-1:ny,-1:nz,3) + &
                           dt*gp(-1:nx,-1:ny,-1:nz,3)/rhohalf(-1:nx,-1:ny,-1:nz)

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) u(-1,:,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) u(nx,:,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) u(:,-1,:,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) u(:,ny,:,:) = ZERO
      if (phys_bc(3,1)==SLIP_WALL .or. phys_bc(3,1)==NO_SLIP_WALL) u(:,:,-1,:) = ZERO
      if (phys_bc(3,2)==SLIP_WALL .or. phys_bc(3,2)==NO_SLIP_WALL) u(:,:,nz,:) = ZERO

    end subroutine create_uvec_3d

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

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

!     Replace (gradient of) pressure by 1/dt * (gradient of) phi.
      gp(0:nx,0:ny,:) = (ONE/dt) * gphi(0:nx,0:ny,:)

    end subroutine mkunew_2d

!   ********************************************************************************************* !

    subroutine mkunew_3d(unew,gp,gphi,rhohalf,ng,dt)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) :: unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::   gp(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: gphi( 0:, 0:, 0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer :: nx,ny,nz,i,j,k

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

!     Replace (gradient of) pressure by 1/dt * (gradient of) phi.
      gp(0:nx,0:ny,0:nz,:) = (ONE/dt) * gphi(0:nx,0:ny,0:nz,:)

    end subroutine mkunew_3d

end subroutine hgproject

subroutine hg_multigrid(mla,unew,rhohalf,phi,dx,the_bc_tower,&
                        divu_verbose,mg_verbose,press_comp)
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

  type(ml_layout), intent(inout) :: mla
  type(multifab ), intent(inout) :: unew(:)
  type(multifab ), intent(in   ) :: rhohalf(:)
  type(multifab ), intent(inout) :: phi(:)
  real(dp_t)     , intent(in)    :: dx(:,:)
  type(bc_tower ), intent(in   ) :: the_bc_tower
  integer        , intent(in   ) :: divu_verbose,mg_verbose
  integer        , intent(in   ) :: press_comp
!
! Local variables
!

  type( multifab) :: ss
  type(imultifab) :: mm
  type(sparse)    :: sparse_object
  type(layout)    :: la
  type(box     )  :: pd
  type(boxarray)  :: pdv

  type(mg_tower), allocatable :: mgt(:)
  type(multifab), allocatable :: coeffs(:),rh(:)

  real(dp_t) :: bottom_solver_eps
  real(dp_t) :: eps
  real(dp_t) :: omega

  integer :: i, dm, nlevs, ns, test
  integer :: bottom_solver, bottom_max_iter
  integer :: max_iter
  integer :: min_width
  integer :: max_nlevel
  integer :: nu1, nu2, gamma, cycle, solver, smoother
  integer :: n, ng, nc
  integer :: max_nlevel_in
  integer :: verbose
  integer :: do_diagnostics

  logical, allocatable :: nodal(:)

  !! Defaults:

  dm    = mla%dim
  nlevs = mla%nlevel

  allocate(nodal(dm))
  nodal = .true.

  allocate(mgt(nlevs))

  test           = 0

  ng                = mgt(nlevs)%ng
  nc                = mgt(nlevs)%nc
  max_nlevel        = mgt(nlevs)%max_nlevel
  max_iter          = mgt(nlevs)%max_iter
  eps               = mgt(nlevs)%eps
  solver            = mgt(nlevs)%solver
  smoother          = mgt(nlevs)%smoother
  nu1               = mgt(nlevs)%nu1
  nu2               = mgt(nlevs)%nu2
  gamma             = mgt(nlevs)%gamma
  omega             = mgt(nlevs)%omega
  cycle             = mgt(nlevs)%cycle
  bottom_solver     = mgt(nlevs)%bottom_solver
  bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
  bottom_max_iter   = mgt(nlevs)%bottom_max_iter
  min_width         = mgt(nlevs)%min_width
  verbose           = mgt(nlevs)%verbose

! Note: put this here to minimize asymmetries - ASA
  eps = 1.d-14

  bottom_solver = 2

  if ( test /= 0 .AND. max_iter == mgt(nlevs)%max_iter ) then
     max_iter = 1000
  end if

! Note: put this here for robustness
  max_iter = 100

  ns = 3**dm

  do n = nlevs, 1, -1

     if (n == 1) then
        max_nlevel_in = max_nlevel
     else
        if ( all(mla%mba%rr(n-1,:) == 2) ) then
           max_nlevel_in = 1
        else if ( all(mla%mba%rr(n-1,:) == 4) ) then
           max_nlevel_in = 2
        else 
           call bl_error("HG_MULTIGRID: confused about ref_ratio")
        end if
     end if

     pd = layout_get_pd(mla%la(n))

     call mg_tower_build(mgt(n), mla%la(n), pd, &
          the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
          dh = dx(n,:), &
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
          max_nlevel = max_nlevel_in, &
          min_width = min_width, &
          eps = eps, &
          verbose = verbose, &
          nodal = nodal)

  end do

  do n = nlevs,1,-1

     allocate(coeffs(mgt(n)%nlevels))

     la = mla%la(n)
     pd = layout_get_pd(la)

     call multifab_build(coeffs(mgt(n)%nlevels), la, 1, 1)
     call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)

     call mkcoeffs(rhohalf(n),coeffs(mgt(n)%nlevels))
     call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

     do i = mgt(n)%nlevels-1, 1, -1
        call multifab_build(coeffs(i), mgt(n)%ss(i)%la, 1, 1)
        call setval(coeffs(i), 0.0_dp_t, 1, all=.true.)
        call coarsen_coeffs(coeffs(i+1),coeffs(i))
        call multifab_fill_boundary(coeffs(i))
     end do

!    NOTE: we define the stencils with the finest dx.
     do i = mgt(n)%nlevels, 1, -1
        pdv = layout_boxarray(mgt(n)%ss(i)%la)
        call stencil_fill_nodal(mgt(n)%ss(i), coeffs(i), &
             mgt(n)%dh(:,i)             , &
             mgt(nlevs)%dh(:,mgt(nlevs)%nlevels), &
             mgt(n)%mm(i), mgt(n)%face_type)
        pd  = coarsen(pd,2)
     end do

     do i = mgt(n)%nlevels, 1, -1
        call multifab_destroy(coeffs(i))
     end do
     deallocate(coeffs)

  end do

  allocate(rh(nlevs))
  do n = 1, nlevs
     call multifab_build(rh(n),mla%la(n),1,1,nodal)
  end do

  call divu(nlevs,mgt,unew,rh,dx,the_bc_tower,mla%mba%rr,divu_verbose,press_comp,nodal)

  if ( mg_verbose >=1 ) then
    do_diagnostics = 1
  else
    do_diagnostics = 0
  end if

! call fabio_ml_write(rh, mla%mba%rr(:,1), "rh_nodal")

  call ml_nd_solve(mla,mgt,rh,phi,mla%mba%rr,do_diagnostics)

! call fabio_ml_write(phi, mla%mba%rr(:,1), "soln_nodal")

  do n = nlevs,1,-1
     call multifab_fill_boundary(phi(n))
  end do

  if ( test == 3 ) then
     call sparse_destroy(sparse_object)
  end if
  if ( test > 0 ) then
     call destroy(ss)
     call destroy(mm)
  end if

  do n = 1, nlevs
     call mg_tower_destroy(mgt(n))
     call multifab_destroy(rh(n))
  end do

  deallocate(mgt)
  deallocate(rh)

end subroutine hg_multigrid

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

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

!   ********************************************************************************************* !

    subroutine divu(nlevs,mgt,unew,rh,dx,the_bc_tower,ref_ratio,verbose,press_comp,nodal)

      integer        , intent(in   ) :: nlevs
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(inout) :: rh(:)
      real(kind=dp_t), intent(in   ) :: dx(:,:)
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: ref_ratio(:,:)
      integer        , intent(in   ) :: verbose
      integer        , intent(in   ) :: press_comp
      logical        , intent(in   ) :: nodal(:)

      type(bc_level) :: bc

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      integer        , pointer ::  mp(:,:,:,:) 

      integer :: i,n,dm,ng
      integer :: mglev_fine,lom(rh(nlevs)%dim)
      real(kind=dp_t) :: rhmax,rhsum
      real(kind=dp_t) :: local_max,local_sum
      type(      box) :: mbox
      type(      box) :: pdc
      type(   layout) :: la_crse,la_fine
      type(bndry_reg) :: brs_flx

      dm = unew(nlevs)%dim
      ng = unew(nlevs)%ng

!     rhsum = ZERO

!     Create the regular single-level divergence.
      do n = 1, nlevs
         mglev_fine = mgt(n)%nlevels
         bc = the_bc_tower%bc_tower_array(n)
         call multifab_fill_boundary(unew(n))
         do i = 1, unew(n)%nboxes
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n), i)
            rhp => dataptr(rh(n)  , i)
            mp   => dataptr(mgt(n)%mm(mglev_fine),i)
            select case (dm)
               case (2)
                 call divu_2d(unp(:,:,1,:), rhp(:,:,1,1), &
                               mp(:,:,1,1),  dx(n,:),  &
                              bc%ell_bc_level_array(i,:,:,press_comp), ng)
               case (3)
                 call divu_3d(unp(:,:,:,:), rhp(:,:,:,1), &
                               mp(:,:,:,1),  dx(n,:), &
                              bc%ell_bc_level_array(i,:,:,press_comp), ng)
            rhmax = max(rhmax,local_max)
!           rhsum = rhsum + local_sum
            end select
         end do
      end do

!     Modify the divu above at coarse-fine interfaces.
      do n = nlevs,2,-1
         la_crse = unew(n-1)%la
         la_fine = unew(n  )%la
         pdc = layout_get_pd(la_crse)
         call bndry_reg_build(brs_flx,la_fine,ref_ratio(n-1,:),pdc,nodal=nodal)
         bc = the_bc_tower%bc_tower_array(n)
         call crse_fine_divu(n,nlevs,rh(n-1),unew,brs_flx,dx,bc,ref_ratio(n-1,:),mgt(n),press_comp)
         call bndry_reg_destroy(brs_flx)
      end do

      rhmax = ZERO
      do n = 1,nlevs
        rhmax = max(rhmax,norm_inf(rh(n)))
      end do

!     if (verbose .eq. 1) print *,'... hgproject: divu max/sum ',rhmax,rhsum
      if (verbose .eq. 1) print *,'... hgproject: divu max ',rhmax

    end subroutine divu

!   ********************************************************************************************* !

    subroutine divu_2d(u,rh,mm,dx,press_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      integer        , intent(inout) :: mm(0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      integer         :: i,j,nx,ny
      real(kind=dp_t) :: vol

      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3

      rh = ZERO

      vol = dx(1) * dx(2)

      do j = 0,ny
      do i = 0,nx
         if (.not. bc_dirichlet(mm(i,j),1,0)) then
           rh(i,j) = (u(i  ,j,1) + u(i  ,j-1,1) &
                     -u(i-1,j,1) - u(i-1,j-1,1)) / dx(1) + &
                     (u(i,j  ,2) + u(i-1,j  ,2) &
                     -u(i,j-1,2) - u(i-1,j-1,2)) / dx(2)
           rh(i,j) = HALF * rh(i,j) * vol
         end if
      end do
      end do

      if (press_bc(1,1) == BC_NEU) rh( 0,:) = TWO*rh( 0,:)
      if (press_bc(1,2) == BC_NEU) rh(nx,:) = TWO*rh(nx,:)
      if (press_bc(2,1) == BC_NEU) rh(:, 0) = TWO*rh(:, 0)
      if (press_bc(2,2) == BC_NEU) rh(:,ny) = TWO*rh(:,ny)

    end subroutine divu_2d

!   ********************************************************************************************* !

    subroutine divu_3d(u,rh,mm,dx,press_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      integer        , intent(inout) :: mm(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      integer         :: i,j,k,nx,ny,nz
      real(kind=dp_t) :: vol

      nx = size(rh,dim=1) - ng
      ny = size(rh,dim=2) - ng
      nz = size(rh,dim=3) - ng

      rh = ZERO

      vol = dx(1) * dx(2) * dx(3)

      do k = 0,nz
      do j = 0,ny
      do i = 0,nx
         if (.not. bc_dirichlet(mm(i,j,k),1,0)) then
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
           rh(i,j,k) = FOURTH * rh(i,j,k) * vol
         end if
      end do
      end do
      end do

      if (press_bc(1,1) == BC_NEU) rh( 0,:,:) = TWO*rh( 0,:,:)
      if (press_bc(1,2) == BC_NEU) rh(nx,:,:) = TWO*rh(nx,:,:)
      if (press_bc(2,1) == BC_NEU) rh(:, 0,:) = TWO*rh(:, 0,:)
      if (press_bc(2,2) == BC_NEU) rh(:,ny,:) = TWO*rh(:,ny,:)
      if (press_bc(3,1) == BC_NEU) rh(:,:, 0) = TWO*rh(:,:, 0)
      if (press_bc(3,2) == BC_NEU) rh(:,:,nz) = TWO*rh(:,:,nz)

    end subroutine divu_3d

!   ********************************************************************************************* !

    subroutine crse_fine_divu(n_fine,nlevs,rh_crse,u,brs_flx,dx,bc_fine,ref_ratio,mgt,press_comp)

      integer        , intent(in   ) :: n_fine,nlevs
      type(multifab) , intent(inout) :: rh_crse
      type(multifab) , intent(inout) :: u(:)
      type(bndry_reg), intent(inout) :: brs_flx
      real(dp_t)     , intent(in   ) :: dx(:,:)
      type(bc_level) , intent(in   ) :: bc_fine
      integer        , intent(in   ) :: ref_ratio(:)
      type(mg_tower) , intent(inout) :: mgt
      integer        , intent(in   ) :: press_comp

      real(kind=dp_t), pointer :: unp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 

      type(multifab) :: temp_rhs
      type(multifab) :: vol_frac
      type(  layout) :: la_crse,la_fine
      type(     box) :: pdc
      integer :: i,dm,n_crse,ng
      integer :: mglev_fine
      logical, allocatable :: nodal(:)

      dm = u(n_fine)%dim
      n_crse = n_fine-1

      ng = u(nlevs)%ng
      allocate(nodal(dm))
      nodal = .true.

      la_crse = u(n_crse)%la
      la_fine = u(n_fine)%la

      call multifab_build(temp_rhs, la_fine, 1, 1, nodal)
      call multifab_build(vol_frac, la_crse, 1, 1, nodal)
      call setval(temp_rhs, ZERO, 1, all=.true.)

!     Zero out the flux registers which will hold the fine contributions
      call bndry_reg_setval(brs_flx, ZERO, all = .true.)

!     Compute the fine contributions at faces, edges and corners.

!     First compute a residual which only takes contributions from the
!        grid on which it is calculated.
       do i = 1, u(n_fine)%nboxes
          if ( multifab_remote(u(n_fine), i) ) cycle
          unp => dataptr(u(n_fine), i)
          rhp => dataptr( temp_rhs, i)
          select case (dm)
             case (2)
               call grid_divu_2d(unp(:,:,1,:), rhp(:,:,1,1), dx(n_fine,:), &
                                 bc_fine%ell_bc_level_array(i,:,:,press_comp), ng)
             case (3)
               call grid_divu_3d(unp(:,:,:,:), rhp(:,:,:,1), dx(n_fine,:), &
                                 bc_fine%ell_bc_level_array(i,:,:,press_comp), ng)
          end select
      end do

      pdc = layout_get_pd(la_crse)
      mglev_fine = mgt%nlevels

      do i = 1,dm
         call ml_fine_contrib(brs_flx%bmf(i,0), &
                              temp_rhs,mgt%mm(mglev_fine),ref_ratio,pdc,-i)
         call ml_fine_contrib(brs_flx%bmf(i,1), &
                              temp_rhs,mgt%mm(mglev_fine),ref_ratio,pdc,+i)
      end do

!     Compute the crse contributions at edges and corners and add to rh(n-1).
      call setval(vol_frac,ONE,all=.true.)
      do i = 1,dm
         call ml_crse_divu_contrib(rh_crse, vol_frac, brs_flx%bmf(i,0), u(n_crse), &
                                   mgt%mm(mglev_fine), dx(n_fine-1,:), &
                                   pdc,ref_ratio, -i)
         call ml_crse_divu_contrib(rh_crse, vol_frac, brs_flx%bmf(i,1), u(n_crse), &
                                   mgt%mm(mglev_fine), dx(n_fine-1,:),  &
                                   pdc,ref_ratio, +i)
      end do

    end subroutine crse_fine_divu

!   ********************************************************************************************* !

    subroutine grid_divu_2d(u,rh,dx,press_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      integer :: i,j,nx,ny
      nx = size(rh,dim=1) - 3
      ny = size(rh,dim=2) - 3

      u(-1,:,:) = ZERO
      u(nx,:,:) = ZERO
      u(:,-1,:) = ZERO
      u(:,ny,:) = ZERO

      do j = 0,ny
      do i = 0,nx
         rh(i,j) = HALF * (u(i  ,j,1) + u(i  ,j-1,1) &
                          -u(i-1,j,1) - u(i-1,j-1,1)) / dx(1) + &
                   HALF * (u(i,j  ,2) + u(i-1,j  ,2) &
                          -u(i,j-1,2) - u(i-1,j-1,2)) / dx(2)
         rh(i,j) = rh(i,j) * dx(1) * dx(2)
      end do
      end do

      if (press_bc(1,1) == BC_NEU) rh( 0,:) = TWO*rh( 0,:)
      if (press_bc(1,2) == BC_NEU) rh(nx,:) = TWO*rh(nx,:)
      if (press_bc(2,1) == BC_NEU) rh(:, 0) = TWO*rh(:, 0)
      if (press_bc(2,2) == BC_NEU) rh(:,ny) = TWO*rh(:,ny)

    end subroutine grid_divu_2d

!   ********************************************************************************************* !

    subroutine grid_divu_3d(u,rh,dx,press_bc,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::  u(-ng:,-ng:,-ng:,1:)
      real(kind=dp_t), intent(inout) :: rh(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      integer :: i,j,k,nx,ny,nz
      real(kind=dp_t) :: rh_sum

      nx = size(rh,dim=1) - ng
      ny = size(rh,dim=2) - ng
      nz = size(rh,dim=3) - ng

      u(-1,:,:,:) = ZERO
      u(nx,:,:,:) = ZERO
      u(:,-1,:,:) = ZERO
      u(:,ny,:,:) = ZERO
      u(:,:,-1,:) = ZERO
      u(:,:,nz,:) = ZERO

      rh_sum = ZERO
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
         rh(i,j,k) = FOURTH*rh(i,j,k) * (dx(1)*dx(2)*dx(3))
         rh_sum = rh_sum + rh(i,j,k)
      end do
      end do
      end do
      print *,'RHSUM ',rh_sum

      if (press_bc(1,1) == BC_NEU) rh( 0,:,:) = TWO*rh( 0,:,:)
      if (press_bc(1,2) == BC_NEU) rh(nx,:,:) = TWO*rh(nx,:,:)
      if (press_bc(2,1) == BC_NEU) rh(:, 0,:) = TWO*rh(:, 0,:)
      if (press_bc(2,2) == BC_NEU) rh(:,ny,:) = TWO*rh(:,ny,:)
      if (press_bc(3,1) == BC_NEU) rh(:,:, 0) = TWO*rh(:,:, 0)
      if (press_bc(3,2) == BC_NEU) rh(:,:,nz) = TWO*rh(:,:,nz)

    end subroutine grid_divu_3d

!   ********************************************************************************************* !

    subroutine ml_crse_divu_contrib(rh, vol, flux, u, mm, dx, crse_domain, ir, side)
     type(multifab), intent(inout) :: rh
     type(multifab), intent(inout) :: vol
     type(multifab), intent(inout) :: flux
     type(multifab), intent(in   ) :: u
     type(imultifab),intent(in   ) :: mm
     real(kind=dp_t),intent(in   ) :: dx(:)
     type(box)      ,intent(in   ) :: crse_domain
     integer        ,intent(in   ) :: ir(:)
     integer        ,intent(in   ) :: side

     type(box) :: rbox, fbox, ubox, mbox
     integer :: lo (rh%dim), hi (rh%dim)
     integer :: lou(rh%dim)
     integer :: lof(rh%dim), hif(rh%dim)
     integer :: lor(rh%dim)
     integer :: lom(rh%dim)
     integer :: lo_dom(rh%dim), hi_dom(rh%dim)
     integer :: dir
     integer :: i, j, n
     logical :: nodal(rh%dim)

     integer :: dm
     real(kind=dp_t), pointer :: rp(:,:,:,:)
     real(kind=dp_t), pointer :: vp(:,:,:,:)
     real(kind=dp_t), pointer :: fp(:,:,:,:)
     real(kind=dp_t), pointer :: up(:,:,:,:)
     integer,         pointer :: mp(:,:,:,:)

     nodal = .true.

     dm  = rh%dim
     dir = iabs(side)

     lo_dom = lwb(crse_domain)
     hi_dom = upb(crse_domain)+1

     do j = 1, u%nboxes

       ubox = get_ibox(u,j)
       lou = lwb(ubox) - u%ng
       ubox = box_nodalize(ubox,nodal)

       rbox = get_ibox(rh,j)
       lor = lwb(rbox) - rh%ng

       do i = 1, flux%nboxes

         fbox = get_ibox(flux,i)
         lof  = lwb(fbox)
         hif  = upb(fbox)

         mbox = get_ibox(mm,i)
         lom  = lwb(mbox) - mm%ng

         if ((.not. (lof(dir) == lo_dom(dir) .or. lof(dir) == hi_dom(dir))) .and. & 
             box_intersects(ubox,fbox)) then
           lo(:) = lwb(box_intersection(ubox,fbox))
           hi(:) = upb(box_intersection(ubox,fbox))

           fp => dataptr(flux,i)
           mp => dataptr(mm  ,i)

           up => dataptr(u  ,j)
           rp => dataptr(rh ,j)
           vp => dataptr(vol,j)

           select case (dm)
           case (2)
               call ml_interface_2d_divu(rp(:,:,1,1), lor, &
                                         vp(:,:,1,1), &
                                         fp(:,:,1,1), lof, hif, &
                                         up(:,:,1,:), lou, &
                                         mp(:,:,1,1), lom, &
                                         lo, hi, ir, side, dx)
!          case (3)
!              call ml_interface_3d_divu(rp(:,:,:,1), lor, &
!                                        vp(:,:,:,1), &
!                                        fp(:,:,:,1), lof, hif, &
!                                        up(:,:,:,:), lou, &
!                                        mp(:,:,:,1), lom, &
!                                        lo, hi, ir, side, dx)
           end select
         end if
       end do
     end do
    end subroutine ml_crse_divu_contrib

!   ********************************************************************************************* !

    subroutine ml_interface_2d_divu(rh, lor, vol_frac, fine_flux, lof, hif, uc, loc, &
                                    mm, lom, lo, hi, ir, side, dx)
    integer, intent(in) :: lor(:)
    integer, intent(in) :: loc(:)
    integer, intent(in) :: lom(:)
    integer, intent(in) :: lof(:), hif(:)
    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(inout) ::        rh(lor(1):,lor(2):)
    real (kind = dp_t), intent(inout) ::  vol_frac(lor(1):,lor(2):)
    real (kind = dp_t), intent(in   ) :: fine_flux(lof(1):,lof(2):)
    real (kind = dp_t), intent(in   ) ::        uc(loc(1):,loc(2):,:)
    integer           , intent(in   ) ::        mm(lom(1):,lom(2):)
    integer           , intent(in   ) :: ir(:)
    integer           , intent(in   ) :: side
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    integer :: i, j
    real (kind = dp_t) :: crse_flux,vol

    i = lo(1)
    j = lo(2)

!   NOTE: THESE STENCILS ONLY WORK FOR DX == DY.

!   NOTE: MM IS ON THE FINE GRID, NOT THE CRSE

    vol = dx(1)*dx(2)

!   Lo i side
    if (side == -1) then

      do j = lo(2),hi(2)

        if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux = FOURTH*(uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(1)))
        else if (j == lof(2) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux =      (uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux = FOURTH*(uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(1)))
        else if (j == hif(2) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux =      (uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        else
          crse_flux = (HALF*(uc(i,j,1) + uc(i,j-1,1))/dx(1) &
                      +HALF*(uc(i,j,2) - uc(i,j-1,2))/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        end if

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
           rh(i,j) = rh(i,j) - crse_flux + fine_flux(i,j)
        end if

      end do

!   Hi i side
    else if (side ==  1) then

      do j = lo(2),hi(2)

        if (j == lof(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux = FOURTH*(-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(1)))
        else if (j == lof(2) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),2,-1)) then
          crse_flux =      (-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        else if (j == hif(2) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux = FOURTH*(-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(1)))
        else if (j == hif(2) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),2,+1)) then
          crse_flux =      (-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        else
          crse_flux = (HALF*(-uc(i-1,j,1)-uc(i-1,j-1,1))/dx(1)  &
                      +HALF*( uc(i-1,j,2)-uc(i-1,j-1,2))/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(1)))
        end if

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
           rh(i,j) = rh(i,j) - crse_flux + fine_flux(i,j)
        end if

      end do

!   Lo j side
    else if (side == -2) then

      do i = lo(1),hi(1)

        if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux = FOURTH*(uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(2)))
        else if (i == lof(1) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux =      (uc(i,j,1)/dx(1) + uc(i,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux = FOURTH*(-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(2)))
        else if (i == hif(1) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux =      (-uc(i-1,j,1)/dx(1) + uc(i-1,j,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        else
          crse_flux = (HALF*(uc(i,j,1)-uc(i-1,j,1))/dx(1)  &
                      +HALF*(uc(i,j,2)+uc(i-1,j,2))/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        end if

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
           rh(i,j) = rh(i,j) - crse_flux + fine_flux(i,j)
        end if

      end do

!   Hi j side
    else if (side ==  2) then

      do i = lo(1),hi(1)

        if (i == lof(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux = FOURTH*(uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(2)))
        else if (i == lof(1) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),1,-1)) then
          crse_flux =      (uc(i,j-1,1)/dx(1) - uc(i,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        else if (i == hif(1) .and. .not. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux = FOURTH*(-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - FOURTH*(ONE-ONE/dble(ir(2)))
        else if (i == hif(1) .and. bc_neumann(mm(ir(1)*i,ir(2)*j),1,+1)) then
          crse_flux =      (-uc(i-1,j-1,1)/dx(1) - uc(i-1,j-1,2)/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        else
          crse_flux = (HALF*( uc(i,j-1,1)-uc(i-1,j-1,1))/dx(1) &
                      +HALF*(-uc(i,j-1,2)-uc(i-1,j-1,2))/dx(2)) * vol
          vol_frac(i,j) = vol_frac(i,j) - HALF*(ONE-ONE/dble(ir(2)))
        end if

        if (bc_dirichlet(mm(ir(1)*i,ir(2)*j),1,0)) then
           rh(i,j) = rh(i,j) - crse_flux + fine_flux(i,j)
        end if

      end do

    end if

!   do j = lor(2),lor(2)+size(rh,dim=2)-1
!   do i = lor(1),lor(1)+size(rh,dim=1)-1
!      rh(i,j) = rh(i,j) * vol_frac(i,j)
!   end do
!   end do

  end subroutine ml_interface_2d_divu

end module hgproject_module
