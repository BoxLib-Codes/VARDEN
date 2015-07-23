module viscous_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module
  use bndry_reg_module

  implicit none

  private

  public :: visc_solve, diff_scalar_solve

contains 

  subroutine visc_solve(mla,unew,lapu,rho,mac_rhs,dx,mu,the_bc_tower)

    use mac_multigrid_module    , only : mac_multigrid
    use probin_module           , only: stencil_order,verbose

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: lapu(:)
    type(multifab ), intent(in   ) :: rho(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),mu
    type(bc_tower ), intent(in   ) :: the_bc_tower
    type(multifab ), intent(in   ) :: mac_rhs(:)

    ! Local  
    type(multifab)  :: rh(mla%nlevel),phi(mla%nlevel)
    type(multifab)  :: alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    integer         :: n,d,nlevs,dm,bc_comp,ng_cell
    real(kind=dp_t) :: nrm1, nrm2, nrm3
    real(kind=dp_t) :: rel_solver_eps
    real(kind=dp_t) :: abs_solver_eps

    type(bl_prof_timer), save :: bpt

    call build(bpt,"visc_solve")

    nlevs = mla%nlevel
    dm    = mla%dim
    ng_cell = nghost(unew(1))

    do n = 1,nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(  phi(n), mla%la(n),  1, 1)
       call multifab_build(alpha(n), mla%la(n),  1, 0)
       do d = 1,dm
          call multifab_build_edge(beta(n,d), mla%la(n), 1, 0, d)
       end do

       call multifab_copy_c(alpha(n),1,rho(n),1,1)
       do d = 1,dm
          call setval( beta(n,d), mu,all=.true.)
       end do
    end do

    if (verbose .ge. 1) then
       if (parallel_IOProcessor()) then
          print *,' '
          print *,'... begin viscous solves  ... '
       end if
       do n = 1,nlevs
          nrm1 = norm_inf(unew(n),1,1)
          nrm2 = norm_inf(unew(n),2,1)
          if ( parallel_IOProcessor() ) then
             print *,'BEFORE: MAX OF U AT LEVEL ',n,nrm1
             print *,'BEFORE: MAX OF V AT LEVEL ',n,nrm2
          end if
          if (dm .eq. 3) then
             nrm3 = norm_inf(unew(n),3,1)
             if ( parallel_IOProcessor() ) then
                print *,'BEFORE: MAX OF W AT LEVEL ',n,nrm3
             end if
          end if
       end do
    endif

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    rel_solver_eps =  1.d-12
    abs_solver_eps = -1.d0

    ! Note that "mu" here is actually (1/2 dt mu) if Crank-Nicolson -- 
    !                              or (    dt mu) if backward Euler
    !      this is set in velocity_advance
    do d = 1,dm
       do n = 1,nlevs
          call mkrhs(rh(n),unew(n),lapu(n),rho(n),phi(n),mac_rhs(n),mu,dx(n,:),d)
       end do
       bc_comp = d
       call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx, &
                          the_bc_tower,bc_comp,stencil_order, &
                          rel_solver_eps,abs_solver_eps)
       do n = 1,nlevs
          call multifab_copy_c(unew(n),d,phi(n),1,1)
       end do
    end do

    call ml_restrict_and_fill(nlevs, unew, mla%mba%rr, the_bc_tower%bc_tower_array)

    if ( verbose .ge. 1 ) then
       do n = 1,nlevs
          nrm1 = norm_inf(unew(n),1,1)
          nrm2 = norm_inf(unew(n),2,1)
          if ( parallel_IOProcessor() ) then
             print *,' AFTER: MAX OF U AT LEVEL ',n,nrm1
             print *,' AFTER: MAX OF V AT LEVEL ',n,nrm2
          end if
          if (dm .eq. 3) then
             nrm3 = norm_inf(unew(n),3,1)
             if ( parallel_IOProcessor() ) then
                print *,' AFTER: MAX OF W AT LEVEL ',n,nrm3
             end if
          end if
       end do

       if ( parallel_IOProcessor() ) then
          print *,'...   end viscous solves  ... '
          print *,' '
       end if
    endif

    do n = 1, nlevs
       call multifab_destroy(rh(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(alpha(n))
       do d = 1,dm
          call multifab_destroy(beta(n,d))
       end do
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    call destroy(bpt)

  contains

    subroutine mkrhs(rh,unew,lapu,rho,phi,mac_rhs,mu,dx,comp)

      type(multifab) , intent(in   ) :: unew,lapu,rho
      type(multifab) , intent(inout) :: rh,phi
      type(multifab ), intent(in   ) :: mac_rhs
      integer        , intent(in   ) :: comp
      real(dp_t)     , intent(in   ) :: mu, dx(:)

      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: rhp(:,:,:,:)
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
      real(kind=dp_t), pointer ::  pp(:,:,:,:)
      real(kind=dp_t), pointer ::  lp(:,:,:,:)
      real(kind=dp_t), pointer ::  mp(:,:,:,:)
      integer :: i,dm,ng_u,ng_r,ng_m

      type(bl_prof_timer), save :: bpt

      call build(bpt,"mkrhs")

      dm     = get_dim(rh)
      ng_u = nghost(unew)
      ng_r = nghost(rho)
      ng_m = nghost(mac_rhs)

      do i = 1, nfabs(unew)
         rhp => dataptr(rh  , i)
         unp => dataptr(unew, i)
         rp => dataptr(rho , i)
         pp => dataptr(phi , i)
         lp => dataptr(lapu, i)
         mp => dataptr(mac_rhs, i)
         select case (dm)
         case (2)
            call mkrhs_2d(rhp(:,:,1,1), unp(:,:,1,comp), lp(:,:,1,comp), rp(:,:,1,1), &
                           pp(:,:,1,1),  mp(:,:,1,1), mu, dx, ng_u, ng_r, ng_m, comp)
         case (3)
            call mkrhs_3d(rhp(:,:,:,1), unp(:,:,:,comp), lp(:,:,:,comp), rp(:,:,:,1), &
                           pp(:,:,:,1),  mp(:,:,:,1), mu, dx, ng_u, ng_r, ng_m, comp)
         end select
      end do

      call destroy(bpt)

    end subroutine mkrhs

    subroutine mkrhs_2d(rh,unew,lapu,rho,phi,mac_rhs,mu,dx,ng_u,ng_r,ng_m,comp)

      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng_u,ng_r,ng_m,comp
      real(kind=dp_t), intent(inout) ::      rh(      :,      :)
      real(kind=dp_t), intent(in   ) ::    unew(1-ng_u:,1-ng_u  :)
      real(kind=dp_t), intent(in   ) ::    lapu(     1:,     1:)
      real(kind=dp_t), intent(in   ) ::     rho(1-ng_r:,1-ng_r:)
      real(kind=dp_t), intent(inout) ::     phi(     0:,     0:)
      real(kind=dp_t), intent(in   ) :: mac_rhs(1-ng_m:,1-ng_m:)
      real(dp_t)     , intent(in   ) :: mu,dx(:)

      integer    :: nx,ny,i,j
      real(dp_t) :: visc_mu_dt

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)

       rh(1:nx  ,1:ny  ) = unew(1:nx  ,1:ny  ) * rho(1:nx,1:ny)
      phi(0:nx+1,0:ny+1) = unew(0:nx+1,0:ny+1)

      ! Note that "mu" here is actually (1/2 dt mu) if Crank-Nicolson -- 
      !                              or (    dt mu) if backward Euler
      !      this is set in velocity_advance
      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny) = rh(1:nx,1:ny) + mu*lapu(1:nx,1:ny)
         visc_mu_dt = TWO *mu
      else if (diffusion_type .eq. 2) then
         visc_mu_dt = mu
      end if

      ! Add (1/3) * mu * dt * grad (divu) to the RHS of the viscous solve
      if (comp.eq.1) then
         do j = 1,ny
         do i = 1,nx
             rh(i,j) = rh(i,j) + THIRD*visc_mu_dt*(mac_rhs(i+1,j) - mac_rhs(i-1,j)) / dx(1)
         end do
         end do
      else
         do j = 1,ny
         do i = 1,nx
             rh(i,j) = rh(i,j) + THIRD*visc_mu_dt*(mac_rhs(i,j+1) - mac_rhs(i,j-1)) / dx(2)
         end do
         end do
      end if

    end subroutine mkrhs_2d

    subroutine mkrhs_3d(rh,unew,lapu,rho,phi,mac_rhs,mu,dx,ng_u,ng_r,ng_m,comp)

      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng_u,ng_r,ng_m,comp
      real(kind=dp_t), intent(inout)    ::   rh(      :,      :,      :)
      real(kind=dp_t), intent(in   ) ::    unew(1-ng_u:,1-ng_u:,1-ng_u:)
      real(kind=dp_t), intent(inout) ::    lapu(     1:,     1:,     1:)
      real(kind=dp_t), intent(in   ) ::     rho(1-ng_r:,1-ng_r:,1-ng_r:)
      real(kind=dp_t), intent(inout) ::     phi(     0:,     0:,     0:)
      real(kind=dp_t), intent(in   ) :: mac_rhs(1-ng_m:,1-ng_m:,1-ng_m:)
      real(dp_t)     , intent(in   ) :: mu,dx(:)

      integer    :: nx,ny,nz,i,j,k
      real(dp_t) :: visc_mu_dt

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)
      nz = size(rh,dim=3)

      phi(0:nx+1,0:ny+1,0:nz+1) = unew(0:nx+1,0:ny+1,0:nz+1)
      rh(1:nx  ,1:ny  ,1:nz  ) = unew(1:nx  ,1:ny  ,1:nz  ) * &
                                  rho(1:nx  ,1:ny  ,1:nz  )

      ! Note that "mu" here is actually (1/2 dt mu) if Crank-Nicolson -- 
      !                              or (    dt mu) if backward Euler
      !      this is set in velocity_advance
      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny,1:nz) = rh(1:nx,1:ny,1:nz) + mu*lapu(1:nx,1:ny,1:nz) 
         visc_mu_dt = TWO *mu
      else if (diffusion_type .eq. 2) then
         visc_mu_dt = mu
      end if

      ! Add (1/3) * mu * dt * grad (divu) to the RHS of the viscous solve
      if (comp.eq.1) then
         do k = 1,nz
         do j = 1,ny
         do i = 1,nx
             rh(i,j,k) = rh(i,j,k) + THIRD*visc_mu_dt*(mac_rhs(i+1,j,k) - mac_rhs(i-1,j,k)) / dx(1)
         end do
         end do
         end do
      else if (comp.eq.2) then
         do k = 1,nz
         do j = 1,ny
         do i = 1,nx
             rh(i,j,k) = rh(i,j,k) + THIRD*visc_mu_dt*(mac_rhs(i,j+1,k) - mac_rhs(i,j-1,k)) / dx(2)
         end do
         end do
         end do
      else 
         do k = 1,nz
         do j = 1,ny
         do i = 1,nx
             rh(i,j,k) = rh(i,j,k) + THIRD*visc_mu_dt*(mac_rhs(i,j,k+1) - mac_rhs(i,j,k-1)) / dx(3)
         end do
         end do
         end do
      end if

    end subroutine mkrhs_3d

  end subroutine visc_solve

  subroutine diff_scalar_solve(mla,snew,laps,dx,mu,the_bc_tower,icomp,bc_comp)

    use mac_multigrid_module    , only : mac_multigrid
    use ml_cc_restriction_module, only : ml_cc_restriction_c
    use probin_module           , only: stencil_order,verbose

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: snew(:)
    type(multifab ), intent(in   ) :: laps(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: mu
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: icomp,bc_comp

    ! Local  
    type(multifab)  :: rh(mla%nlevel),phi(mla%nlevel)
    type(multifab)  :: alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    integer         :: n,d,nlevs,dm,ng_cell
    real(kind=dp_t) :: rel_solver_eps
    real(kind=dp_t) :: abs_solver_eps
    real(kind=dp_t) :: nrm1

    type(bl_prof_timer), save :: bpt

    call build(bpt,"diff_scalar_solve")

    nlevs = mla%nlevel
    dm    = mla%dim
    ng_cell = nghost(snew(1))

    do n = 1,nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(  phi(n), mla%la(n),  1, 1)
       call multifab_build(alpha(n), mla%la(n),  1, 1)
       do d = 1,dm
          call multifab_build_edge(beta(n,d), mla%la(n), 1, 1, d)
       end do
       call setval(alpha(n),ONE,all=.true.)
       do d = 1,dm
          call setval( beta(n,d), mu,all=.true.)
       end do
    end do

    if ( verbose .ge. 1 ) then
       if (parallel_IOProcessor()) then
          print *,' '
          print *,'... begin diffusive solve  ... '
       end if
       do n = 1,nlevs
          nrm1 = norm_inf(snew(n),icomp,1)
          if ( parallel_IOProcessor() ) print *,'BEFORE: MAX OF S AT LEVEL ',n,nrm1
       end do
    end if

    do n = 1,nlevs
       call mkrhs(rh(n),snew(n),laps(n),phi(n),mu,icomp)
    end do

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    rel_solver_eps =  1.d-12
    abs_solver_eps = -1.d0

    call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx, &
                       the_bc_tower,bc_comp,stencil_order, &
                       rel_solver_eps,abs_solver_eps)

    do n = 1,nlevs
       call multifab_copy_c(snew(n),icomp,phi(n),1,1)
    end do

    do n = 1, nlevs
       call multifab_fill_boundary_c(snew(n),icomp,1)
       call multifab_physbc(snew(n),icomp,bc_comp,1,the_bc_tower%bc_tower_array(n))
    enddo

    do n = nlevs, 2, -1
       call ml_cc_restriction_c(snew(n-1),icomp,snew(n),icomp,mla%mba%rr(n-1,:),1)
    end do

    do n = 2, nlevs
       call multifab_fill_ghost_cells(snew(n),snew(n-1), &
                                      ng_cell,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n  ), &
                                      icomp,bc_comp,1)
    end do

    if ( verbose .ge. 1 ) then
       do n = 1,nlevs
          nrm1 = norm_inf(snew(n),icomp,1)
          if ( parallel_IOProcessor() ) print *,'AFTER: MAX OF S AT LEVEL ',n,nrm1
       end do
       if ( parallel_IOProcessor() ) then
          print *,' '
          print *,'...   end diffusive solve  ... '
       end if
    endif

    do n = 1, nlevs
       call multifab_destroy(rh(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(alpha(n))
       do d = 1,dm
          call multifab_destroy(beta(n,d))
       end do
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    call destroy(bpt)

  contains

    subroutine mkrhs(rh,snew,laps,phi,mu,comp)

      type(multifab) , intent(in   ) :: snew,laps
      type(multifab) , intent(inout) :: rh,phi
      real(dp_t)     , intent(in   ) :: mu
      integer        , intent(in   ) :: comp

      real(kind=dp_t), pointer :: sp(:,:,:,:)
      real(kind=dp_t), pointer :: rp(:,:,:,:)
      real(kind=dp_t), pointer :: pp(:,:,:,:)
      real(kind=dp_t), pointer :: lp(:,:,:,:)
      integer :: i,dm,ng

      type(bl_prof_timer), save :: bpt

      call build(bpt,"diff_scalar_solve/mkrhs")

      dm   = get_dim(rh)
      ng   = nghost(snew)

      do i = 1, nfabs(snew)
         rp => dataptr(rh  , i)
         pp => dataptr(phi , i)
         sp => dataptr(snew, i)
         lp => dataptr(laps, i)
         select case (dm)
         case (2)
            call mkrhs_2d(rp(:,:,1,1), sp(:,:,1,comp), lp(:,:,1,comp), pp(:,:,1,1), mu, ng)
         case (3)
            call mkrhs_3d(rp(:,:,:,1), sp(:,:,:,comp), lp(:,:,:,comp), pp(:,:,:,1), mu, ng)
         end select
      end do

      call destroy(bpt)

    end subroutine mkrhs

    subroutine mkrhs_2d(rh,snew,laps,phi,mu,ng)
      
      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::   rh(    :,    :)
      real(kind=dp_t), intent(in   ) :: snew(1-ng:,1-ng:)
      real(kind=dp_t), intent(in   ) :: laps(1   :,   1:)
      real(kind=dp_t), intent(inout) ::  phi(0   :,   0:)
      real(dp_t)     , intent(in   ) :: mu

      integer :: nx,ny

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)

      rh(1:nx,1:ny) = snew(1:nx,1:ny)
      phi(0:nx+1,0:ny+1) = snew(0:nx+1,0:ny+1)

      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny) = rh(1:nx,1:ny) + mu*laps(1:nx,1:ny)
      end if

    end subroutine mkrhs_2d

    subroutine mkrhs_3d(rh,snew,laps,phi,mu,ng)

      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::   rh(    :,    :,    :)
      real(kind=dp_t), intent(in   ) :: snew(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(in   ) :: laps(1   :,   1:,   1:)
      real(kind=dp_t), intent(inout) ::  phi(0   :,   0:,   0:)
      real(dp_t)     , intent(in   ) :: mu

      integer :: nx,ny,nz

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)
      nz = size(rh,dim=3)

      rh(1:nx,1:ny,1:nz) = snew(1:nx,1:ny,1:nz)
      phi(0:nx+1,0:ny+1,0:nz+1) = snew(0:nx+1,0:ny+1,0:nz+1)

      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny,1:nz) = rh(1:nx,1:ny,1:nz) + mu*laps(1:nx,1:ny,1:nz)
      end if

    end subroutine mkrhs_3d

  end subroutine diff_scalar_solve

end module viscous_module
