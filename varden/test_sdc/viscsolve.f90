module viscous_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: visc_solve, diff_scalar_solve

contains 

  subroutine visc_solve(mla,unew,lapu,rho,dx,t,mu,the_bc_tower)

    use bl_constants_module
    use bndry_reg_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use macproject_module,     only: mac_multigrid
    use ml_restriction_module, only: ml_cc_restriction
    use probin_module, only: stencil_order,verbose

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: lapu(:)
    type(multifab ), intent(in   ) :: rho(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),t,mu
    type(bc_tower ), intent(in   ) :: the_bc_tower
 
    ! Local  
    type(multifab)  :: rh(mla%nlevel),phi(mla%nlevel)
    type(multifab)  ::alpha(mla%nlevel),beta(mla%nlevel)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    integer                     :: n,nlevs,d,dm
    integer                     :: bc_comp,ng_cell
    real(kind=dp_t)             :: nrm1, nrm2

    nlevs = mla%nlevel
    dm    = mla%dim
    ng_cell = unew(1)%ng

    do n = 1,nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(  phi(n), mla%la(n),  1, 1)
       call multifab_build(alpha(n), mla%la(n),  1, 1)
       call multifab_build( beta(n), mla%la(n), dm, 1)

       call multifab_copy_c(alpha(n),1,rho(n),1,1)
       call setval(beta(n),mu,all=.true.)
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
       end do
    endif

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do d = 1,dm
       do n = 1,nlevs
          call mkrhs(rh(n),unew(n),lapu(n),rho(n),phi(n),mu,d)
       end do
       bc_comp = d
       call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx, &
                          the_bc_tower,bc_comp,stencil_order,mla%mba%rr)
       do n = 1,nlevs
          call multifab_copy_c(unew(n),d,phi(n),1,1)
       end do
    end do

    do n = 1, nlevs
       call multifab_fill_boundary(unew(n))
       call multifab_physbc(unew(n),1,1,dm,the_bc_tower%bc_tower_array(n),dx(n,:),t)
    enddo

    do n = nlevs, 2, -1
       call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(unew(n),unew(n-1), &
                                      ng_cell,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n  ), &
                                      1,1,dm,dx(n-1:n,:),t)
    end do

    if ( verbose .ge. 1 ) then
       do n = 1,nlevs
          nrm1 = norm_inf(unew(n),1,1)
          nrm2 = norm_inf(unew(n),2,1)
          if ( parallel_IOProcessor() ) then
             print *,' AFTER: MAX OF U AT LEVEL ',n,nrm1
             print *,' AFTER: MAX OF V AT LEVEL ',n,nrm2
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
       call multifab_destroy(beta(n))
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

  contains

    subroutine mkrhs(rh,unew,lapu,rho,phi,mu,comp)

      type(multifab) , intent(in   ) :: unew,lapu,rho
      type(multifab) , intent(inout) :: rh,phi
      integer        , intent(in   ) :: comp
      real(dp_t)     , intent(in   ) :: mu

      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: rhp(:,:,:,:)
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
      real(kind=dp_t), pointer ::  pp(:,:,:,:)
      real(kind=dp_t), pointer ::  lp(:,:,:,:)
      integer :: i,dm,ng_u,ng_rho

      dm     = rh%dim
      ng_u   = unew%ng
      ng_rho = rho%ng

      do i = 1, unew%nboxes
         if ( multifab_remote(unew, i) ) cycle
         rhp => dataptr(rh  , i)
         unp => dataptr(unew, i)
         rp => dataptr(rho , i)
         pp => dataptr(phi , i)
         lp => dataptr(lapu, i)
         select case (dm)
         case (2)
            call mkrhs_2d(rhp(:,:,1,1), unp(:,:,1,comp), lp(:,:,1,comp), rp(:,:,1,1), &
                          pp(:,:,1,1), mu, ng_u, ng_rho)
         case (3)
            call mkrhs_3d(rhp(:,:,:,1), unp(:,:,:,comp), lp(:,:,:,comp), rp(:,:,:,1), &
                          pp(:,:,:,1), mu, ng_u, ng_rho)
         end select
      end do

    end subroutine mkrhs

    subroutine mkrhs_2d(rh,unew,lapu,rho,phi,mu,ng_u,ng_rho)

      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(inout) ::   rh(        :,        :)
      real(kind=dp_t), intent(in   ) :: unew(1-ng_u  :,1-ng_u  :)
      real(kind=dp_t), intent(in   ) :: lapu(       1:,       1:)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(inout) ::  phi(       0:,       0:)
      real(dp_t)     , intent(in   ) :: mu

      integer :: nx,ny

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)

      rh(1:nx  ,1:ny  ) = unew(1:nx  ,1:ny  ) * rho(1:nx,1:ny)
      phi(0:nx+1,0:ny+1) = unew(0:nx+1,0:ny+1)

      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny) = rh(1:nx,1:ny) + mu*lapu(1:nx,1:ny)
      end if

    end subroutine mkrhs_2d

    subroutine mkrhs_3d(rh,unew,lapu,rho,phi,mu,ng_u,ng_rho)

      use probin_module, only: diffusion_type

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(inout) ::   rh(        :,        :,        :)
      real(kind=dp_t), intent(in   ) :: unew(1-ng_u  :,1-ng_u  :,1-ng_u  :)
      real(kind=dp_t), intent(inout) :: lapu(       1:,       1:,       1:)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng_rho:,1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(inout) ::  phi(       0:,       0:,       0:)
      real(dp_t)     , intent(in   ) :: mu

      integer :: nx,ny,nz

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)
      nz = size(rh,dim=3)

      phi(0:nx+1,0:ny+1,0:nz+1) = unew(0:nx+1,0:ny+1,0:nz+1)
      rh(1:nx  ,1:ny  ,1:nz  ) = unew(1:nx  ,1:ny  ,1:nz  ) * &
           rho(1:nx  ,1:ny  ,1:nz  )

      if (diffusion_type .eq. 1) then
         rh(1:nx,1:ny,1:nz) = rh(1:nx,1:ny,1:nz) + mu*lapu(1:nx,1:ny,1:nz)
      end if

    end subroutine mkrhs_3d

  end subroutine visc_solve

  subroutine diff_scalar_solve(mla,snew,laps,dx,t,mu,the_bc_tower,icomp,bc_comp)

    use bndry_reg_module
    use bl_constants_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use macproject_module,     only: mac_multigrid
    use ml_restriction_module, only: ml_cc_restriction_c
    use probin_module, only: stencil_order, verbose, mass_fractions, nscal

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: snew(:)
    type(multifab ), intent(in   ) :: laps(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    real(dp_t)     , intent(in   ) :: t
    real(dp_t)     , intent(in   ) :: mu
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: icomp,bc_comp

    ! Local  
    type(multifab)  :: rh(mla%nlevel),phi(mla%nlevel)
    type(multifab)  :: alpha(mla%nlevel),beta(mla%nlevel)
    type(bndry_reg) :: fine_flx(mla%nlevel)
    integer         :: n,nlevs,dm
    integer         :: ng_cell
    real(kind=dp_t) :: nrm1

    nlevs = mla%nlevel
    dm    = mla%dim
    ng_cell = snew(1)%ng
    
    if (mass_fractions) then
       do n = 1,nlevs
          call multifab_build(   rh(n), mla%la(n),  1, 0)
          call multifab_build(  phi(n), mla%la(n),  1, 1)
          call multifab_build(alpha(n), mla%la(n),  1, 1)
          call multifab_build( beta(n), mla%la(n), dm, 1)

          ! not filling ghost cells here.  change?
          call multifab_copy_c(alpha(n),1,snew(n),1,1)
          call setval(beta(n),mu,all=.true.)
          call multifab_mult_mult_c(beta(n),1,snew(n),1,1)
       end do
      
    else
       do n = 1,nlevs
          call multifab_build(   rh(n), mla%la(n),  1, 0)
          call multifab_build(  phi(n), mla%la(n),  1, 1)
          call multifab_build(alpha(n), mla%la(n),  1, 1)
          call multifab_build( beta(n), mla%la(n), dm, 1)
        
          call setval(alpha(n),ONE,all=.true.)
          call setval( beta(n), mu,all=.true.)
       end do
    end if

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

    call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx, &
                       the_bc_tower,bc_comp,stencil_order,mla%mba%rr)

    do n = 1,nlevs
       call multifab_copy_c(snew(n),icomp,phi(n),1,1)
    end do

    do n = 1, nlevs
       call multifab_fill_boundary_c(snew(n),icomp,1)
       call multifab_physbc(snew(n),icomp,bc_comp,1,the_bc_tower%bc_tower_array(n),dx(n,:),t)
    enddo

    do n = nlevs, 2, -1
       call ml_cc_restriction_c(snew(n-1),icomp,snew(n),icomp,mla%mba%rr(n-1,:),1)
       call multifab_fill_ghost_cells(snew(n),snew(n-1), &
                                      ng_cell,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1), &
                                      the_bc_tower%bc_tower_array(n  ), &
                                      icomp,bc_comp,1,dx(n-1:n,:),t)
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
       call multifab_destroy(beta(n))
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

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

      dm   = rh%dim
      ng   = snew%ng

      do i = 1, snew%nboxes
         if ( multifab_remote(snew, i) ) cycle
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
