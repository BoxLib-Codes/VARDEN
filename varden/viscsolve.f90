module viscous_module

  use bl_types
  use bc_module
  use multifab_module
  use boxarray_module
  use stencil_module
  use macproject_module

  implicit none

  real (kind = dp_t), private, parameter :: TWO   = 2.0_dp_t

contains 

subroutine visc_solve(unew,rho,dx,mu,visc_bc)

  type(multifab), intent(inout) :: unew
  type(multifab), intent(in   ) :: rho
  real(dp_t)    , intent(in   ) :: dx(:),mu
  integer       , intent(in   ) :: visc_bc(:,:)

! Local  
  type(multifab) :: rh,phi,alpha,beta
  type(layout) :: la
  integer  :: n,dm,stencil_order

  la = unew%la
  dm = rho%dim

  call multifab_build(   rh, la,  1, 0)
  call multifab_build(  phi, la,  1, 1)
  call multifab_build(alpha, la,  1, 1)
  call multifab_build( beta, la, dm, 1)

! call multifab_copy_c(alpha,1,rho,1,1)
! call setval(beta,mu,all=.true.)
! HACK HACK HACK 
  call setval(alpha,0.0_dp_t,all=.true.)
  call setval(beta,1.0_dp_t,all=.true.)

  stencil_order = 1

  do n = 1,dm
     call mkrhs(rh,unew,rho,phi,n)
     call mac_multigrid(la,rh,phi,alpha,beta,dx,visc_bc,stencil_order)
     call multifab_copy_c(unew,n,phi,1,1)
  end do

  call multifab_destroy(rh)
  call multifab_destroy(phi)
  call multifab_destroy(alpha)
  call multifab_destroy(beta)

  contains

    subroutine mkrhs(rh,unew,rho,phi,comp)

      type(multifab) , intent(in   ) :: unew,rho
      type(multifab) , intent(inout) :: rh,phi
      integer        , intent(in   ) :: comp

      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: rhp(:,:,:,:)
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
      real(kind=dp_t), pointer ::  pp(:,:,:,:)
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
         select case (dm)
            case (2)
              call mkrhs_2d(rhp(:,:,1,1), unp(:,:,1,comp), rp(:,:,1,1), &
                            pp(:,:,1,1), ng_u, ng_rho)
            case (3)
              call mkrhs_3d(rhp(:,:,:,1), unp(:,:,:,comp), rp(:,:,:,1), &
                            pp(:,:,:,1), ng_u, ng_rho)
         end select
      end do

    end subroutine mkrhs

    subroutine mkrhs_2d(rh,unew,rho,phi,ng_u,ng_rho)

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(inout) ::   rh(        :,        :)
      real(kind=dp_t), intent(in   ) :: unew(1-ng_u  :,1-ng_u  :)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(inout) ::  phi(       0:,       0:)

      integer :: i,j,nx,ny

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)

       rh(1:nx  ,1:ny  ) = unew(1:nx  ,1:ny  ) * rho(1:nx,1:ny)
      phi(0:nx+1,0:ny+1) = unew(0:nx+1,0:ny+1)

    end subroutine mkrhs_2d

    subroutine mkrhs_3d(rh,unew,rho,phi,ng_u,ng_rho)

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(inout) ::   rh(        :,        :,        :)
      real(kind=dp_t), intent(in   ) :: unew(1-ng_u  :,1-ng_u  :,1-ng_u  :)
      real(kind=dp_t), intent(in   ) ::  rho(1-ng_rho:,1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(inout) ::  phi(       0:,       0:,       0:)

      integer :: i,j,k,nx,ny,nz

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)
      nz = size(rh,dim=3)

      phi(0:nx+1,0:ny+1,0:nz+1) = unew(0:nx+1,0:ny+1,0:nz+1)
       rh(1:nx  ,1:ny  ,1:nz  ) = unew(1:nx  ,1:ny  ,1:nz  ) * &
                                   rho(1:nx  ,1:ny  ,1:nz  )

    end subroutine mkrhs_3d

end subroutine visc_solve

subroutine diff_scalar_solve(snew,dx,mu,visc_bc,comp)

  type(multifab), intent(inout) :: snew
  real(dp_t)    , intent(in   ) :: dx(:),mu
  integer       , intent(in   ) :: visc_bc(:,:)
  integer       , intent(in   ) :: comp

! Local  
  type(multifab) :: rh,phi,alpha,beta
  type(layout) :: la
  integer  :: n,dm,stencil_order

  la = snew%la
  dm = snew%dim

  call multifab_build(   rh, la,  1, 0)
  call multifab_build(  phi, la,  1, 1)
  call multifab_build(alpha, la,  1, 1)
  call multifab_build( beta, la, dm, 1)
 
  call setval(alpha,1.0_dp_t,all=.true.)
  call setval(beta,mu,all=.true.)
  
  call mkrhs(rh,snew,phi,comp)

  stencil_order = 2

  call mac_multigrid(la,rh,phi,alpha,beta,dx,visc_bc,stencil_order)

  call multifab_plus_plus_c(snew,comp,phi,1,1)

  call multifab_destroy(rh)
  call multifab_destroy(phi)
  call multifab_destroy(alpha)
  call multifab_destroy(beta)

  contains

    subroutine mkrhs(rh,snew,phi,comp)

      type(multifab) , intent(in   ) :: snew
      type(multifab) , intent(inout) :: rh,phi
      integer        , intent(in   ) :: comp

      real(kind=dp_t), pointer :: sp(:,:,:,:)
      real(kind=dp_t), pointer :: rp(:,:,:,:)
      real(kind=dp_t), pointer :: pp(:,:,:,:)
      integer :: i,dm,ng

      dm   = rh%dim
      ng   = snew%ng

      do i = 1, snew%nboxes
         if ( multifab_remote(snew, i) ) cycle
         rp => dataptr(rh  , i)
         pp => dataptr(phi , i)
         sp => dataptr(snew, i)
         select case (dm)
            case (2)
              call mkrhs_2d(rp(:,:,1,1), sp(:,:,1,comp), pp(:,:,1,1), ng)
            case (3)
              call mkrhs_3d(rp(:,:,:,1), sp(:,:,:,comp), pp(:,:,:,1), ng)
         end select
      end do

    end subroutine mkrhs

    subroutine mkrhs_2d(rh,snew,phi,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::   rh(    :,    :)
      real(kind=dp_t), intent(in   ) :: snew(1-ng:,1-ng:)
      real(kind=dp_t), intent(inout) ::  phi(0   :,   0:)

      integer :: i,j,nx,ny

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)

      rh(1:nx,1:ny) = snew(1:nx,1:ny)
      phi(0:nx+1,0:ny+1) = snew(0:nx+1,0:ny+1)

    end subroutine mkrhs_2d

    subroutine mkrhs_3d(rh,snew,phi,ng)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::   rh(    :,    :,    :)
      real(kind=dp_t), intent(in   ) :: snew(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(inout) ::  phi(0   :,   0:,   0:)

      integer :: i,j,k,nx,ny,nz

      nx = size(rh,dim=1)
      ny = size(rh,dim=2)
      nz = size(rh,dim=3)

      rh(1:nx,1:ny,1:nz) = snew(1:nx,1:ny,1:nz)
      phi(0:nx+1,0:ny+1,0:nz+1) = snew(0:nx+1,0:ny+1,0:nz+1)

    end subroutine mkrhs_3d

end subroutine diff_scalar_solve

end module viscous_module
