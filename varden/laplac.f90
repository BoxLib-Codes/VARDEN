module laplac_module

   use bl_types
   use bc_module

   implicit none

contains
      subroutine laplac_2d(u,diff,dx,ng,bc)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(in   ) :: u(1-ng:,1-ng:)
      real(kind=dp_t), intent(inout) :: diff(0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      real(kind=dp_t) ::  ux_lft, ux_lft_wall
      real(kind=dp_t) ::  ux_rgt, ux_rgt_wall
      real(kind=dp_t) ::  uy_bot, uy_bot_wall
      real(kind=dp_t) ::  uy_top, uy_top_wall
      integer         ::  i,j,nx,ny

! NOTE: This assumes Dirichlet boundary conditions on walls and inflow.
!       The value is imposed at the face, not the ghost cell.
      nx = size(diff,dim=1)
      ny = size(diff,dim=2)

      do j = 0,ny-1
      do i = 0,nx-1
            ux_lft = (u(i,j) - u(i-1,j))
            ux_lft_wall = (-16.0_dp_t * u(0,j) + 20.0_dp_t * u(1,j) &
                            -5.0_dp_t * u(2,j) + u(3,j) ) * 0.2_dp_t
            ux_lft = merge(ux_lft_wall, ux_lft, i .eq. 1 .and. &
                          (bc(1,1) .eq. BC_DIR)) 
            ux_lft = ux_lft / dx(1)
 
            ux_rgt = (u(i+1,j) - u(i,j))
            ux_rgt_wall = -(-16.0_dp_t * u(nx+1,j) + 20.0_dp_t * u(nx,j) &
                             -5.0_dp_t * u(nx-1,j) + u(nx-2,j) ) * 0.2_dp_t
            ux_rgt = merge(ux_rgt_wall, ux_rgt, i .eq. nx .and. &
                            (bc(1,2) .eq. BC_DIR) )
            ux_rgt = ux_rgt / dx(1)
 
            uy_bot = (u(i,j) - u(i,j-1))
            uy_bot_wall = (-16.0_dp_t * u(i,0) + 20.0_dp_t * u(i,1) &
                            -5.0_dp_t * u(i,2) + u(i,3) ) * 0.2_dp_t
            uy_bot = merge(uy_bot_wall, uy_bot, j .eq. 1 .and. &
                            (bc(2,1) .eq. BC_DIR) )
            uy_bot =  uy_bot / dx(2)
 
            uy_top = (u(i,j+1) - u(i,j))
            uy_top_wall = -(-16.0_dp_t * u(i,ny+1) + 20.0_dp_t * u(i,ny) &
                             -5.0_dp_t * u(i,ny-1) + u(i,ny-2) ) * 0.2_dp_t
            uy_top = merge(uy_top_wall, uy_top, j .eq. ny .and. &
                            (bc(2,2) .eq. BC_DIR) )
            uy_top = uy_top / dx(2)
 
            diff(i,j) = (ux_rgt-ux_lft)/dx(1) + (uy_top-uy_bot)/dx(2)
      end do
      end do

     end subroutine laplac_2d

      subroutine laplac_3d(u,diff,dx,ng,bc)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(in   ) :: u(1-ng:,1-ng:,1-ng:)
      real(kind=dp_t), intent(inout) :: diff(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      real(kind=dp_t) ::  ux_lft, ux_lft_wall
      real(kind=dp_t) ::  ux_rgt, ux_rgt_wall
      real(kind=dp_t) ::  uy_bot, uy_bot_wall
      real(kind=dp_t) ::  uy_top, uy_top_wall
      real(kind=dp_t) ::  uz_bot, uz_bot_wall
      real(kind=dp_t) ::  uz_top, uz_top_wall
      integer         ::  i,j,k,nx,ny,nz

! NOTE: This assumes Dirichlet boundary conditions on walls and inflow.
!       The value is imposed at the face, not the ghost cell.
      nx = size(diff,dim=1)
      ny = size(diff,dim=2)
      nz = size(diff,dim=3)

      do k = 0,nz-1
      do j = 0,ny-1
      do i = 0,nx-1
            ux_lft = (u(i,j,k) - u(i-1,j,k))
            ux_lft_wall = (-16.0_dp_t * u(0,j,k) + 20.0_dp_t * u(1,j,k) &
                            -5.0_dp_t * u(2,j,k) + u(3,j,k) ) * 0.2_dp_t
            ux_lft = merge(ux_lft_wall, ux_lft, i .eq. 1 .and. &
                          (bc(1,1) .eq. BC_DIR)) 
            ux_lft = ux_lft / dx(1)
 
            ux_rgt = (u(i+1,j,k) - u(i,j,k))
            ux_rgt_wall = -(-16.0_dp_t * u(nx+1,j,k) + 20.0_dp_t * u(nx,j,k) &
                             -5.0_dp_t * u(nx-1,j,k) + u(nx-2,j,k) ) * 0.2_dp_t
            ux_rgt = merge(ux_rgt_wall, ux_rgt, i .eq. nx .and. &
                          (bc(1,2) .eq. BC_DIR)) 
            ux_rgt = ux_rgt / dx(1)
 
            uy_bot = (u(i,j,k) - u(i,j-1,k))
            uy_bot_wall = (-16.0_dp_t * u(i,0,k) + 20.0_dp_t * u(i,1,k) &
                            -5.0_dp_t * u(i,2,k) + u(i,3,k) ) * 0.2_dp_t
            uy_bot = merge(uy_bot_wall, uy_bot, j .eq. 1 .and. &
                          (bc(2,1) .eq. BC_DIR)) 
            uy_bot =  uy_bot / dx(2)
 
            uy_top = (u(i,j+1,k) - u(i,j,k))
            uy_top_wall = -(-16.0_dp_t * u(i,ny+1,k) + 20.0_dp_t * u(i,ny,k) &
                             -5.0_dp_t * u(i,ny-1,k) + u(i,ny-2,k) ) * 0.2_dp_t
            uy_top = merge(uy_top_wall, uy_top, j .eq. ny .and. &
                          (bc(2,2) .eq. BC_DIR)) 
            uy_top = uy_top / dx(2)
 
            uz_bot = (u(i,j,k) - u(i,j,k-1))
            uz_bot_wall = (-16.0_dp_t * u(i,j,0) + 20.0_dp_t * u(i,j,1) &
                            -5.0_dp_t * u(i,j,2) + u(i,j,3) ) * 0.2_dp_t
            uz_bot = merge(uz_bot_wall, uz_bot, k .eq. 1 .and. &
                          (bc(3,1) .eq. BC_DIR)) 
            uz_bot =  uz_bot / dx(3)
 
            uz_top = (u(i,j,k+1) - u(i,j,k))
            uz_top_wall = -(-16.0_dp_t * u(i,j,nz+1) + 20.0_dp_t * u(i,j,nz) &
                             -5.0_dp_t * u(i,j,nz-1) + u(i,j,nz-2) ) * 0.2_dp_t
            uz_top = merge(uz_top_wall, uz_top, k .eq. nz .and. &
                          (bc(3,2) .eq. BC_DIR)) 
            uz_top = uz_top / dx(3)
 
            diff(i,j,k) = (ux_rgt-ux_lft)/dx(1) + (uy_top-uy_bot)/dx(2) &
                         +(uz_top-uz_bot)/dx(3) 
      end do
      end do
      end do

     end subroutine laplac_3d

end module laplac_module
