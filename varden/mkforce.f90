module mkforce_module

  use bl_types
  use laplac_module

  implicit none

contains

      subroutine mkforce_2d(force,ext_force,gp,rho,u,ng_u,ng_rho,dx,&
                            norm_vel_bc,tang_vel_bc,visc_coef,visc_fac)

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(  out) :: force(0:,0:,:)
      real(kind=dp_t), intent(in   ) :: ext_force(0:,0:,:)
      real(kind=dp_t), intent(in   ) ::    gp(0:,0:,:)
      real(kind=dp_t), intent(in   ) ::   rho(1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(in   ) ::     u(1-ng_u  :,1-ng_u  :,:)
      real(kind=dp_t), intent(in   ) ::    dx(2)
      integer        , intent(in   ) :: norm_vel_bc(:,:)
      integer        , intent(in   ) :: tang_vel_bc(:,:)
      real(kind=dp_t), intent(in   ) :: visc_coef, visc_fac

      real(kind=dp_t) :: diff(1:size(force,dim=1)-2,1:size(force,dim=2)-2)
      real(kind=dp_t) :: lapu
      integer :: i,j,n

!     Note: we can get away with sending only norm_vel_bc because we
!           know norm_vel_bc and tang_vel_bc are the same if viscous flow
 

      force = 0.0_dp_t
      do n = 1,2
         diff = 0.0_dp_t
         if (visc_coef * visc_fac > 0.0_dp_t) &
            call laplac_2d(u(:,:,n),diff,dx,ng_u,norm_vel_bc)
         do j = 1,size(force,dim=2)-2
         do i = 1,size(force,dim=1)-2
           lapu = visc_coef * visc_fac * diff(i,j)
           force(i,j,n) = ext_force(i,j,n) + (lapu - gp(i,j,n)) / rho(i,j)
         end do
         end do
      end do

      end subroutine mkforce_2d

      subroutine mkforce_3d(force,ext_force,gp,rho,u,ng_u,ng_rho,dx, &
                            norm_vel_bc,tang_vel_bc,visc_coef,visc_fac)

      integer        , intent(in   ) :: ng_u,ng_rho
      real(kind=dp_t), intent(  out) :: force(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) :: ext_force(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) ::    gp(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) ::   rho(1-ng_rho:,1-ng_rho:,1-ng_rho:)
      real(kind=dp_t), intent(in   ) ::     u(1-ng_u  :,1-ng_u  :,1-ng_u  :,:)
      real(kind=dp_t), intent(in   ) ::    dx(3)
      integer        , intent(in   ) :: norm_vel_bc(:,:)
      integer        , intent(in   ) :: tang_vel_bc(:,:)
      real(kind=dp_t), intent(in   ) :: visc_coef, visc_fac

      real(kind=dp_t) :: diff(1:size(force,dim=1)-2,1:size(force,dim=2)-2,&
                              1:size(force,dim=3)-2)
      real(kind=dp_t) :: lapu
      integer :: i,j,k,n

!     Note: we can get away with sending only norm_vel_bc because we
!           know norm_vel_bc and tang_vel_bc are the same if viscous flow

      force = 0.0_dp_t
      do n = 1,3
         diff = 0.0_dp_t
         if (visc_coef * visc_fac > 0.0_dp_t) &
            call laplac_3d(u(:,:,:,n),diff,dx,ng_u,norm_vel_bc)
         do k = 1,size(force,dim=3)-2
         do j = 1,size(force,dim=2)-2
         do i = 1,size(force,dim=1)-2
           lapu = visc_coef * visc_fac * diff(i,j,k)
           force(i,j,k,n) = ext_force(i,j,k,n) + (lapu - gp(i,j,k,n)) / rho(i,j,k)
         end do
         end do
         end do
      end do

      end subroutine mkforce_3d

      subroutine mkscalforce_2d(force,ext_force,s,ng,dx,bc,diff_coef,visc_fac)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(  out) :: force(0:,0:,:)
      real(kind=dp_t), intent(in   ) :: ext_force(0:,0:,:)
      real(kind=dp_t), intent(in   ) ::     s(1-ng  :,1-ng  :,:)
      real(kind=dp_t), intent(in   ) ::    dx(:)
      integer        , intent(in   ) :: bc(:,:) 
      real(kind=dp_t), intent(in   ) :: diff_coef, visc_fac

      real(kind=dp_t) :: diff(1:size(force,dim=1)-2,1:size(force,dim=2)-2)
      real(kind=dp_t) :: lapu
      integer :: i,j,n
 
      force = 0.0_dp_t
       diff = 0.0_dp_t
      do n = 1,size(s,dim=3)
         if (n > 1 .and. diff_coef * visc_fac > 0.0_dp_t) &
            call laplac_2d(s(:,:,n),diff,dx,ng,bc)
         do j = 1,size(force,dim=2)-2
         do i = 1,size(force,dim=1)-2
           lapu = diff_coef * visc_fac * diff(i,j)
           force(i,j,n) = ext_force(i,j,n) + lapu
         end do
         end do
      end do

      end subroutine mkscalforce_2d

      subroutine mkscalforce_3d(force,ext_force,s,ng,dx,bc,diff_coef,visc_fac)

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(  out) :: force(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) :: ext_force(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) :: s(1-ng  :,1-ng  :,1-ng  :,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:) 
      real(kind=dp_t), intent(in   ) :: diff_coef, visc_fac

      real(kind=dp_t) :: diff(1:size(force,dim=1)-2,1:size(force,dim=2)-2,&
                              1:size(force,dim=3)-2)
      real(kind=dp_t) :: lapu
      integer :: i,j,k,n

      force = 0.0_dp_t
       diff = 0.0_dp_t
      do n = 1,size(s,dim=4)
         diff = 0.0_dp_t
         if (n > 1 .and. diff_coef * visc_fac > 0.0_dp_t) &
            call laplac_3d(s(:,:,:,n),diff,dx,ng,bc)
         do k = 1,size(force,dim=3)-2
         do j = 1,size(force,dim=2)-2
         do i = 1,size(force,dim=1)-2
           lapu = diff_coef * visc_fac * diff(i,j,k)
           force(i,j,k,n) = ext_force(i,j,k,n) + lapu
         end do
         end do
         end do
      end do

      end subroutine mkscalforce_3d

end module mkforce_module
