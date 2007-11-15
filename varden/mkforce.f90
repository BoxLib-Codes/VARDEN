module mkforce_module

  use bl_types
  use laplac_module

  implicit none

contains

  subroutine mkvelforce_2d(vel_force,ext_vel_force,gp,rho,u,ng_u,ng_rho,dx,&
                           bc,visc_coef,visc_fac)

    integer        , intent(in   ) :: ng_u,ng_rho
    real(kind=dp_t), intent(  out) :: vel_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: gp(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho(1-ng_rho:,1-ng_rho:)
    real(kind=dp_t), intent(in   ) :: u(1-ng_u:,1-ng_u:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: visc_coef, visc_fac
    
    real(kind=dp_t) :: diff(1:size(vel_force,dim=1)-2,1:size(vel_force,dim=2)-2)
    real(kind=dp_t) :: lapu
    integer :: i,j,n
    
    vel_force = 0.0_dp_t
    do n = 1,2
       diff = 0.0_dp_t
       if (visc_coef * visc_fac > 0.0_dp_t) &
            call laplac_2d(u(:,:,n),diff,dx,ng_u,bc(:,:,n))
       do j = 1,size(vel_force,dim=2)-2
       do i = 1,size(vel_force,dim=1)-2
          lapu = visc_coef * visc_fac * diff(i,j)
          vel_force(i,j,n) = ext_vel_force(i,j,n) + (lapu - gp(i,j,n)) / rho(i,j)
       end do
       end do
    end do
    
  end subroutine mkvelforce_2d
  
  subroutine mkvelforce_3d(vel_force,ext_vel_force,gp,rho,u,ng_u,ng_rho,dx, &
                           bc,visc_coef,visc_fac)
    
    integer        , intent(in   ) :: ng_u,ng_rho
    real(kind=dp_t), intent(  out) :: vel_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: gp(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho(1-ng_rho:,1-ng_rho:,1-ng_rho:)
    real(kind=dp_t), intent(in   ) :: u(1-ng_u:,1-ng_u:,1-ng_u:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)
    real(kind=dp_t), intent(in   ) :: visc_coef, visc_fac
    
    real(kind=dp_t) :: diff(1:size(vel_force,dim=1)-2, &
                            1:size(vel_force,dim=2)-2, &
                            1:size(vel_force,dim=3)-2)
    
    real(kind=dp_t) :: lapu
    integer :: i,j,k,n
    
    vel_force = 0.0_dp_t
    do n = 1,3
       diff = 0.0_dp_t
       if (visc_coef * visc_fac > 0.0_dp_t) &
            call laplac_3d(u(:,:,:,n),diff,dx,ng_u,bc(:,:,n))
       do k = 1,size(vel_force,dim=3)-2
       do j = 1,size(vel_force,dim=2)-2
       do i = 1,size(vel_force,dim=1)-2
          lapu = visc_coef * visc_fac * diff(i,j,k)
          vel_force(i,j,k,n) = ext_vel_force(i,j,k,n) + (lapu - gp(i,j,k,n)) / rho(i,j,k)
       end do
       end do
       end do
    end do
    
  end subroutine mkvelforce_3d
  
  subroutine mkscalforce_2d(scal_force,ext_scal_force,s,ng,dx,bc,diff_coef,diff_fac)
    
    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: scal_force    (0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:) 
    real(kind=dp_t), intent(in   ) :: diff_coef, diff_fac
    
    real(kind=dp_t) :: diff(1:size(scal_force,dim=1)-2,1:size(scal_force,dim=2)-2)
    real(kind=dp_t) :: lapu
    integer :: i,j,n
    
    scal_force = 0.0_dp_t
    diff = 0.0_dp_t
    do n = 1,size(s,dim=3)
       if (n > 1 .and. diff_coef * diff_fac > 0.0_dp_t) &
            call laplac_2d(s(:,:,n),diff,dx,ng,bc(:,:,n))
       do j = 1,size(scal_force,dim=2)-2
       do i = 1,size(scal_force,dim=1)-2
          lapu = diff_coef * diff_fac * diff(i,j)
          scal_force(i,j,n) = ext_scal_force(i,j,n) + lapu
       end do
       end do
    end do
    
  end subroutine mkscalforce_2d
  
  subroutine mkscalforce_3d(scal_force,ext_scal_force,s,ng,dx,bc,diff_coef,diff_fac)
    
    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: scal_force    (0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,1-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:) 
    real(kind=dp_t), intent(in   ) :: diff_coef, diff_fac
    
    real(kind=dp_t) :: diff(1:size(scal_force,dim=1)-2, &
                            1:size(scal_force,dim=2)-2, &
                            1:size(scal_force,dim=3)-2)

    real(kind=dp_t) :: lapu
    integer :: i,j,k,n
    
    scal_force = 0.0_dp_t
    diff = 0.0_dp_t
    do n = 1,size(s,dim=4)
       diff = 0.0_dp_t
       if (n > 1 .and. diff_coef * diff_fac > 0.0_dp_t) &
            call laplac_3d(s(:,:,:,n),diff,dx,ng,bc(:,:,n))
       do k = 1,size(scal_force,dim=3)-2
       do j = 1,size(scal_force,dim=2)-2
       do i = 1,size(scal_force,dim=1)-2
          lapu = diff_coef * diff_fac * diff(i,j,k)
          scal_force(i,j,k,n) = ext_scal_force(i,j,k,n) + lapu
       end do
       end do
       end do
    end do
    
  end subroutine mkscalforce_3d
  
end module mkforce_module
