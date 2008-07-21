module mkforce_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: mkvelforce, mkscalforce, mksource

contains

  subroutine mkvelforce(nlevs,vel_force,ext_vel_force,s,gp,lapu,visc_fac)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: vel_force(:)
    type(multifab) , intent(in   ) :: ext_vel_force(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: gp(:)
    type(multifab) , intent(in   ) :: lapu(:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    ! local
    integer :: i,n,ng,dm
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: lp(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: gpp(:,:,:,:)

    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs
       do i = 1, vel_force(n)%nboxes
          if ( remote(vel_force(n),i) ) cycle
          fp  => dataptr(vel_force(n),i)
          ep  => dataptr(ext_vel_force(n),i)
          gpp => dataptr(gp(n),i)
          sp  => dataptr(s(n),i)
          lp => dataptr(lapu(n),i)
          select case (dm)
          case (2)
             call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), gpp(:,:,1,:), &
                                sp(:,:,1,:), lp(:,:,1,:), &
                                ng, visc_fac)
          case (3)
             call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), gpp(:,:,:,:), &
                                sp(:,:,:,:), lp(:,:,:,:), &
                                ng, visc_fac)
          end select
       end do

       call multifab_fill_boundary(vel_force(n))
    enddo

  end subroutine mkvelforce

  subroutine mkvelforce_2d(vel_force,ext_vel_force,gp,s,lapu,ng,visc_fac)

    use probin_module, only: boussinesq, visc_coef
 
    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: vel_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: gp(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,:)
    real(kind=dp_t), intent(in   ) :: lapu(1:,1:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local
    integer :: i,j,is,ie,js,je,n

    is = 1
    js = 1
    ie = size(vel_force,dim=1)-2
    je = size(vel_force,dim=2)-2

    do n = 1,2
       if (boussinesq .eq. 1) then
          do j = js, je
          do i = is, ie
             vel_force(i,j,n) = s(i,j,2) * ext_vel_force(i,j,n)
          end do
          end do
       else
          do j = js, je
          do i = is, ie
             vel_force(i,j,n) = ext_vel_force(i,j,n)
          end do
          end do
       end if

       do j = js, je
       do i = is, ie
          lapu_local = visc_coef * visc_fac * lapu(i,j,n)
          vel_force(i,j,n) = vel_force(i,j,n) + (lapu_local - gp(i,j,n)) / s(i,j,1)
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i = is, ie
          lapu_local = visc_coef * visc_fac * lapu(i,js,n)
          vel_force(i,js-1,n) = ext_vel_force(i,js-1,n) &
               + (lapu_local - gp(i,js-1,n)) / s(i,js-1,1)
       enddo
       do i = is, ie
          lapu_local = visc_coef * visc_fac * lapu(i,je,n)
          vel_force(i,je+1,n) = ext_vel_force(i,je+1,n) &
               + (lapu_local - gp(i,je+1,n)) / s(i,je+1,1)
       enddo
       do j = js, je
          lapu_local = visc_coef * visc_fac * lapu(is,j,n)
          vel_force(is-1,j,n) = ext_vel_force(is-1,j,n) &
               + (lapu_local - gp(is-1,j,n)) / s(is-1,j,1)
       enddo
       do j = js, je
          lapu_local = visc_coef * visc_fac * lapu(ie,j,n)
          vel_force(ie+1,j,n) = ext_vel_force(ie+1,j,n) &
               + (lapu_local - gp(ie+1,j,n)) / s(ie+1,j,1)
       enddo
    end do

  end subroutine mkvelforce_2d

  subroutine mkvelforce_3d(vel_force,ext_vel_force,gp,s,lapu,ng,visc_fac)

    use probin_module, only: boussinesq, visc_coef

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: vel_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: gp(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,1-ng:,:)
    real(kind=dp_t), intent(in   ) :: lapu(1:,1:,1:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local
    integer :: i,j,k,is,ie,js,je,ks,ke,n

    is = 1
    js = 1
    ks = 1
    ie = size(vel_force,dim=1)-2
    je = size(vel_force,dim=2)-2
    ke = size(vel_force,dim=3)-2

    do n = 1,3
       if (boussinesq .eq. 1) then
          do k = ks, ke
          do j = js, je
          do i = is, ie
             vel_force(i,j,k,n) = s(i,j,k,2) * ext_vel_force(i,j,k,n)
          end do
          end do
          end do
       else
          do k = ks, ke
          do j = js, je
          do i = is, ie
             vel_force(i,j,k,n) = ext_vel_force(i,j,k,n)
          end do
          end do
          end do
       end if
       do k = ks, ke
       do j = js, je
       do i = is, ie
          lapu_local = visc_coef * visc_fac * lapu(i,j,k,n)
          vel_force(i,j,k,n) = vel_force(i,j,k,n) + &
               (lapu_local - gp(i,j,k,n)) / s(i,j,k,1)
       end do
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i=is,ie
          do j=js,je
             lapu_local = visc_coef * visc_fac * lapu(i,j,ks,n)
             vel_force(i,j,ks-1,n) = ext_vel_force(i,j,ks-1,n) + &
                  (lapu_local - gp(i,j,ks-1,n)) / s(i,j,ks-1,1)
          enddo
       enddo
       do i=is,ie
          do j=js,je
             lapu_local = visc_coef * visc_fac * lapu(i,j,ke,n)
             vel_force(i,j,ke+1,n) = ext_vel_force(i,j,ke+1,n) + &
                  (lapu_local - gp(i,j,ke+1,n)) / s(i,j,ke+1,1)
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             lapu_local = visc_coef * visc_fac * lapu(i,js,k,n)
             vel_force(i,js-1,k,n) = ext_vel_force(i,js-1,k,n) + &
                  (lapu_local - gp(i,js-1,k,n)) / s(i,js-1,k,1)
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             lapu_local = visc_coef * visc_fac * lapu(i,je,k,n)
             vel_force(i,je+1,k,n) = ext_vel_force(i,je+1,k,n) + &
                  (lapu_local - gp(i,je+1,k,n)) / s(i,je+1,k,1)
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             lapu_local = visc_coef * visc_fac * lapu(is,j,k,n)
             vel_force(i,j,ks-1,n) = ext_vel_force(i,j,ks-1,n) + &
                  (lapu_local - gp(i,j,ks-1,n)) / s(i,j,ks-1,1)
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             lapu_local = visc_coef * visc_fac * lapu(ie,j,k,n)
             vel_force(i,j,ke+1,n) = ext_vel_force(i,j,ke+1,n) + &
                  (lapu_local - gp(i,j,ke+1,n)) / s(i,j,ke+1,1)
          enddo
       enddo
    end do

  end subroutine mkvelforce_3d

  subroutine mkscalforce(nlevs,scal_force,ext_scal_force,laps,diff_fac)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    type(multifab) , intent(in   ) :: laps(:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    ! local
    integer :: i,n,dm
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: lp(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)

    dm = scal_force(1)%dim

    do n=1,nlevs
       do i = 1, scal_force(n)%nboxes
          if ( remote(scal_force(n),i) ) cycle
          fp => dataptr(scal_force(n),i)
          lp => dataptr(laps(n),i)
          ep => dataptr(ext_scal_force(n),i)
          select case (dm)
          case (2)
             call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), lp(:,:,1,:), diff_fac)
          case (3)
             call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), lp(:,:,:,:), diff_fac)
          end select
       end do

       call multifab_fill_boundary(scal_force(n))
    enddo

  end subroutine mkscalforce

  subroutine mkscalforce_2d(scal_force,ext_scal_force,laps,diff_fac)

    use probin_module, only : nscal, diff_coef

    real(kind=dp_t), intent(  out) :: scal_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: laps(1:,1:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local
    integer :: i,j,is,ie,js,je,n

    is = 1
    js = 1
    ie = size(scal_force,dim=1)-2
    je = size(scal_force,dim=2)-2

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do n = 2,nscal

       do j = js, je
          do i = is, ie
             laps_local = diff_coef * diff_fac * laps(i,j,n)
             scal_force(i,j,n) = ext_scal_force(i,j,n) + laps_local
          enddo
       enddo

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i = is, ie
          laps_local = diff_coef * diff_fac * laps(i,js,n)
          scal_force(i,js-1,n) = ext_scal_force(i,js-1,n) + laps_local
       enddo
       do i = is, ie
          laps_local = diff_coef * diff_fac * laps(i,je,n)
          scal_force(i,je+1,n) = ext_scal_force(i,je+1,n) + laps_local
       enddo
       do j = js, je
          laps_local = diff_coef * diff_fac * laps(is,j,n)
          scal_force(is-1,j,n) = ext_scal_force(is-1,j,n) + laps_local
       enddo
       do j = js, je
          laps_local = diff_coef * diff_fac * laps(ie,j,n)
          scal_force(ie+1,j,n) = ext_scal_force(ie+1,j,n) + laps_local
       enddo

    end do

  end subroutine mkscalforce_2d

  subroutine mkscalforce_3d(scal_force,ext_scal_force,laps,diff_fac)

    use probin_module, only : nscal, diff_coef

    real(kind=dp_t), intent(  out) :: scal_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: laps(1:,1:,1:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local
    integer :: i,j,k,is,ie,js,je,ks,ke,n

    is = 1
    js = 1
    ks = 1
    ie = size(scal_force,dim=1)-2
    je = size(scal_force,dim=2)-2
    ke = size(scal_force,dim=3)-2

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do n = 2, nscal
       do k = ks, ke
          do j = js, je
             do i = is, ie
                laps_local = diff_coef * diff_fac * laps(i,j,k,n)
                scal_force(i,j,k,n) = ext_scal_force(i,j,k,n) + laps_local
             end do
          end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i=is,ie
          do j=js,je
             laps_local = diff_coef * diff_fac * laps(i,j,ks,n)
             scal_force(i,j,ks-1,n) = ext_scal_force(i,j,ks-1,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do j=js,je
             laps_local = diff_coef * diff_fac * laps(i,j,ke,n)
             scal_force(i,j,ke+1,n) = ext_scal_force(i,j,ke+1,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             laps_local = diff_coef * diff_fac * laps(i,js,k,n)
             scal_force(i,js-1,k,n) = ext_scal_force(i,js-1,k,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             laps_local = diff_coef * diff_fac * laps(i,je,k,n)
             scal_force(i,je+1,k,n) = ext_scal_force(i,je+1,k,n) + laps_local
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             laps_local = diff_coef * diff_fac * laps(is,j,k,n)
             scal_force(is-1,j,k,n) = ext_scal_force(is-1,j,k,n) + laps_local
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             laps_local = diff_coef * diff_fac * laps(ie,j,k,n)
             scal_force(ie+1,j,k,n) = ext_scal_force(ie+1,j,k,n) + laps_local
          enddo
       enddo
    end do

  end subroutine mkscalforce_3d

  subroutine mksource(nlevs,scal_force,ext_scal_force,laps,diff_fac)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    type(multifab) , intent(in   ) :: laps(:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    ! local
    integer :: i,n,dm
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: lp(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)

    dm = scal_force(1)%dim

    do n=1,nlevs
       do i = 1, scal_force(n)%nboxes
          if ( remote(scal_force(n),i) ) cycle
          fp => dataptr(scal_force(n),i)
          lp => dataptr(laps(n),i)
          ep => dataptr(ext_scal_force(n),i)
          select case (dm)
          case (2)
             call mksource_2d(fp(:,:,1,:), ep(:,:,1,:), lp(:,:,1,:), diff_fac)
          case (3)
             call mksource_3d(fp(:,:,:,:), ep(:,:,:,:), lp(:,:,:,:), diff_fac)
          end select
       end do

       call multifab_fill_boundary(scal_force(n))
    enddo

  end subroutine mksource

 subroutine mksource_2d(scal_force,ext_scal_force,laps,diff_fac)

    use probin_module, only : nscal

    real(kind=dp_t), intent(  out) :: scal_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,:)
    real(kind=dp_t), intent(in   ) :: laps(1:,1:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local
    integer :: i,j,is,ie,js,je,n

    is = 1
    js = 1
    ie = size(scal_force,dim=1)-2
    je = size(scal_force,dim=2)-2

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do n = 2,nscal

       do j = js, je
          do i = is, ie
             laps_local = diff_fac * laps(i,j,n-1)
             scal_force(i,j,n) = ext_scal_force(i,j,n) + laps_local
          enddo
       enddo

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i = is, ie
          laps_local = diff_fac * laps(i,js,n-1)
          scal_force(i,js-1,n) = ext_scal_force(i,js-1,n) + laps_local
       enddo
       do i = is, ie
          laps_local = diff_fac * laps(i,je,n-1)
          scal_force(i,je+1,n) = ext_scal_force(i,je+1,n) + laps_local
       enddo
       do j = js, je
          laps_local = diff_fac * laps(is,j,n-1)
          scal_force(is-1,j,n) = ext_scal_force(is-1,j,n) + laps_local
       enddo
       do j = js, je
          laps_local = diff_fac * laps(ie,j,n-1)
          scal_force(ie+1,j,n) = ext_scal_force(ie+1,j,n) + laps_local
       enddo

    end do

  end subroutine mksource_2d

  subroutine mksource_3d(scal_force,ext_scal_force,laps,diff_fac)

    use probin_module, only : nscal

    real(kind=dp_t), intent(  out) :: scal_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(0:,0:,0:,:)
    real(kind=dp_t), intent(in   ) :: laps(1:,1:,1:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local
    integer :: i,j,k,is,ie,js,je,ks,ke,n

    is = 1
    js = 1
    ks = 1
    ie = size(scal_force,dim=1)-2
    je = size(scal_force,dim=2)-2
    ke = size(scal_force,dim=3)-2

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do n = 2, nscal
       do k = ks, ke
          do j = js, je
             do i = is, ie
                laps_local = diff_fac * laps(i,j,k,n-1)
                scal_force(i,j,k,n) = ext_scal_force(i,j,k,n) + laps_local
             end do
          end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do i=is,ie
          do j=js,je
             laps_local = diff_fac * laps(i,j,ks,n-1)
             scal_force(i,j,ks-1,n) = ext_scal_force(i,j,ks-1,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do j=js,je
             laps_local = diff_fac * laps(i,j,ke,n-1)
             scal_force(i,j,ke+1,n) = ext_scal_force(i,j,ke+1,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             laps_local = diff_fac * laps(i,js,k,n-1)
             scal_force(i,js-1,k,n) = ext_scal_force(i,js-1,k,n) + laps_local
          enddo
       enddo
       do i=is,ie
          do k=ks,ke
             laps_local = diff_fac * laps(i,je,k,n-1)
             scal_force(i,je+1,k,n) = ext_scal_force(i,je+1,k,n) + laps_local
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             laps_local = diff_fac * laps(is,j,k,n-1)
             scal_force(is-1,j,k,n) = ext_scal_force(is-1,j,k,n) + laps_local
          enddo
       enddo
       do j=js,je
          do k=ks,ke
             laps_local = diff_fac * laps(ie,j,k,n-1)
             scal_force(ie+1,j,k,n) = ext_scal_force(ie+1,j,k,n) + laps_local
          enddo
       enddo
    end do

  end subroutine mksource_3d

end module mkforce_module
