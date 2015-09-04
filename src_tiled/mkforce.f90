module mkforce_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use ml_restrict_fill_module
  use define_bc_module

  implicit none

  private

  public :: mkvelforce, mkscalforce

contains

  subroutine mkvelforce(mla,vel_force,ext_vel_force,s,gp,lapu,visc_fac,the_bc_tower)

    use probin_module, only: extrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_force(:)
    type(multifab) , intent(in   ) :: ext_vel_force(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: gp(:)
    type(multifab) , intent(in   ) :: lapu(:)
    real(kind=dp_t), intent(in   ) :: visc_fac
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,n,dm,nlevs,lo(mla%dim),hi(mla%dim)
    integer :: ng_f,ng_e,ng_g,ng_s,ng_l
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: lp(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: gpp(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt
    call build(bpt,"mkvelforce")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_f = vel_force(1)%ng
    ng_e = ext_vel_force(1)%ng
    ng_g = gp(1)%ng
    ng_s = s(1)%ng
    ng_l = lapu(1)%ng

    do n=1,nlevs
       call setval(vel_force(n),ZERO,all=.true.)
    end do


    !$omp parallel private(mfi,n,i,tilebox,tlo,thi) &
    !$omp private(fp,ep,gpp,sp,lp,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi,vel_force(n),tiling=.true.)

       do while(next_tile(mfi,i))
          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          fp  => dataptr(vel_force(n),i)
          ep  => dataptr(ext_vel_force(n),i)
          gpp => dataptr(gp(n),i)
          sp  => dataptr(s(n),i)
          lp => dataptr(lapu(n),i)
          lo = lwb(get_box(vel_force(n),i))
          hi = upb(get_box(vel_force(n),i))
          select case (dm)
          case (2)
             call mkvelforce_2d(fp(:,:,1,:), ep(:,:,1,:), gpp(:,:,1,:), &
                                sp(:,:,1,:), lp(:,:,1,:), &
                                ng_f,ng_e,ng_g,ng_s,ng_l, visc_fac, lo,hi,tlo,thi)
          case (3)
             call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), gpp(:,:,:,:), &
                                sp(:,:,:,:), lp(:,:,:,:), &
                                ng_f,ng_e,ng_g,ng_s,ng_l, visc_fac, lo,hi,tlo,thi)
          end select
       end do
    enddo
    !$omp end parallel

    call ml_restrict_and_fill(nlevs, vel_force, mla%mba%rr, the_bc_tower%bc_tower_array, &
         icomp=1,bcomp=extrap_comp,nc=dm,same_boundary=.true.)

    call destroy(bpt)

  end subroutine mkvelforce

  subroutine mkvelforce_2d(vel_force,ext_vel_force,gp,s,lapu, &
                           ng_f,ng_e,ng_g,ng_s,ng_l,visc_fac,glo,ghi,tlo,thi)

    use probin_module, only: boussinesq, visc_coef
 
    integer        , intent(in   ) :: glo(:),ghi(:),ng_f,ng_e,ng_g,ng_s,ng_l
    integer        , intent(in   ) :: tlo(:),thi(:)
    real(kind=dp_t), intent(  out) ::     vel_force(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(glo(1)-ng_e:,glo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::            gp(glo(1)-ng_g:,glo(2)-ng_g:,:)
    real(kind=dp_t), intent(in   ) ::             s(glo(1)-ng_s:,glo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::          lapu(glo(1)-ng_l:,glo(2)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local(2)
    integer :: i,j

       if (boussinesq .eq. 1) then
          do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             vel_force(i,j,1:2) = s(i,j,2) * ext_vel_force(i,j,1:2)
          end do
          end do
       else
          do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             vel_force(i,j,1:2) = ext_vel_force(i,j,1:2)
          end do
          end do
       end if

       do j = tlo(2), thi(2)
       do i = tlo(1), thi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,j,1:2)
          vel_force(i,j,1:2) = vel_force(i,j,1:2) + (lapu_local(1:2) - gp(i,j,1:2)) / s(i,j,1)
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       if (tlo(1) .eq. glo(1)) then
       do j = tlo(2), thi(2)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(tlo(1),j,1:2)
          vel_force(tlo(1)-1,j,1:2) = ext_vel_force(tlo(1)-1,j,1:2) &
               + (lapu_local(1:2) - gp(tlo(1)-1,j,1:2)) / s(tlo(1)-1,j,1)
       enddo
       end if

       if (thi(1) .eq. ghi(1)) then
       do j = tlo(2), thi(2)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(thi(1),j,1:2)
          vel_force(thi(1)+1,j,1:2) = ext_vel_force(thi(1)+1,j,1:2) &
               + (lapu_local(1:2) - gp(thi(1)+1,j,1:2)) / s(thi(1)+1,j,1)
       enddo
       end if

       if (tlo(2) .eq. glo(2)) then
       do i = tlo(1), thi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,tlo(2),1:2)
          vel_force(i,tlo(2)-1,1:2) = ext_vel_force(i,tlo(2)-1,1:2) &
               + (lapu_local(1:2) - gp(i,tlo(2)-1,1:2)) / s(i,tlo(2)-1,1)
       enddo
       end if

       if (thi(2) .eq. ghi(2)) then
       do i = tlo(1), thi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,thi(2),1:2)
          vel_force(i,thi(2)+1,1:2) = ext_vel_force(i,thi(2)+1,1:2) &
               + (lapu_local(1:2) - gp(i,thi(2)+1,1:2)) / s(i,thi(2)+1,1)
       enddo
       end if

  end subroutine mkvelforce_2d

  subroutine mkvelforce_3d(vel_force,ext_vel_force,gp,s,lapu, &
                           ng_f,ng_e,ng_g,ng_s,ng_l,visc_fac,glo,ghi,tlo,thi)

    use probin_module, only: boussinesq, visc_coef

    integer        , intent(in   ) :: glo(:),ghi(:),ng_f,ng_e,ng_g,ng_s,ng_l
    integer        , intent(in   ) :: tlo(:),thi(:)
    real(kind=dp_t), intent(  out) ::     vel_force(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::            gp(glo(1)-ng_g:,glo(2)-ng_g:,glo(3)-ng_g:,:)
    real(kind=dp_t), intent(in   ) ::             s(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::          lapu(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local(3)
    integer :: i,j,k

       if (boussinesq .eq. 1) then
          do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             vel_force(i,j,k,1:3) = s(i,j,k,2) * ext_vel_force(i,j,k,1:3)
          end do
          end do
          end do
       else
          do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             vel_force(i,j,k,1:3) = ext_vel_force(i,j,k,1:3)
          end do
          end do
          end do
       end if
       do k = tlo(3), thi(3)
       do j = tlo(2), thi(2)
       do i = tlo(1), thi(1)
          lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,k,1:3)
          vel_force(i,j,k,1:3) = vel_force(i,j,k,1:3) + &
               (lapu_local(1:3) - gp(i,j,k,1:3)) / s(i,j,k,1)
       end do
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       if (tlo(1) .eq. glo(1)) then
       do k=tlo(3),thi(3)
          do j=tlo(2),thi(2)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(tlo(1),j,k,1:3)
             vel_force(tlo(1)-1,j,k,1:3) = ext_vel_force(tlo(1)-1,j,k,1:3) + &
                  (lapu_local(1:3) - gp(tlo(1)-1,j,k,1:3)) / s(tlo(1)-1,j,k,1)
          enddo
       enddo
       end if

       if (thi(1) .eq. ghi(1)) then
       do k=tlo(3),thi(3)
          do j=tlo(2),thi(2)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(thi(1),j,k,1:3)
             vel_force(thi(1)+1,j,k,1:3) = ext_vel_force(thi(1)+1,j,k,1:3) + &
                  (lapu_local(1:3) - gp(thi(1)+1,j,k,1:3)) / s(thi(1)+1,j,k,1)
          enddo
       enddo      
       end if
       
       if (tlo(2) .eq. glo(2)) then
       do k=tlo(3),thi(3)
          do i=tlo(1),thi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,tlo(2),k,1:3)
             vel_force(i,tlo(2)-1,k,1:3) = ext_vel_force(i,tlo(2)-1,k,1:3) + &
                  (lapu_local(1:3) - gp(i,tlo(2)-1,k,1:3)) / s(i,tlo(2)-1,k,1)
          enddo
       enddo
       end if

       if (thi(2) .eq. ghi(2)) then
       do k=tlo(3),thi(3)
          do i=tlo(1),thi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,thi(2),k,1:3)
             vel_force(i,thi(2)+1,k,1:3) = ext_vel_force(i,thi(2)+1,k,1:3) + &
                  (lapu_local(1:3) - gp(i,thi(2)+1,k,1:3)) / s(i,thi(2)+1,k,1)
          enddo
       enddo
       end if

       if (tlo(3) .eq. glo(3)) then
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,tlo(3),1:3)
             vel_force(i,j,tlo(3)-1,1:3) = ext_vel_force(i,j,tlo(3)-1,1:3) + &
                  (lapu_local(1:3) - gp(i,j,tlo(3)-1,1:3)) / s(i,j,tlo(3)-1,1)
          enddo
       enddo
       end if

       if (thi(3) .eq. ghi(3)) then
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,thi(3),1:3)
             vel_force(i,j,thi(3)+1,1:3) = ext_vel_force(i,j,thi(3)+1,1:3) + &
                  (lapu_local(1:3) - gp(i,j,thi(3)+1,1:3)) / s(i,j,thi(3)+1,1)
          enddo
       enddo
       end if

  end subroutine mkvelforce_3d

  subroutine mkscalforce(mla,scal_force,ext_scal_force,laps,diff_fac,the_bc_tower)

    use probin_module, only: nscal,extrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    type(multifab) , intent(in   ) :: laps(:)
    real(kind=dp_t), intent(in   ) :: diff_fac
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,n,ng_f,ng_e,ng_l,dm,nlevs,lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: lp(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"mkscalforce")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_f = scal_force(1)%ng
    ng_e = ext_scal_force(1)%ng
    ng_l = laps(1)%ng

    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do


    !$omp parallel private(mfi,n,i,tilebox,tlo,thi) &
    !$omp private(fp,ep,lp,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi,scal_force(n),tiling=.true.)

       do while(next_tile(mfi,i))
          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          fp => dataptr(scal_force(n),i)
          ep => dataptr(ext_scal_force(n),i)
          lp => dataptr(laps(n),i)
          lo = lwb(get_box(scal_force(n),i))
          hi = upb(get_box(scal_force(n),i))
          select case (dm)
          case (2)
             call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), lp(:,:,1,:), ng_f, ng_e, ng_l, diff_fac, &
                  lo, hi, tlo, thi)
          case (3)
             call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), lp(:,:,:,:), ng_f, ng_e, ng_l, diff_fac, &
                  lo, hi, tlo, thi)
          end select
       end do
    enddo
    !$omp end parallel

    call ml_restrict_and_fill(nlevs, scal_force, mla%mba%rr, the_bc_tower%bc_tower_array, &
         icomp=1,bcomp=extrap_comp,nc=nscal,same_boundary=.true.)

    call destroy(bpt)

  end subroutine mkscalforce

  subroutine mkscalforce_2d(scal_force,ext_scal_force,laps,ng_f,ng_e,ng_l,diff_fac, &
       glo,ghi,tlo,thi)

    use probin_module, only : nscal, diff_coef

    integer        , intent(in   ) :: ng_f,ng_e,ng_l,glo(:),ghi(:),tlo(:),thi(:)
    real(kind=dp_t), intent(  out) ::     scal_force(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(glo(1)-ng_e:,glo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::           laps(glo(1)-ng_l:,glo(2)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local(2:nscal)
    integer :: i,j


    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do j = tlo(2), thi(2)
       do i = tlo(1), thi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,2:nscal)
          scal_force(i,j,2:nscal) = ext_scal_force(i,j,2:nscal) + laps_local(2:nscal)
       enddo
    enddo
    
    ! we use 0th order extrapolation for laplacian term in ghost cells
    if (tlo(1) .eq. glo(1)) then
       do j = tlo(2), thi(2)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(tlo(1),j,2:nscal)
          scal_force(tlo(1)-1,j,2:nscal) = ext_scal_force(tlo(1)-1,j,2:nscal) + laps_local(2:nscal)
       enddo
    end if

    if (thi(1) .eq. ghi(1)) then
       do j = tlo(2), thi(2)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(thi(1),j,2:nscal)
          scal_force(thi(1)+1,j,2:nscal) = ext_scal_force(thi(1)+1,j,2:nscal) + laps_local(2:nscal)
       enddo
    end if

    if (tlo(2) .eq. glo(2)) then
       do i = tlo(1), thi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,tlo(2),2:nscal)
          scal_force(i,tlo(2)-1,2:nscal) = ext_scal_force(i,tlo(2)-1,2:nscal) + laps_local(2:nscal)
       enddo
    end if

    if (thi(2) .eq. ghi(2)) then
       do i = tlo(1), thi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,thi(2),2:nscal)
          scal_force(i,thi(2)+1,2:nscal) = ext_scal_force(i,thi(2)+1,2:nscal) + laps_local(2:nscal)
       enddo
    end if

  end subroutine mkscalforce_2d

  subroutine mkscalforce_3d(scal_force,ext_scal_force,laps,ng_f,ng_e,ng_l,diff_fac, &
       glo,ghi,tlo,thi)

    use probin_module, only : nscal, diff_coef

    integer        , intent(in   ) :: ng_f,ng_e,ng_l,glo(:),ghi(:),tlo(:),thi(:)
    real(kind=dp_t), intent(  out) ::     scal_force(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::           laps(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local(2:nscal)
    integer :: i,j,k


    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do k = tlo(3), thi(3)
       do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,k,2:nscal)
             scal_force(i,j,k,2:nscal) = ext_scal_force(i,j,k,2:nscal) + laps_local(2:nscal)
          end do
       end do
    end do
    
    ! we use 0th order extrapolation for laplacian term in ghost cells
    if (tlo(1) .eq. glo(1)) then
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,tlo(3),2:nscal)
             scal_force(i,j,tlo(3)-1,2:nscal) = ext_scal_force(i,j,tlo(3)-1,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

    if (thi(1) .eq. ghi(1)) then
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,thi(3),2:nscal)
             scal_force(i,j,thi(3)+1,2:nscal) = ext_scal_force(i,j,thi(3)+1,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

    if (tlo(2) .eq. glo(2)) then
       do k=tlo(3),thi(3)
          do i=tlo(1),thi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,tlo(2),k,2:nscal)
             scal_force(i,tlo(2)-1,k,2:nscal) = ext_scal_force(i,tlo(2)-1,k,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

    if (thi(2) .eq. ghi(2)) then
       do k=tlo(3),thi(3)
          do i=tlo(1),thi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,thi(2),k,2:nscal)
             scal_force(i,thi(2)+1,k,2:nscal) = ext_scal_force(i,thi(2)+1,k,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

    if (tlo(3) .eq. glo(3)) then
       do k=tlo(3),thi(3)
          do j=tlo(2),thi(2)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(tlo(1),j,k,2:nscal)
             scal_force(tlo(1)-1,j,k,2:nscal) = ext_scal_force(tlo(1)-1,j,k,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

    if (thi(3) .eq. ghi(3)) then
       do k=tlo(3),thi(3)
          do j=tlo(2),thi(2)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(thi(1),j,k,2:nscal)
             scal_force(thi(1)+1,j,k,2:nscal) = ext_scal_force(thi(1)+1,j,k,2:nscal) + laps_local(2:nscal)
          enddo
       enddo
    end if

  end subroutine mkscalforce_3d

end module mkforce_module
