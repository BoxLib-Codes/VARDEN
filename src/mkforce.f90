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
       do i = 1, nfabs(vel_force(n))
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
                                ng_f,ng_e,ng_g,ng_s,ng_l, visc_fac, lo, hi)
          case (3)
             call mkvelforce_3d(fp(:,:,:,:), ep(:,:,:,:), gpp(:,:,:,:), &
                                sp(:,:,:,:), lp(:,:,:,:), &
                                ng_f,ng_e,ng_g,ng_s,ng_l, visc_fac, lo, hi)
          end select
       end do
    enddo

    call ml_restrict_and_fill(nlevs, vel_force, mla%mba%rr, the_bc_tower%bc_tower_array, &
         icomp=1,bcomp=extrap_comp,nc=dm,same_boundary=.true.)

    call destroy(bpt)

  end subroutine mkvelforce

  subroutine mkvelforce_2d(vel_force,ext_vel_force,gp,s,lapu, &
                           ng_f,ng_e,ng_g,ng_s,ng_l,visc_fac,lo,hi)

    use probin_module, only: boussinesq, visc_coef
 
    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_e,ng_g,ng_s,ng_l
    real(kind=dp_t), intent(  out) ::     vel_force(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::            gp(lo(1)-ng_g:,lo(2)-ng_g:,:)
    real(kind=dp_t), intent(in   ) ::             s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::          lapu(lo(1)-ng_l:,lo(2)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local(2)
    integer :: i,j,n

       if (boussinesq .eq. 1) then
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel_force(i,j,1:2) = s(i,j,2) * ext_vel_force(i,j,1:2)
          end do
          end do
       else
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel_force(i,j,1:2) = ext_vel_force(i,j,1:2)
          end do
          end do
       end if

       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,j,1:2)
          vel_force(i,j,1:2) = vel_force(i,j,1:2) + (lapu_local(1:2) - gp(i,j,1:2)) / s(i,j,1)
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do j = lo(2), hi(2)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(lo(1),j,1:2)
          vel_force(lo(1)-1,j,1:2) = ext_vel_force(lo(1)-1,j,1:2) &
               + (lapu_local(1:2) - gp(lo(1)-1,j,1:2)) / s(lo(1)-1,j,1)
       enddo
       do j = lo(2), hi(2)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(hi(1),j,1:2)
          vel_force(hi(1)+1,j,1:2) = ext_vel_force(hi(1)+1,j,1:2) &
               + (lapu_local(1:2) - gp(hi(1)+1,j,1:2)) / s(hi(1)+1,j,1)
       enddo
       do i = lo(1), hi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,lo(2),1:2)
          vel_force(i,lo(2)-1,1:2) = ext_vel_force(i,lo(2)-1,1:2) &
               + (lapu_local(1:2) - gp(i,lo(2)-1,1:2)) / s(i,lo(2)-1,1)
       enddo
       do i = lo(1), hi(1)
          lapu_local(1:2) = visc_coef * visc_fac * lapu(i,hi(2),1:2)
          vel_force(i,hi(2)+1,1:2) = ext_vel_force(i,hi(2)+1,1:2) &
               + (lapu_local(1:2) - gp(i,hi(2)+1,1:2)) / s(i,hi(2)+1,1)
       enddo


  end subroutine mkvelforce_2d

  subroutine mkvelforce_3d(vel_force,ext_vel_force,gp,s,lapu, &
                           ng_f,ng_e,ng_g,ng_s,ng_l,visc_fac,lo,hi)

    use probin_module, only: boussinesq, visc_coef

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_e,ng_g,ng_s,ng_l
    real(kind=dp_t), intent(  out) ::     vel_force(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_vel_force(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::            gp(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
    real(kind=dp_t), intent(in   ) ::             s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) ::          lapu(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: visc_fac

    real(kind=dp_t) :: lapu_local(3)
    integer :: i,j,k,n

       if (boussinesq .eq. 1) then
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel_force(i,j,k,1:3) = s(i,j,k,2) * ext_vel_force(i,j,k,1:3)
          end do
          end do
          end do
       else
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             vel_force(i,j,k,1:3) = ext_vel_force(i,j,k,1:3)
          end do
          end do
          end do
       end if
       do k = lo(3), hi(3)
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,k,1:3)
          vel_force(i,j,k,1:3) = vel_force(i,j,k,1:3) + &
               (lapu_local(1:3) - gp(i,j,k,1:3)) / s(i,j,k,1)
       end do
       end do
       end do

       ! we use 0th order extrapolation for laplacian term in ghost cells
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(lo(1),j,k,1:3)
             vel_force(lo(1)-1,j,k,1:3) = ext_vel_force(lo(1)-1,j,k,1:3) + &
                  (lapu_local(1:3) - gp(lo(1)-1,j,k,1:3)) / s(lo(1)-1,j,k,1)
          enddo
       enddo

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(hi(1),j,k,1:3)
             vel_force(hi(1)+1,j,k,1:3) = ext_vel_force(hi(1)+1,j,k,1:3) + &
                  (lapu_local(1:3) - gp(hi(1)+1,j,k,1:3)) / s(hi(1)+1,j,k,1)
          enddo
       enddo      

       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,lo(2),k,1:3)
             vel_force(i,lo(2)-1,k,1:3) = ext_vel_force(i,lo(2)-1,k,1:3) + &
                  (lapu_local(1:3) - gp(i,lo(2)-1,k,1:3)) / s(i,lo(2)-1,k,1)
          enddo
       enddo

       do k=lo(3),hi(3)
          do i=lo(1),hi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,hi(2),k,1:3)
             vel_force(i,hi(2)+1,k,1:3) = ext_vel_force(i,hi(2)+1,k,1:3) + &
                  (lapu_local(1:3) - gp(i,hi(2)+1,k,1:3)) / s(i,hi(2)+1,k,1)
          enddo
       enddo

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,lo(3),1:3)
             vel_force(i,j,lo(3)-1,1:3) = ext_vel_force(i,j,lo(3)-1,1:3) + &
                  (lapu_local(1:3) - gp(i,j,lo(3)-1,1:3)) / s(i,j,lo(3)-1,1)
          enddo
       enddo

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             lapu_local(1:3) = visc_coef * visc_fac * lapu(i,j,hi(3),1:3)
             vel_force(i,j,hi(3)+1,1:3) = ext_vel_force(i,j,hi(3)+1,1:3) + &
                  (lapu_local(1:3) - gp(i,j,hi(3)+1,1:3)) / s(i,j,hi(3)+1,1)
          enddo
       enddo

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

    type(bl_prof_timer), save :: bpt

    call build(bpt,"mkscalforce")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_f = scal_force(1)%ng
    ng_e = ext_scal_force(1)%ng
    ng_l = laps(1)%ng

    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
       do i = 1, nfabs(scal_force(n))
          fp => dataptr(scal_force(n),i)
          ep => dataptr(ext_scal_force(n),i)
          lp => dataptr(laps(n),i)
          lo = lwb(get_box(scal_force(n),i))
          hi = upb(get_box(scal_force(n),i))
          select case (dm)
          case (2)
             call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), lp(:,:,1,:), ng_f, ng_e, ng_l, diff_fac, lo, hi)
          case (3)
             call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), lp(:,:,:,:), ng_f, ng_e, ng_l, diff_fac, lo, hi)
          end select
       end do
    enddo

    call ml_restrict_and_fill(nlevs, scal_force, mla%mba%rr, the_bc_tower%bc_tower_array, &
         icomp=1,bcomp=extrap_comp,nc=nscal,same_boundary=.true.)

    call destroy(bpt)

  end subroutine mkscalforce

  subroutine mkscalforce_2d(scal_force,ext_scal_force,laps,ng_f,ng_e,ng_l,diff_fac,lo,hi)

    use probin_module, only : nscal, diff_coef

    integer        , intent(in   ) :: ng_f,ng_e,ng_l,lo(:),hi(:)
    real(kind=dp_t), intent(  out) ::     scal_force(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::           laps(lo(1)-ng_l:,lo(2)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local(2:nscal)
    integer :: i,j,n

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,2:nscal)
          scal_force(i,j,2:nscal) = ext_scal_force(i,j,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    ! we use 0th order extrapolation for laplacian term in ghost cells
    do j = lo(2), hi(2)
       laps_local(2:nscal) = diff_coef * diff_fac * laps(lo(1),j,2:nscal)
       scal_force(lo(1)-1,j,2:nscal) = ext_scal_force(lo(1)-1,j,2:nscal) + laps_local(2:nscal)
    enddo
    do j = lo(2), hi(2)
       laps_local(2:nscal) = diff_coef * diff_fac * laps(hi(1),j,2:nscal)
       scal_force(hi(1)+1,j,2:nscal) = ext_scal_force(hi(1)+1,j,2:nscal) + laps_local(2:nscal)
    enddo 
   do i = lo(1), hi(1)
       laps_local(2:nscal) = diff_coef * diff_fac * laps(i,lo(2),2:nscal)
       scal_force(i,lo(2)-1,2:nscal) = ext_scal_force(i,lo(2)-1,2:nscal) + laps_local(2:nscal)
    enddo
    do i = lo(1), hi(1)
       laps_local(2:nscal) = diff_coef * diff_fac * laps(i,hi(2),2:nscal)
       scal_force(i,hi(2)+1,2:nscal) = ext_scal_force(i,hi(2)+1,2:nscal) + laps_local(2:nscal)
    enddo

  end subroutine mkscalforce_2d

  subroutine mkscalforce_3d(scal_force,ext_scal_force,laps,ng_f,ng_e,ng_l,diff_fac,lo,hi)

    use probin_module, only : nscal, diff_coef

    integer        , intent(in   ) :: ng_f,ng_e,ng_l,lo(:),hi(:)
    real(kind=dp_t), intent(  out) ::     scal_force(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: ext_scal_force(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t), intent(in   ) ::           laps(lo(1)-ng_l:,lo(2)-ng_l:,lo(3)-ng_l:,:)
    real(kind=dp_t), intent(in   ) :: diff_fac

    real(kind=dp_t) :: laps_local(2:nscal)
    integer :: i,j,k,n

    scal_force = 0.0_dp_t

    ! NOTE: component 1 is density which doesn't diffuse, so we start with component 2 
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,k,2:nscal)
             scal_force(i,j,k,2:nscal) = ext_scal_force(i,j,k,2:nscal) + laps_local(2:nscal)
          end do
       end do
    end do

    ! we use 0th order extrapolation for laplacian term in ghost cells

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,lo(3),2:nscal)
          scal_force(i,j,lo(3)-1,2:nscal) = ext_scal_force(i,j,lo(3)-1,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,j,hi(3),2:nscal)
          scal_force(i,j,hi(3)+1,2:nscal) = ext_scal_force(i,j,hi(3)+1,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,lo(2),k,2:nscal)
          scal_force(i,lo(2)-1,k,2:nscal) = ext_scal_force(i,lo(2)-1,k,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(i,hi(2),k,2:nscal)
          scal_force(i,hi(2)+1,k,2:nscal) = ext_scal_force(i,hi(2)+1,k,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(lo(1),j,k,2:nscal)
          scal_force(lo(1)-1,j,k,2:nscal) = ext_scal_force(lo(1)-1,j,k,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          laps_local(2:nscal) = diff_coef * diff_fac * laps(hi(1),j,k,2:nscal)
          scal_force(hi(1)+1,j,k,2:nscal) = ext_scal_force(hi(1)+1,j,k,2:nscal) + laps_local(2:nscal)
       enddo
    enddo

  end subroutine mkscalforce_3d

end module mkforce_module
