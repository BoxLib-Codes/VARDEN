module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bl_types
  use bl_error_module
  use probin_module, only: rho_bc, trac_bc, u_bc, v_bc, w_bc

  implicit none

  private

  public :: multifab_physbc

contains

  subroutine multifab_physbc(s,start_scomp,start_bccomp,num_comp,the_bc_level, &
                             time_in,dx_in,prob_lo_in,prob_hi_in)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ), optional :: time_in,dx_in(:),prob_lo_in(:),prob_hi_in(:)

    ! Local
    integer                  :: lo(get_dim(s))
    integer                  :: i,ng,dm,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    
    type(mfiter) :: mfi
    type(box) :: growntilebox
    integer :: gtlo(get_dim(s)), gthi(get_dim(s))

    type(bl_prof_timer), save :: bpt
    
    call build(bpt,"multifab_physbc")

    ng = nghost(s)
    dm = get_dim(s)
    
    !$omp parallel private(mfi,i,growntilebox,gtlo,gthi) &
    !$omp private(sp,lo,scomp,bccomp)

    call mfiter_build(mfi,s,tiling=.true.)

    do while(more_tile(mfi))
       i = get_fab_index(mfi)

       growntilebox = get_growntilebox(mfi,ng)
       gtlo = lwb(growntilebox)
       gthi = upb(growntilebox)

       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, &
                            gtlo,gthi)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, &
                            gtlo,gthi)
          end do
       end select
    end do
    !$omp end parallel
 
    call destroy(bpt)

  end subroutine multifab_physbc

  subroutine physbc_2d(s,glo,ng,bc,icomp,gtlo,gthi)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: glo(:),ng,gtlo(:),gthi(:)
    real(kind=dp_t), intent(inout) :: s(glo(1)-ng:,glo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j,ghi(2)

    ghi(1) = glo(1) + size(s,dim=1) - (2*ng+1)
    ghi(2) = glo(2) + size(s,dim=2) - (2*ng+1)

    if (gtlo(1)+ng .eq. glo(1)) then
       if (bc(1,1) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2)) = u_bc(1,1)
          if (icomp.eq.2) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2)) = v_bc(1,1)
          if (icomp.eq.3) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2)) = rho_bc(1,1)
          if (icomp.eq.4) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2)) = trac_bc(1,1)
       else if (bc(1,1) .eq. FOEXTRAP) then
          do j = gtlo(2), gthi(2)
             s(gtlo(1):gtlo(1)+ng-1,j) = s(gtlo(1)+ng,j)
          end do
       else if (bc(1,1) .eq. HOEXTRAP) then
          do j = gtlo(2), gthi(2)
             s(gtlo(1):gtlo(1)+ng-1,j) = EIGHTH* &
                  ( 15.0_dp_t*s(gtlo(1)+ng  ,j) &
                   -10.0_dp_t*s(gtlo(1)+ng+1,j) &
                   + 3.0_dp_t*s(gtlo(1)+ng+2,j) )
          end do
       else if (bc(1,1) .eq. REFLECT_EVEN) then
          do j = gtlo(2), gthi(2)
             do i = 1, ng
                s(gtlo(1)+ng-i,j) = s(gtlo(1)+ng+i-1,j)
             end do
          end do
       else if (bc(1,1) .eq. REFLECT_ODD) then
          do j = gtlo(2), gthi(2)
             do i = 1, ng
                s(gtlo(1)+ng-i,j) = -s(gtlo(1)+ng+i-1,j)
             end do
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(1,1) = ',bc(1,1)
          call bl_error('BC(1,1) = NOT YET SUPPORTED')
       end if
    end if

    if (gthi(1)-ng .eq. ghi(1)) then
       if (bc(1,2) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2)) = u_bc(1,2)
          if (icomp.eq.2) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2)) = v_bc(1,2)
          if (icomp.eq.3) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2)) = rho_bc(1,2)
          if (icomp.eq.4) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2)) = trac_bc(1,2)
       else if (bc(1,2) .eq. FOEXTRAP) then
          do j = gtlo(2), gthi(2)
             s(gthi(1)-ng+1:gthi(1),j) = s(gthi(1)-ng,j)
          end do
       else if (bc(1,2) .eq. HOEXTRAP) then
          do j = gtlo(2), gthi(2)
             s(gthi(1)-ng+1:gthi(1),j) = EIGHTH* &
                  ( 15.0_dp_t * s(gthi(1)-ng  ,j) &
                   -10.0_dp_t * s(gthi(1)-ng-1,j) &
                   + 3.0_dp_t * s(gthi(1)-ng-2,j) )
          end do
       else if (bc(1,2) .eq. REFLECT_EVEN) then
          do j = gtlo(2), gthi(2)
             do i = 1, ng
                s(gthi(1)-ng+i,j) = s(gthi(1)-ng-i+1,j)
             end do
          end do
       else if (bc(1,2) .eq. REFLECT_ODD) then
          do j = gtlo(2), gthi(2)
             do i = 1, ng
                s(gthi(1)-ng+i,j) = -s(gthi(1)-ng-i+1,j)
             end do
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(1,2) = ',bc(1,2)
          call bl_error('BC(1,2) = NOT YET SUPPORTED')
       end if
    end if

    if (gtlo(2)+ng .eq. glo(2)) then
       if (bc(2,1) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1) = u_bc(2,1)
          if (icomp.eq.2) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1) = v_bc(2,1)
          if (icomp.eq.3) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1) = rho_bc(2,1)
          if (icomp.eq.4) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1) = trac_bc(2,1)
       else if (bc(2,1) .eq. FOEXTRAP) then
          do i = gtlo(1), gthi(1)
             s(i,gtlo(2):gtlo(2)+ng-1) = s(i,gtlo(2)+ng)
          end do
       else if (bc(2,1) .eq. HOEXTRAP) then
          do i = gtlo(1), gthi(1)
             s(i,gtlo(2):gtlo(2)+ng-1) = EIGHTH* &
                 ( 15.0_dp_t*s(i,gtlo(2)+ng) &
                  -10.0_dp_t*s(i,gtlo(2)+ng+1) &
                  + 3.0_dp_t*s(i,gtlo(2)+ng+2) )
          end do
       else if (bc(2,1) .eq. REFLECT_EVEN) then
          do i = gtlo(1), gthi(1)
             do j = 1, ng
                s(i,gtlo(2)+ng-j) = s(i,gtlo(2)+ng+j-1)
             end do
          end do
       else if (bc(2,1) .eq. REFLECT_ODD) then
          do i = gtlo(1), gthi(1)
             do j = 1, ng
                s(i,gtlo(2)+ng-j) = -s(i,gtlo(2)+ng+j-1)
             end do
          end do
       else if (bc(2,1) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(2,1) = ',bc(2,1)
          call bl_error('BC(2,1) = NOT YET SUPPORTED')
       end if
    end if

    if (gthi(2)-ng .eq. ghi(2)) then
       if (bc(2,2) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2)) = u_bc(2,2)
          if (icomp.eq.2) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2)) = v_bc(2,2)
          if (icomp.eq.3) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2)) = rho_bc(2,2)
          if (icomp.eq.4) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2)) = trac_bc(2,2)
       else if (bc(2,2) .eq. FOEXTRAP) then
          do i = gtlo(1), gthi(1)
             s(i,gthi(2)-ng+1:gthi(2)) = s(i,gthi(2)-ng)
          end do
       else if (bc(2,2) .eq. HOEXTRAP) then
          do i = gtlo(1), gthi(1)
             s(i,gthi(2)-ng+1:gthi(2)) = EIGHTH* &
                 ( 15.0_dp_t*s(i,gthi(2)-ng) &
                  -10.0_dp_t*s(i,gthi(2)-ng-1) &
                  + 3.0_dp_t*s(i,gthi(2)-ng-2) )
          end do
       else if (bc(2,2) .eq. REFLECT_EVEN) then
          do i = gtlo(1), gthi(1)
             do j = 1, ng
                s(i,gthi(2)-ng+j) = s(i,gthi(2)-ng-j+1)
             end do
          end do
       else if (bc(2,2) .eq. REFLECT_ODD) then
          do i = gtlo(1), gthi(1)
             do j = 1, ng
                s(i,gthi(2)-ng+j) = -s(i,gthi(2)-ng-j+1)
             end do
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(2,2) = ',bc(2,2)
          call bl_error('BC(2,2) = NOT YET SUPPORTED')
       end if
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,glo,ng,bc,icomp,gtlo,gthi)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: glo(:),ng,gtlo(:),gthi(:)
    real(kind=dp_t), intent(inout) :: s(glo(1)-ng:, glo(2)-ng:, glo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j,k,ghi(3)

    ghi(1) = glo(1) + size(s,dim=1) - (2*ng+1)
    ghi(2) = glo(2) + size(s,dim=2) - (2*ng+1)
    ghi(3) = glo(3) + size(s,dim=3) - (2*ng+1)

    if (gtlo(1)+ng .eq. glo(1)) then
       if (bc(1,1) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2),gtlo(3):gthi(3)) = u_bc(1,1)
          if (icomp.eq.2) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2),gtlo(3):gthi(3)) = v_bc(1,1)
          if (icomp.eq.3) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2),gtlo(3):gthi(3)) = w_bc(1,1)
          if (icomp.eq.4) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2),gtlo(3):gthi(3)) = rho_bc(1,1)
          if (icomp.eq.5) s(gtlo(1):gtlo(1)+ng-1,gtlo(2):gthi(2),gtlo(3):gthi(3)) = trac_bc(1,1)
       else if (bc(1,1) .eq. FOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                s(gtlo(1):gtlo(1)+ng-1,j,k) = s(gtlo(1)+ng,j,k)
             end do
          end do
       else if (bc(1,1) .eq. HOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                s(gtlo(1):gtlo(1)+ng-1,j,k) = &
                     ( 15.0_dp_t * s(gtlo(1)+ng  ,j,k) &
                     -10.0_dp_t * s(gtlo(1)+ng+1,j,k) &
                     + 3.0_dp_t * s(gtlo(1)+ng+2,j,k) ) * EIGHTH
             end do
          end do
       else if (bc(1,1) .eq. REFLECT_EVEN) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                do i = 1,ng
                   s(gtlo(1)+ng-i,j,k) = s(gtlo(1)+ng+i-1,j,k)
                end do
             end do
          end do
       else if (bc(1,1) .eq. REFLECT_ODD) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                do i = 1,ng
                   s(gtlo(1)+ng-i,j,k) = -s(gtlo(1)+ng+i-1,j,k)
                end do
             end do
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(1,1) = ',bc(1,1)
          call bl_error('BC(1,1) = NOT YET SUPPORTED')
       end if
    end if

    if (gthi(1)-ng .eq. ghi(1)) then
       if (bc(1,2) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2),gtlo(3):gthi(3)) = u_bc(1,2)
          if (icomp.eq.2) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2),gtlo(3):gthi(3)) = v_bc(1,2)
          if (icomp.eq.3) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2),gtlo(3):gthi(3)) = w_bc(1,2)
          if (icomp.eq.4) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2),gtlo(3):gthi(3)) = rho_bc(1,2)
          if (icomp.eq.5) s(gthi(1)-ng+1:gthi(1),gtlo(2):gthi(2),gtlo(3):gthi(3)) = trac_bc(1,2)
       else if (bc(1,2) .eq. FOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                s(gthi(1)-ng+1:gthi(1),j,k) = s(gthi(1)-ng,j,k)
             end do
          end do
       else if (bc(1,2) .eq. HOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                s(gthi(1)-ng+1:gthi(1),j,k) = &
                     ( 15.0_dp_t * s(gthi(1)-ng  ,j,k) &
                     -10.0_dp_t * s(gthi(1)-ng-1,j,k) &
                     + 3.0_dp_t * s(gthi(1)-ng-2,j,k) ) * EIGHTH
             end do
          end do
       else if (bc(1,2) .eq. REFLECT_EVEN) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                do i = 1,ng
                   s(gthi(1)-ng+i,j,k) = s(gthi(1)-ng-i+1,j,k)
                end do
             end do
          end do
       else if (bc(1,2) .eq. REFLECT_ODD) then
          do k = gtlo(3),gthi(3)
             do j = gtlo(2),gthi(2)
                do i = 1,ng
                   s(gthi(1)-ng+i,j,k) = -s(gthi(1)-ng-i+1,j,k)
                end do
             end do
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(1,2) = ',bc(1,2)
          call bl_error('BC(1,2) = NOT YET SUPPORTED')
       end if
    end if

   if (gtlo(2)+ng .eq. glo(2)) then
      if (bc(2,1) .eq. EXT_DIR) then
         if (icomp.eq.1) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1,gtlo(3):gthi(3)) = u_bc(2,1)
         if (icomp.eq.2) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1,gtlo(3):gthi(3)) = v_bc(2,1)
         if (icomp.eq.3) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1,gtlo(3):gthi(3)) = w_bc(2,1)
         if (icomp.eq.4) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1,gtlo(3):gthi(3)) = rho_bc(2,1)
         if (icomp.eq.5) s(gtlo(1):gthi(1),gtlo(2):gtlo(2)+ng-1,gtlo(3):gthi(3)) = trac_bc(2,1)
      else if (bc(2,1) .eq. FOEXTRAP) then
         do k = gtlo(3),gthi(3)
            do i = gtlo(1),gthi(1)
               s(i,gtlo(2):gtlo(2)+ng-1,k) = s(i,gtlo(2)+ng,k)
            end do
         end do
      else if (bc(2,1) .eq. HOEXTRAP) then
         do k = gtlo(3),gthi(3)
            do i = gtlo(1),gthi(1)
               s(i,gtlo(2):gtlo(2)+ng-1,k) = &
                    ( 15.0_dp_t * s(i,gtlo(2)+ng  ,k) &
                    -10.0_dp_t * s(i,gtlo(2)+ng+1,k) &
                    + 3.0_dp_t * s(i,gtlo(2)+ng+2,k) ) * EIGHTH
            end do
         end do
      else if (bc(2,1) .eq. REFLECT_EVEN) then
         do k = gtlo(3),gthi(3)
            do i = gtlo(1),gthi(1)
               do j = 1,ng
                  s(i,gtlo(2)+ng-j,k) = s(i,gtlo(2)+ng+j-1,k)
               end do
            end do
         end do
      else if (bc(2,1) .eq. REFLECT_ODD) then
         do k = gtlo(3),gthi(3)
            do i = gtlo(1),gthi(1)
               do j = 1,ng
                  s(i,gtlo(2)+ng-j,k) = -s(i,gtlo(2)+ng+j-1,k)
               end do
            end do
         end do
      else if (bc(2,1) .eq. INTERIOR) then
         ! do nothing
      else 
         print *,'bc(2,1) = ',bc(2,1)
         call bl_error('BC(2,1) = NOT YET SUPPORTED')
      end if
    end if

    if (gthi(2)-ng .eq. ghi(2)) then
       if (bc(2,2) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2),gtlo(3):gthi(3)) = u_bc(2,2)
          if (icomp.eq.2) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2),gtlo(3):gthi(3)) = v_bc(2,2)
          if (icomp.eq.3) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2),gtlo(3):gthi(3)) = w_bc(2,2)
          if (icomp.eq.4) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2),gtlo(3):gthi(3)) = rho_bc(2,2)
          if (icomp.eq.5) s(gtlo(1):gthi(1),gthi(2)-ng+1:gthi(2),gtlo(3):gthi(3)) = trac_bc(2,2)
       else if (bc(2,2) .eq. FOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do i = gtlo(1),gthi(1)
                s(i,gthi(2)-ng+1:gthi(2),k) = s(i,gthi(2)-ng,k)
             end do
          end do
       else if (bc(2,2) .eq. HOEXTRAP) then
          do k = gtlo(3),gthi(3)
             do i = gtlo(1),gthi(1)
                s(i,gthi(2)-ng+1:gthi(2),k) = &
                    ( 15.0_dp_t * s(i,gthi(2)-ng  ,k) &
                     -10.0_dp_t * s(i,gthi(2)-ng-1,k) &
                     + 3.0_dp_t * s(i,gthi(2)-ng-2,k) ) * EIGHTH
             end do
          end do
       else if (bc(2,2) .eq. REFLECT_EVEN) then
          do k = gtlo(3),gthi(3)
             do i = gtlo(1),gthi(1)
                do j = 1,ng
                   s(i,gthi(2)-ng+j,k) = s(i,gthi(2)-ng-j+1,k)
                end do
             end do
          end do
       else if (bc(2,2) .eq. REFLECT_ODD) then
          do k = gtlo(3),gthi(3)
             do i = gtlo(1),gthi(1)
                do j = 1,ng
                   s(i,gthi(2)-ng+j,k) = -s(i,gthi(2)-ng-j+1,k)
                end do
             end do
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(2,2) = ',bc(2,2)
          call bl_error('BC(2,2) = NOT YET SUPPORTED')
       end if
    end if

    if (gtlo(3)+ng .eq. glo(3)) then
       if (bc(3,1) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gtlo(3):gtlo(3)+ng-1) = u_bc(3,1)
          if (icomp.eq.2) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gtlo(3):gtlo(3)+ng-1) = v_bc(3,1)
          if (icomp.eq.3) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gtlo(3):gtlo(3)+ng-1) = w_bc(3,1)
          if (icomp.eq.4) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gtlo(3):gtlo(3)+ng-1) = rho_bc(3,1)
          if (icomp.eq.5) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gtlo(3):gtlo(3)+ng-1) = trac_bc(3,1)
       else if (bc(3,1) .eq. FOEXTRAP) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                s(i,j,gtlo(3):gtlo(3)+ng-1) = s(i,j,gtlo(3)+ng)
             end do
          end do
       else if (bc(3,1) .eq. HOEXTRAP) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                s(i,j,gtlo(3):gtlo(3)+ng-1) = &
                    ( 15.0_dp_t * s(i,j,gtlo(3)+ng  ) &
                     -10.0_dp_t * s(i,j,gtlo(3)+ng+1) &
                     + 3.0_dp_t * s(i,j,gtlo(3)+ng+2) ) * EIGHTH
             end do
          end do
       else if (bc(3,1) .eq. REFLECT_EVEN) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                do k = 1,ng
                   s(i,j,gtlo(3)+ng-k) = s(i,j,gtlo(3)+ng+k-1)
                end do
             end do
          end do
       else if (bc(3,1) .eq. REFLECT_ODD) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                do k = 1,ng
                   s(i,j,gtlo(3)+ng-k) = -s(i,j,gtlo(3)+ng+k-1)
                end do
             end do
          end do
       else if (bc(3,1) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(3,1) = ',bc(3,1)
          call bl_error('BC(3,1) = NOT YET SUPPORTED')
       end if
    end if

    if (gthi(3)-ng .eq. ghi(3)) then
       if (bc(3,2) .eq. EXT_DIR) then
          if (icomp.eq.1) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gthi(3)-ng+1:gthi(3)) = u_bc(3,2)
          if (icomp.eq.2) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gthi(3)-ng+1:gthi(3)) = v_bc(3,2)
          if (icomp.eq.3) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gthi(3)-ng+1:gthi(3)) = w_bc(3,2)
          if (icomp.eq.4) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gthi(3)-ng+1:gthi(3)) = rho_bc(3,2)
          if (icomp.eq.5) s(gtlo(1):gthi(1),gtlo(2):gthi(2),gthi(3)-ng+1:gthi(3)) = trac_bc(3,2)
       else if (bc(3,2) .eq. FOEXTRAP) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                s(i,j,gthi(3)-ng+1:gthi(3)) = s(i,j,gthi(3)-ng)
             end do
          end do
       else if (bc(3,2) .eq. HOEXTRAP) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                s(i,j,gthi(3)-ng+1:gthi(3)) = &
                    ( 15.0_dp_t * s(i,j,gthi(3)-ng  ) &
                     -10.0_dp_t * s(i,j,gthi(3)-ng-1) &
                     + 3.0_dp_t * s(i,j,gthi(3)-ng-2) ) * EIGHTH
             end do
          end do
       else if (bc(3,2) .eq. REFLECT_EVEN) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                do k = 1,ng
                   s(i,j,gthi(3)-ng+k) = s(i,j,gthi(3)-ng-k+1)
                end do
             end do
          end do
       else if (bc(3,2) .eq. REFLECT_ODD) then
          do j = gtlo(2),gthi(2)
             do i = gtlo(1),gthi(1)
                do k = 1,ng
                   s(i,j,gthi(3)-ng+k) = -s(i,j,gthi(3)-ng-k+1)
                end do
             end do
          end do
       else if (bc(3,2) .eq. INTERIOR) then
          ! do nothing
       else 
          print *,'bc(3,2) = ',bc(3,2)
          call bl_error('BC(3,2) = NOT YET SUPPORTED')
       end if
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
