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
    
    type(bl_prof_timer), save :: bpt
    
    call build(bpt,"multifab_physbc")

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp)
          end do
       end select
    end do
 
    call destroy(bpt)

  end subroutine multifab_physbc

  subroutine physbc_2d(s,lo,ng,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j,hi(2)
    integer :: ngylo, ngyhi

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

    if (bc(2,1) .eq. INTERIOR) then  ! y-lo direction
       ngylo = ng
    else
       ngylo = 0  ! avoid corner ghost cells because we don't have good data in y-lo faces
    end if

    if (bc(2,2) .eq. INTERIOR) then  ! y-hi direction
       ngyhi = ng
    else
       ngyhi = 0
    end if

    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = u_bc(1,1)
       if (icomp.eq.2) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = v_bc(1,1)
       if (icomp.eq.3) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = rho_bc(1,1)
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = trac_bc(1,1)
    else if (bc(1,1) .eq. FOEXTRAP) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          s(lo(1)-ng:lo(1)-1,j) = EIGHTH* &
               (15.0_dp_t*s(lo(1),j) - 10.0_dp_t*s(lo(1)+1,j) + 3.0_dp_t*s(lo(1)+2,j))
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          do i = 1, ng
             s(lo(1)-i,j) = s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          do i = 1, ng
             s(lo(1)-i,j) = -s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = u_bc(1,2)
       if (icomp.eq.2) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = v_bc(1,2)
       if (icomp.eq.3) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = rho_bc(1,2)
       if (icomp.eq.4) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = trac_bc(1,2)
    else if (bc(1,2) .eq. FOEXTRAP) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          s(hi(1)+1:hi(1)+ng,j) = EIGHTH* &
               (15.0_dp_t * s(hi(1)  ,j) - 10.0_dp_t*s(hi(1)-1,j) + 3.0_dp_t*s(hi(1)-2,j))
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          do i = 1, ng
             s(hi(1)+i,j) = s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do j = lo(2)-ngylo, hi(2)+ngyhi
          do i = 1, ng
             s(hi(1)+i,j) = -s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = u_bc(2,1)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = v_bc(2,1)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = rho_bc(2,1)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = trac_bc(2,1)
    else if (bc(2,1) .eq. FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2))
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = EIGHTH* &
               (15.0_dp_t*s(i,lo(2)) - 10.0_dp_t*s(i,lo(2)+1) + 3.0_dp_t*s(i,lo(2)+2))
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do j = 1, ng
          do i = lo(1)-ng, hi(1)+ng
             s(i,lo(2)-j) = s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do j = 1, ng
          do i = lo(1)-ng, hi(1)+ng
             s(i,lo(2)-j) = -s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = u_bc(2,2)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = v_bc(2,2)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = rho_bc(2,2)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = trac_bc(2,2)
    else if (bc(2,2) .eq. FOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2))
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do i = lo(1)-ng, hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = EIGHTH* &
               (15.0_dp_t*s(i,hi(2)) - 10.0_dp_t*s(i,hi(2)-1) + 3.0_dp_t*s(i,hi(2)-2))
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do j = 1, ng
          do i = lo(1)-ng, hi(1)+ng
             s(i,hi(2)+j) = s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do j = 1, ng
          do i = lo(1)-ng, hi(1)+ng
             s(i,hi(2)+j) = -s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,ng,bc,icomp)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp

    ! Local variables
    integer :: i,j,k,hi(3)
    integer :: ngylo, ngyhi, ngzlo, ngzhi

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    if (bc(2,1) .eq. INTERIOR) then  ! y-lo direction
       ngylo = ng
    else
       ngylo = 0  ! avoid corner ghost cells because we don't have good data in y-lo faces
    end if

    if (bc(2,2) .eq. INTERIOR) then  ! y-hi direction
       ngyhi = ng
    else
       ngyhi = 0
    end if

    if (bc(3,1) .eq. INTERIOR) then  ! z-lo direction
       ngzlo = ng
    else
       ngzlo = 0
    end if

    if (bc(3,2) .eq. INTERIOR) then  ! z-hi direction
       ngzhi = ng
    else
       ngzhi = 0
    end if

    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = u_bc(1,1)
       if (icomp.eq.2) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = v_bc(1,1)
       if (icomp.eq.3) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = w_bc(1,1)
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = rho_bc(1,1)
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = trac_bc(1,1)
    else if (bc(1,1) .eq. FOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k)
          end do
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             s(lo(1)-ng:lo(1)-1,j,k) = &
                  ( 15.0_dp_t * s(lo(1)  ,j,k) &
                  -10.0_dp_t * s(lo(1)+1,j,k) &
                  + 3.0_dp_t * s(lo(1)+2,j,k) ) * EIGHTH
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             do i = 1,ng
                s(lo(1)-i,j,k) = s(lo(1)+i-1,j,k)
             end do
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             do i = 1,ng
                s(lo(1)-i,j,k) = -s(lo(1)+i-1,j,k)
             end do
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,1) = ',bc(1,1)
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = u_bc(1,2)
       if (icomp.eq.2) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = v_bc(1,2)
       if (icomp.eq.3) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = w_bc(1,2)
       if (icomp.eq.4) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = rho_bc(1,2)
       if (icomp.eq.5) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = trac_bc(1,2)
    else if (bc(1,2) .eq. FOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             s(hi(1)+1:hi(1)+ng,j,k) = &
                  ( 15.0_dp_t * s(hi(1)  ,j,k) &
                  -10.0_dp_t * s(hi(1)-1,j,k) &
                  + 3.0_dp_t * s(hi(1)-2,j,k) ) * EIGHTH
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             do i = 1,ng
                s(hi(1)+i,j,k) = s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = lo(2)-ngylo,hi(2)+ngyhi
             do i = 1,ng
                s(hi(1)+i,j,k) = -s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(1,2) = ',bc(1,2)
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = u_bc(2,1)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = v_bc(2,1)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = w_bc(2,1)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = rho_bc(2,1)
       if (icomp.eq.5) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = trac_bc(2,1)
    else if (bc(2,1) .eq. FOEXTRAP) then 
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = &
                  ( 15.0_dp_t * s(i,lo(2)  ,k) &
                  -10.0_dp_t * s(i,lo(2)+1,k) &
                  + 3.0_dp_t * s(i,lo(2)+2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = 1,ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,lo(2)-j,k) = s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = 1,ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,lo(2)-j,k) = -s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,1) = ',bc(2,1)
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = u_bc(2,2)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = v_bc(2,2)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = w_bc(2,2)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = rho_bc(2,2)
       if (icomp.eq.5) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = trac_bc(2,2)
    else if (bc(2,2) .eq. FOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = &
                  ( 15.0_dp_t * s(i,hi(2)  ,k) &
                  -10.0_dp_t * s(i,hi(2)-1,k) &
                  + 3.0_dp_t * s(i,hi(2)-2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = 1,ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,hi(2)+j,k) = s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ngzlo,hi(3)+ngzhi
          do j = 1,ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,hi(2)+j,k) = -s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(2,2) = ',bc(2,2)
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

    if (bc(3,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = u_bc(3,1)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = v_bc(3,1)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = w_bc(3,1)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = rho_bc(3,1)
       if (icomp.eq.5) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = trac_bc(3,1)
    else if (bc(3,1) .eq. FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3))
          end do
       end do
    else if (bc(3,1) .eq. HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = &
                  ( 15.0_dp_t * s(i,j,lo(3)  ) &
                  -10.0_dp_t * s(i,j,lo(3)+1) &
                  + 3.0_dp_t * s(i,j,lo(3)+2) ) * EIGHTH
          end do
       end do
    else if (bc(3,1) .eq. REFLECT_EVEN) then
       do k = 1,ng
          do j = lo(2)-ng,hi(2)+ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,j,lo(3)-k) = s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. REFLECT_ODD) then
       do k = 1,ng
          do j = lo(2)-ng,hi(2)+ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,j,lo(3)-k) = -s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(3,1) = ',bc(3,1)
       call bl_error('BC(3,1) = NOT YET SUPPORTED')
    end if

    if (bc(3,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = u_bc(3,2)
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = v_bc(3,2)
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = w_bc(3,2)
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = rho_bc(3,2)
       if (icomp.eq.5) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = trac_bc(3,2)
    else if (bc(3,2) .eq. FOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3))
          end do
       end do
    else if (bc(3,2) .eq. HOEXTRAP) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = &
                  ( 15.0_dp_t * s(i,j,hi(3)  ) &
                  -10.0_dp_t * s(i,j,hi(3)-1) &
                  + 3.0_dp_t * s(i,j,hi(3)-2) ) * EIGHTH
          end do
       end do
    else if (bc(3,2) .eq. REFLECT_EVEN) then
       do k = 1,ng
          do j = lo(2)-ng,hi(2)+ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,j,hi(3)+k) = s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. REFLECT_ODD) then
       do k = 1,ng
          do j = lo(2)-ng,hi(2)+ng
             do i = lo(1)-ng,hi(1)+ng
                s(i,j,hi(3)+k) = -s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! do nothing
    else 
       print *,'bc(3,2) = ',bc(3,2)
       call bl_error('BC(3,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
