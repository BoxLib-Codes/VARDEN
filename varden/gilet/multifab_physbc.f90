module multifab_physbc_module

  use multifab_module
  use define_bc_module
  use bl_types
  use bl_error_module
  use probin_module, only : prob_hi_y

  implicit none

  private

  public :: multifab_physbc

contains

  subroutine multifab_physbc(s,start_scomp,start_bccomp,num_comp,the_bc_level,dx,t)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: t

    ! Local
    integer                  :: lo(s%dim)
    integer                  :: i,ng,dm,scomp,bccomp
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    
    ng = s%ng
    dm = s%dim
    
    do i=1,s%nboxes
       if ( multifab_remote(s,i) ) cycle
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       select case (dm)
       case (2)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_2d(sp(:,:,1,scomp), lo, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp,dx,t)
          end do
       case (3)
          do scomp = start_scomp,start_scomp+num_comp-1
             bccomp = start_bccomp + scomp - start_scomp
             call physbc_3d(sp(:,:,:,scomp), lo, ng, &
                           the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp,dx,t)
          end do
       end select
    end do
 
  end subroutine multifab_physbc

  subroutine physbc_2d(s,lo,ng,bc,icomp,dx,t)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: t

    ! Local variables
    integer, parameter  :: INLET_VX = zero
    integer, parameter  :: INLET_VY = zero
    integer, parameter  :: INLET_DEN = one
    integer, parameter  :: INLET_TRA_A = zero
    integer, parameter  :: INLET_TRA_B = zero
    integer, parameter  :: INLET_TRA_C = zero
    integer, parameter  :: INLET_TRA_D = zero
 
    real (kind = dp_t), parameter  :: lambda = fourth
    real (kind = dp_t), parameter  :: Um = one
    real (kind = dp_t), parameter  :: delta = 2.d-2


    integer :: i,j,hi(2)


    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

!    *** 2-D ***
! last index represents
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : density
!   COMP = 4 - (4+nscal-1) : tracer_i
!   last COMP   : pressure

!    *** 3-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : z-velocity
!   COMP = 4  : density
!   COMP = 5 - (5+nscal-1) : tracer_i
!   last COMP   : pressure

! bdry conditions for x at low end
    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp.eq.1) then   
          do j = lo(2)-ng,hi(2)+ng
             if (j*dx(2) < half*prob_hi_y) then
                s(lo(1)-ng:lo(1)-1,j) = one
             else
                s(lo(1)-ng:lo(1)-1,j) = 10.d0
             end if
             !s(lo(1)-ng:lo(1)-1,j) = Um*(one + lambda * tanh(two*(dx(2)*j-half)/delta))&
              !                       * (one + 1.d-2 * sin(306.6*t) + &
               !                               5.5d-3 * sin(102.2*t))
          enddo
       endif
       if (icomp.eq.2) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = zero
       if (icomp.eq.3) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = one
       if (icomp.eq.4) then
          do j = lo(2)-ng,hi(2)+ng
             if (j*dx(2) <= half*prob_hi_y) then
                s(lo(1)-ng:lo(1)-1,j) = one
             else
                s(lo(1)-ng:lo(1)-1,j) = zero
             end if
          enddo
       endif
       if (icomp.eq.5) then
          do j = lo(2)-ng,hi(2)+ng
             if (j*dx(2) <= half*prob_hi_y) then
                s(lo(1)-ng:lo(1)-1,j) = zero
             else
                s(lo(1)-ng:lo(1)-1,j) = one
             end if
          enddo
       endif
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = zero
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = zero
   else if (bc(1,1) .eq. FOEXTRAP) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2) 
          s(lo(1)-ng:lo(1)-1,j) = EIGHTH* &
               (15.0_dp_t*s(lo(1),j) - 10.0_dp_t*s(lo(1)+1,j) + 3.0_dp_t*s(lo(1)+2,j))
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          do i = 1, ng
             s(lo(1)-i,j) = s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          do i = 1, ng
             s(lo(1)-i,j) = -s(lo(1)+i-1,j)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if

! bdry conditions for x at high end
    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = INLET_VX
       if (icomp.eq.2) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = INLET_VY
       if (icomp.eq.3) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng) = INLET_DEN
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_A
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_B
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_C
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_D
    else if (bc(1,2) .eq. FOEXTRAP) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          s(hi(1)+1:hi(1)+ng,j) = EIGHTH* &
               (15.0_dp_t * s(hi(1)  ,j) - 10.0_dp_t*s(hi(1)-1,j) + 3.0_dp_t*s(hi(1)-2,j))
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          do i = 1, ng
             s(hi(1)+i,j) = s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
!       do j = lo(2)-ng, hi(2)+ng
       do j = lo(2), hi(2)
          do i = 1, ng
             s(hi(1)+i,j) = -s(hi(1)-i+1,j)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

! bdry conditions for y at lo end
    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1) = INLET_DEN
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_A
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_B
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_C
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_D
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
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,lo(2)-j) = -s(i,lo(2)+j-1)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

! bdry conditions for y at hi end
    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng) = INLET_DEN
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_A
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_B
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_C
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng) = INLET_TRA_D
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
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do i = lo(1)-ng, hi(1)+ng
          do j = 1, ng
             s(i,hi(2)+j) = -s(i,hi(2)-j+1)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,ng,bc,icomp,dx,t)

    use bl_constants_module
    use bc_module

    integer        , intent(in   ) :: lo(:),ng
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng:, lo(2)-ng:, lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: icomp
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: t

    ! Local variables
    integer, parameter  :: INLET_VX = zero
    integer, parameter  :: INLET_VY = zero
    integer, parameter  :: INLET_VZ = zero
    integer, parameter  :: INLET_DEN = one
    integer, parameter  :: INLET_TRA_A = zero
    integer, parameter  :: INLET_TRA_B = zero
    integer, parameter  :: INLET_TRA_C = zero
    integer, parameter  :: INLET_TRA_D = zero

    integer :: i,j,k,hi(3)

    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

    if (bc(1,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VZ
       if (icomp.eq.4) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(1,1) .eq. FOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k)
          end do
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = &
                  ( 15.0_dp_t * s(lo(1)  ,j,k) &
                  -10.0_dp_t * s(lo(1)+1,j,k) &
                  + 3.0_dp_t * s(lo(1)+2,j,k) ) * EIGHTH
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(lo(1)-i,j,k) = s(lo(1)+i-1,j,k)
             end do
          end do
       end do
    else if (bc(1,1) .eq. REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(lo(1)-i,j,k) = -s(lo(1)+1-i,j,k)
             end do
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(1,1) = NOT YET SUPPORTED')
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VX
       if (icomp.eq.2) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VY
       if (icomp.eq.3) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VZ
       if (icomp.eq.4) s(hi(1)+1:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(1,2) .eq. FOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k)
          end do
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = &
                  ( 15.0_dp_t * s(hi(1)  ,j,k) &
                  -10.0_dp_t * s(hi(1)-1,j,k) &
                  + 3.0_dp_t * s(hi(1)-2,j,k) ) * EIGHTH
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do j = lo(2)-ng,hi(2)+ng
             do i = 1,ng
                s(hi(1)+i,j,k) = -s(hi(1)-i+1,j,k)
             end do
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(1,2) = NOT YET SUPPORTED')
    end if

    if (bc(2,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = INLET_VZ
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:lo(2)-1,lo(3)-ng:hi(3)+ng) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k)
          end do
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = &
                  ( 15.0_dp_t * s(i,lo(2)  ,k) &
                  -10.0_dp_t * s(i,lo(2)+1,k) &
                  + 3.0_dp_t * s(i,lo(2)+2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. REFLECT_ODD) then
       do k = lo(3)-1,hi(3)+1
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,lo(2)-j,k) = -s(i,lo(2)+j-1,k)
             end do
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(2,1) = NOT YET SUPPORTED')
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_VZ
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,hi(2)+1:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k)
          end do
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = &
                  ( 15.0_dp_t * s(i,hi(2)  ,k) &
                  -10.0_dp_t * s(i,hi(2)-1,k) &
                  + 3.0_dp_t * s(i,hi(2)-2,k) ) * EIGHTH
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_EVEN) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. REFLECT_ODD) then
       do k = lo(3)-ng,hi(3)+ng
          do i = lo(1)-ng,hi(1)+ng
             do j = 1,ng
                s(i,hi(2)+j,k) = -s(i,hi(2)-j+1,k)
             end do
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(2,2) = NOT YET SUPPORTED')
    end if

    if (bc(3,1) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_VZ
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,lo(3)-ng:lo(3)-1) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. REFLECT_EVEN) then
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
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,lo(3)-k) = -s(i,j,lo(3)+k-1)
             end do
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(3,1) = NOT YET SUPPORTED')
    end if

    if (bc(3,2) .eq. EXT_DIR) then
       if (icomp.eq.1) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = INLET_VX
       if (icomp.eq.2) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = INLET_VY
       if (icomp.eq.3) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = INLET_VZ
       if (icomp.eq.4) s(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,hi(3)+1:hi(3)+ng) = INLET_DEN
       if (icomp.eq.5) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_A
       if (icomp.eq.6) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_B
       if (icomp.eq.7) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_C
       if (icomp.eq.8) s(lo(1)-ng:lo(1)-1,lo(2)-ng:hi(2)+ng,lo(3)-ng:hi(3)+ng) = INLET_TRA_D
    else if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. REFLECT_EVEN) then
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
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. REFLECT_ODD) then
       do j = lo(2)-ng,hi(2)+ng
          do i = lo(1)-ng,hi(1)+ng
             do k = 1,ng
                s(i,j,hi(3)+k) = -s(i,j,hi(3)-k+1)
             end do
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! do nothing
    else 
       call bl_error('BC(3,2) = NOT YET SUPPORTED')
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
