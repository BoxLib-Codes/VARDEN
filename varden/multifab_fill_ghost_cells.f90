module multifab_fill_ghost_module

  use layout_module
  use fab_module
  use bl_mem_stat_module
  use multifab_module
  use bc_module
  use setbc_module
  use interp_module

  implicit none

contains

  subroutine multifab_fill_ghost_cells(fine,crse,fine_domain,ng,ir,bc,icomp,bc_comp,nc)

    type(multifab), intent(inout) :: fine
    type(multifab), intent(inout) :: crse
    type(box     ), intent(in   ) :: fine_domain
    integer       , intent(in   ) :: ng
    integer       , intent(in   ) :: ir(:)
    integer       , intent(in   ) :: bc(:,:,:)
    integer       , intent(in   ) :: icomp,bc_comp,nc

    type(box) :: fbox,fstrip,cstrip
    type(box) :: crse_domain
    integer   :: i,j,n,d,dm
    integer   :: iface,ishift
    integer   :: lo(4),hi(4),interior_lo(4)
    logical   :: lim_slope, lin_limit

    real(kind=dp_t), allocatable :: cp(:,:,:,:)
    real(kind=dp_t), pointer     :: fp(:,:,:,:)
    real(kind=dp_t), pointer     :: fine_fab(:,:,:,:)
    real(kind=dp_t), allocatable :: fvcx(:),fvcy(:),fvcz(:)
    real(kind=dp_t), allocatable :: cvcx(:),cvcy(:),cvcz(:)
    integer                      :: local_bc(fine%dim,2,nc)

    real(kind=dp_t) :: dx(2)
    dx = ONE

    dm = fine%dim

    lo = 1
    hi = 1

    hi(4) = nc

    lim_slope = .true.
    lin_limit = .false.

    local_bc = INTERIOR

    crse_domain = box_coarsen_v(fine_domain,ir)

    do i = 1, fine%nboxes
      fp => dataptr(fine,i)

        do iface = -1,1,2
          do d = 1,dm
             fbox   = get_ibox(fine,i)
             fstrip  = grow(fbox,ng)
             ishift = iface * (ng + box_extent_d(fbox,d))
             fstrip  = box_intersection(fstrip,box_shift_d(fstrip,ishift,d))
             fstrip  = box_intersection(fstrip,fine_domain)

             if (.not. box_empty(fstrip)) then

                cstrip = box_coarsen_v(fstrip,ir)
                interior_lo(1:dm) = lwb(cstrip)

                if (cstrip%lo(1) == crse_domain%lo(1)) then
                   local_bc(1,1,1:nc) = bc(1,1,bc_comp:bc_comp+nc-1)
                end if
                if (cstrip%hi(1) == crse_domain%hi(1)) then
                   local_bc(1,2,1:nc) = bc(1,2,bc_comp:bc_comp+nc-1)
                end if
                if (dm > 1) then
                   if (cstrip%lo(2) == crse_domain%lo(2)) then
                      local_bc(2,1,1:nc) = bc(2,1,bc_comp:bc_comp+nc-1)
                   end if
                   if (cstrip%hi(2) == crse_domain%hi(2)) then
                      local_bc(2,2,1:nc) = bc(2,2,bc_comp:bc_comp+nc-1)
                   end if
                end if
                if (dm > 2) then
                   if (cstrip%lo(3) == crse_domain%lo(3)) then
                      local_bc(3,1,1:nc) = bc(3,1,bc_comp:bc_comp+nc-1)
                   end if
                   if (cstrip%hi(3) == crse_domain%hi(3)) then
                      local_bc(3,2,1:nc) = bc(3,2,bc_comp:bc_comp+nc-1)
                   end if
                end if

                lo(1:dm) = lwb(fstrip)
                hi(1:dm) = upb(fstrip)
                allocate(fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4)))
                allocate(fvcx(lo(1):hi(1)+1))
                do j = lo(1),hi(1)+1
                   fvcx(j) = dble(j)
                end do

                if (dm > 1) then
                   allocate(fvcy(lo(2):hi(2)+1))
                   do j = lo(2),hi(2)+1
                      fvcy(j) = dble(j) 
                   end do
                end if
                if (dm > 2) then
                   allocate(fvcz(lo(3):hi(3)+1))
                   do j = lo(3),hi(3)+1
                      fvcy(j) = dble(j) 
                   end do
                end if

                cstrip = grow(cstrip,1)

                lo(1:dm) = lwb(cstrip)
                hi(1:dm) = upb(cstrip)
                allocate(cp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),lo(4):hi(4)))
                allocate(cvcx(lo(1):hi(1)+1))
                do j = lo(1),hi(1)+1
                   cvcx(j) = dble(j) * TWO
                end do
                if (dm > 1) then
                   allocate(cvcy(lo(2):hi(2)+1))
                   do j = lo(2),hi(2)+1
                      cvcy(j) = dble(j) * TWO
                   end do
                end if
                if (dm > 2) then
                   allocate(cvcz(lo(3):hi(3)+1))
                   do j = lo(3),hi(3)+1
                      cvcy(j) = dble(j) * TWO
                   end do
                end if
 
                call multifab_fab_copy(cp,lo,1,crse,icomp,nc)

!               Note: these calls assume that the crse array is one larger than
!                     the coarsened fine array in the lo- AND hi-directions in each dimension.
                select case (dm)
                case (2)
                   do n = 1,nc
                      call setbc_2d(cp(:,:,1,n),interior_lo,1,local_bc(:,:,n),dx,bc_comp+n-1)
                   end do
                   call lin_cc_interp_2d(fp(:,:,1,:), cp(:,:,1,:), ir, local_bc, &
                                         fvcx, fvcy, cvcx, cvcy, &
                                         lim_slope, lin_limit)
                case (3)
                   do n = 1,nc
                      call setbc_3d(cp(:,:,:,n),interior_lo,1,local_bc(:,:,n),dx,bc_comp+n-1)
                   end do
                   call lin_cc_interp_3d(fp(:,:,:,:), cp(:,:,:,:), ir, local_bc, &
                                         fvcx, fvcy, fvcz, cvcx, cvcy, cvcz, &
                                         lim_slope, lin_limit)
                end select
   
                deallocate(cp)
                deallocate(cvcx)
                deallocate(fvcx)
                if (dm > 1) deallocate(cvcy)
                if (dm > 1) deallocate(fvcy)
                if (dm > 2) deallocate(cvcz)
                if (dm > 2) deallocate(fvcz)

                if (dm < 3) fstrip%lo(3) = 1
                if (dm < 3) fstrip%hi(3) = 1

                fine_fab => dataptr(fine,i)
                fine_fab(fstrip%lo(1):fstrip%hi(1),fstrip%lo(2):fstrip%hi(2), &
                         fstrip%lo(3):fstrip%hi(3),1:nc) = &
                      fp(fstrip%lo(1):fstrip%hi(1),fstrip%lo(2):fstrip%hi(2), &
                         fstrip%lo(3):fstrip%hi(3),1:nc)

             end if
         end do
       end do
     end do

  end subroutine multifab_fill_ghost_cells

end module multifab_fill_ghost_module
