module sdc_interpolation

  use bl_types
  use multifab_module
  use bl_constants_module
  use probin_module, only : nscal, nspec

  implicit none
  
  ! number of terms at various times/quadrature points that
  ! need to be saved for advection and diffusion
  integer :: n_adv = 1
  integer :: n_diff = 4


  type dpt_pntr
     real(kind=dp_t), dimension(:), pointer :: p(:,:,:,:)
  end type dpt_pntr

  public ::  provisional_intgrl, intgrl, get_interp_pts_array

contains

  subroutine provisional_intgrl(I_ADR,sold,snew,adv,diff,dt,nlevs)
    !*******************************************************
    ! this subroutine computes:
    ! I/dt + Adv     
    ! where
    ! I/dt = {intgrl_rxns + intgrl(interpol of AD)}/dt 
    ! intgrl_rxns = s_new-s_old - dt*(provis approx to AD
    !
    ! used as a source in computing Adv in first iter of SDC
    !*******************************************************

    type(multifab) , intent(inout) :: I_ADR(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    type(multifab) , intent(in   ) :: adv(0:,:), diff(0:,:)
    real(dp_t),      intent(in   ) :: dt
    integer,         intent(in   ) :: nlevs

    ! local
    integer                 :: n,i,dm,ng,comp
    integer                 :: ix,iy,iz,A_last,D_last
    integer                 :: lo(sold(1)%dim),hi(sold(1)%dim)
    real(dp_t), pointer     :: soop(:,:,:,:),snop(:,:,:,:),Iop(:,:,:,:)
    type(dpt_pntr) :: A(0:size(adv,dim=2)-1), D(0:size(diff,dim=2)-1) 


    dm    = sold(1)%dim

    A_last = size(adv,dim=2)-1
    D_last = size(diff,dim=2)-1
    
    ng = sold(1)%ng

    do n=1,nlevs
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle      
          call get_interp_pts_array(adv,A,n,i,A_last)
          call get_interp_pts_array(diff,D,n,i,D_last)
          soop    => dataptr(sold(n), i)         
          snop    => dataptr(snew(n), i)
          Iop     => dataptr(I_ADR(n), i)
          lo = lwb(get_box(sold(n), i))
          hi = upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)

                   do comp = 2, nscal
                      Iop(ix,iy,1,comp-1) = (snop(ix,iy,1,comp) - soop(ix,iy,1,comp))/dt&
                           + half*(D(2)%p(ix,iy,1,comp-1) + D(0)%p(ix,iy,1,comp-1))&
                           - D(1)%p(ix,iy,1,comp-1) - A(0)%p(ix,iy,1,comp-1)
                   enddo
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 2, nscal
                         Iop(ix,iy,iz,comp-1) = (snop(ix,iy,iz,comp) - soop(ix,iy,iz,comp))/dt&
                              + half*(D(2)%p(ix,iy,iz,comp-1) + D(0)%p(ix,iy,iz,comp-1))&
                              - D(1)%p(ix,iy,iz,comp-1) - A(0)%p(ix,iy,iz,comp-1)
                      enddo

                   end do
                end do
             end do

          end select
       end do

    end do

  end subroutine provisional_intgrl


  subroutine intgrl(I_ADR,sold,snew,adv,diff,dt,nlevs)
    !*******************************************************
    ! this subroutine computes:
    ! I/dt + Adv
    ! where
    ! I/dt  = [intgrl_rxns + intgrl(interpol of AD)]/dt 
    !  
    ! used as a source for computing advection in SDC
    !*******************************************************

    type(multifab) , intent(inout) :: I_ADR(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    type(multifab) , intent(in   ) :: adv(0:,:), diff(0:,:)
    real(dp_t),      intent(in   ) :: dt
    integer,         intent(in   ) :: nlevs

    ! local
    integer                 :: n,i,dm,ng,comp
    integer                 :: ix,iy,iz,A_last,D_last
    integer                 :: lo(sold(1)%dim),hi(sold(1)%dim)
    real(dp_t), pointer     :: soop(:,:,:,:),snop(:,:,:,:),Iop(:,:,:,:)
    type(dpt_pntr) :: A(0:size(adv,dim=2)-1), D(0:size(diff,dim=2)-1) 


    dm    = sold(1)%dim

    A_last = size(adv,dim=2)-1
    D_last = size(diff,dim=2)-1
    
    ng = sold(1)%ng

    do n=1,nlevs
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle    
          call get_interp_pts_array(adv,A,n,i,A_last)
          call get_interp_pts_array(diff,D,n,i,D_last)
          soop    => dataptr(sold(n), i)         
          snop    => dataptr(snew(n), i)
          Iop     => dataptr(I_ADR(n), i)
          lo = lwb(get_box(sold(n), i))
          hi = upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)                     
                do ix = lo(1),hi(1)

                   do comp = 2, nscal
                      Iop(ix,iy,1,comp-1) = (snop(ix,iy,1,comp) - soop(ix,iy,1,comp))/dt&
                           + half*(D(3)%p(ix,iy,1,comp-1) + D(2)%p(ix,iy,1,comp-1))&
                           - D(1)%p(ix,iy,1,comp-1) - A(0)%p(ix,iy,1,comp-1)
                   enddo
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)

                      do comp = 2, nscal
                         Iop(ix,iy,iz,comp-1) = (snop(ix,iy,iz,comp) - soop(ix,iy,iz,comp))/dt&
                           + half*(D(3)%p(ix,iy,iz,comp-1) + D(2)%p(ix,iy,iz,comp-1))&
                           - D(1)%p(ix,iy,iz,comp-1) - A(0)%p(ix,iy,iz,comp-1)
                      enddo

                   end do
                end do
             end do

          end select
       end do

    end do

  end subroutine intgrl

  subroutine get_interp_pts_array(full_array,times_array,lev,box,hi)
    !*********************************
    ! this subroutine takes the full array of adv (or diff) data 
    ! and creates an array that holds adv (or diff) at all the various 
    ! times (i.e. adv at all the quadrature points) 
    ! associated with a single spatial location
    !*********************************

    type(multifab),  intent(in   ) :: full_array(:,0:)
    type(dpt_pntr),  intent(  out) :: times_array(0:)
    integer,         intent(in   ) :: lev,box,hi

    integer             :: t,s
    real(kind=dp_t), pointer :: aop(:,:,:,:)


    do t = 0, hi
       aop    => dataptr(full_array(lev,t), box)
       times_array(t)%p => aop
    end do

  end subroutine get_interp_pts_array

end module sdc_interpolation
