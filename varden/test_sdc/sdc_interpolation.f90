module sdc_interpolation

  use bl_types
  use multifab_module
  use bl_constants_module
  use probin_module, only : nscal

  implicit none
  
  ! number of terms at various times we need to save
  integer :: n_adv = 1
  integer :: n_diff = 4

  public :: AD_provisional, sdc_interpolant, provisional_intgrl,&
            intgrl,get_single_location_array

contains

  subroutine sdc_interpolant(soln,adv,diff,t,dt)
    !*******************************************
    ! this subroutine computes 
    ! (interpolation of A+D) + (D_new-D_old) + (A_new-A_old)
    ! used in integrating rxns
    !*******************************************

    real(dp_t)    , intent(inout)  :: soln(:)
    real(dp_t)    , intent(in   )  :: adv(0:,:)
    real(dp_t)    , intent(in   )  :: diff(0:,:)
    real(dp_t)    , intent(in   )  :: t,dt

    ! local
    integer             :: i

    do i = 2, nscal
       soln(i) = adv(0,i) + diff(0,i) - diff(2,i) + diff(1,i) &
                 + t*(diff(2,i)-diff(0,i))/dt
    end do

  end subroutine sdc_interpolant

  subroutine provisional_intgrl(I_ADR,sold,snew,adv,diff,dt,nlevs)
    !*******************************************************
    ! this subroutine computes:
    ! I  {= intgrl_rxns (= s_new-s_old - dt*(provis approx to AD))
    !       + intgrl(interpol of AD)} 
    !   + dt*Adv   
    ! used as a source for computing advection in first iter of SDC
    !*******************************************************

    !   type(ml_layout), intent(in   ) :: mla   
    !   type(bc_tower) , intent(in   ) :: the_bc_tower
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
    real(dp_t), allocatable :: A(:,:), D(:,:)
                               

    ! diff(0,:) = D(s_n)
    ! diff(1,:) = D(s_n+1) w/o rxns
    ! diff(2,:) = D(s_n+1) (k-1)th sdc iteration
    ! diff(3,:) = D(s_n+1) kth sdc iteration

    dm    = sold(1)%dim

    A_last = size(adv,dim=1)-1
    allocate(A(0:A_last,nscal-1))
    D_last = size(diff,dim=1)-1
    allocate(D(0:D_last,nscal-1))
    
    ng = sold(1)%ng

    do n=1,nlevs
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle
          soop    => dataptr(sold(n), i)         
          snop    => dataptr(snew(n), i)
          Iop     => dataptr(I_ADR(n), i)
          lo = lwb(get_box(sold(n), i))
          hi = upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             do ix = lo(1),hi(1)
                do iy = lo(2), hi(2)          
                   call get_single_location_array(adv,A,n,i,ix,iy,1,A_last)
                   call get_single_location_array(diff,D,n,i,ix,iy,1,D_last)

                   do comp = 2, nscal
                      Iop(ix,iy,1,comp-1) = snop(ix,iy,1,comp) - soop(ix,iy,1,comp)&
                           + half*dt*(D(2,comp-1) + D(0,comp-1)) &
                           - dt*(D(1,comp-1) + A(0,comp-1))
                   enddo
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call get_single_location_array(adv,A,n,i,ix,iy,iz,A_last)
                      call get_single_location_array(diff,D,n,i,ix,iy,iz,D_last)
    
                      do comp = 2, nscal
                         Iop(ix,iy,1,comp-1) = snop(ix,iy,iz,comp) - soop(ix,iy,iz,comp)&
                              + half*dt*(D(2,comp-1) + D(0,comp-1)) &
                              - dt*(D(1,comp-1) + A(0,comp-1))
                      enddo

                   end do
                end do
             end do

          end select
       end do

       ! shouldn't have to worry about ghost cells b/c get filled
       ! in mkscalforce sudroutine?
       ! fill ghost cells  
       !call multifab_fill_boundary(s(n))       
       !call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n))

    end do

!    do n = nlevs,2, -1
!       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
!       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
!                                      the_bc_tower%bc_tower_array(n-1),&
!                                      the_bc_tower%bc_tower_array(n),&
!                                      1,dm+1,nscal)
!    end do

    deallocate(A)
    deallocate(D)

  end subroutine provisional_intgrl


  subroutine intgrl(I_ADR,sold,snew,adv,diff,dt,nlevs)
    !*******************************************************
    ! this subroutine computes:
    ! I  {= intgrl_rxns + intgrl(interpol of AD)} 
    !   + dt*Adv   
    ! used as a source for computing advection in SDC
    !*******************************************************

    !   type(ml_layout), intent(in   ) :: mla   
    !   type(bc_tower) , intent(in   ) :: the_bc_tower
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
    real(dp_t), allocatable :: A(:,:), D(:,:)
                               


    ! diff(0,:) = D(s_n)
    ! diff(1,:) = D(s_n+1) w/o rxns
    ! diff(2,:) = D(s_n+1) (k-1)th sdc iteration
    ! diff(3,:) = D(s_n+1) kth sdc iteration

    dm    = sold(1)%dim

    A_last = size(adv,dim=1)-1
    allocate(A(0:A_last,nscal-1))
    D_last = size(diff,dim=1)-1
    allocate(D(0:D_last,nscal-1))
    
    ng = sold(1)%ng

    do n=1,nlevs
       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n), i) ) cycle
          soop    => dataptr(sold(n), i)         
          snop    => dataptr(snew(n), i)
          Iop     => dataptr(I_ADR(n), i)
          lo = lwb(get_box(sold(n), i))
          hi = upb(get_box(sold(n), i))
          select case (dm)
          case (2)
             do ix = lo(1),hi(1)
                do iy = lo(2), hi(2)          
                   call get_single_location_array(adv,A,n,i,ix,iy,1,A_last)
                   call get_single_location_array(diff,D,n,i,ix,iy,1,D_last)

                   do comp = 2, nscal
                      Iop(ix,iy,1,comp-1) = snop(ix,iy,1,comp) - soop(ix,iy,1,comp)&
                           + half*dt*(D(1,comp-1) - D(2,comp-1)) &
                           - dt*(D(1,comp-1) + A(0,comp-1))
                   enddo
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call get_single_location_array(adv,A,n,i,ix,iy,iz,A_last)
                      call get_single_location_array(diff,D,n,i,ix,iy,iz,D_last)

                      do comp = 2, nscal
                         Iop(ix,iy,1,comp-1) = snop(ix,iy,iz,comp) - soop(ix,iy,iz,comp)&
                              + half*dt*(D(3,comp-1) + D(2,comp-1)) &
                              - dt*(D(1,comp-1) + A(0,comp-1))
                      enddo

                   end do
                end do
             end do

          end select
       end do

       ! shouldn't have to worry about ghost cells b/c get filled
       ! in mkscalforce sudroutine?
       ! fill ghost cells  
       !call multifab_fill_boundary(s(n))       
       !call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n))

    end do

!    do n = nlevs,2, -1
!       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
!       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
!                                      the_bc_tower%bc_tower_array(n-1),&
!                                      the_bc_tower%bc_tower_array(n),&
!                                      1,dm+1,nscal)
!    end do


    deallocate(A)
    deallocate(D)

  end subroutine intgrl


  subroutine AD_provisional(soln,adv,diff,t)
    !***************************************
    ! this subroutine is used in the last step of calculating 
    ! the provisional solution
    ! depends only on choice of provisional method
    !***************************************

    real(dp_t)    , intent(inout)  :: soln(:)
    real(dp_t)    , intent(in   )  :: adv(0:,:)
    real(dp_t)    , intent(in   )  :: diff(0:,:)
    real(dp_t)    , intent(in   )  :: t

    ! local
    integer             :: i

    do i = 2, nscal
       soln(i) = adv(0,i) + diff(1,i)
    end do

  end subroutine AD_provisional

  subroutine get_single_location_array(full_array,sngl_pt,lev,box,x,y,z,hi)
    !*********************************
    ! this subroutine takes the full array of adv (or diff) data 
    ! and creates an array that holds adv (diff) at all the various times 
    ! (i.e. adv at all the quadrature points) 
    ! associated with a single spatial location
    !*********************************

    type(multifab),   intent(in   ) :: full_array(0:,:)
    real(kind=dp_t),  intent(  out) :: sngl_pt(0:,:)
    integer,          intent(in   ) :: lev,box,hi
    integer,          intent(in   ) :: x,y,z

    integer             :: t,s
    real(kind=dp_t), pointer :: aop(:,:,:,:)


    do t = 0, hi
       aop    => dataptr(full_array(t,lev), box)
       do s = 1, nscal-1
          sngl_pt(t,s) = aop(x,y,z,s)
       end do
    end do

  end subroutine get_single_location_array

end module sdc_interpolation
