module sdc_interpolation

  use bl_types
  use multifab_module
  use bl_constants_module
  use probin_module, only : nscal, nspec

  implicit none
  
  ! number of terms at various times we need to save
  integer, parameter :: n_adv = 1
  integer, parameter :: n_diff = 4
  integer, parameter :: n_interp_pts = 2
  real(kind=dp_t), parameter :: sdc_dt_fac(n_interp_pts) = (/THIRD,TWO3RD/)
  real(kind=dp_t), parameter :: pts(0:n_interp_pts) = (/ZERO,THIRD,ONE/)


  type dpt_pntr
     real(kind=dp_t), dimension(:), pointer :: p(:,:,:,:)
  end type dpt_pntr

  public ::  mk_I_AD, get_single_location_array,&
             get_interp_pts_array

contains
  
  subroutine mk_I_AD(I_AD,adv,diff,nlevs)
    !*************************************************************** 
    ! This takes is really I_AD/(sdc_fac*dt) with A canceled with A's from I_R
    !***************************************************************
    type(multifab) , intent(inout) :: I_AD(:,:)
    type(multifab) , intent(in   ) :: adv(:,:), diff(:,:)
    integer,         intent(in   ) :: nlevs

     ! local
    integer                 :: n,i,j,ng,comp,dm
    integer                 :: ix,iy,iz
    integer                 :: lo(I_AD(1,1)%dim),hi(I_AD(1,1)%dim)
    type(dpt_pntr) :: A(n_adv), D(n_diff), Intgrl(n_interp_pts+1) 


    ng = I_AD(1,1)%ng
    dm = I_AD(1,1)%dim

    do n=1,nlevs
       do i = 1, I_AD(n,1)%nboxes
          if ( multifab_remote(I_AD(n,1), i) ) cycle      
          call get_interp_pts_array(adv, A     ,n,i,n_adv)
          call get_interp_pts_array(diff,D     ,n,i,n_diff)
          call get_interp_pts_array(I_AD,Intgrl,n,i,n_interp_pts+1)
          lo = lwb(get_box(I_AD(n,1), i))
          hi = upb(get_box(I_AD(n,1), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)

                      do comp = 1, nspec
                         Intgrl(1)%p(ix,iy,1,comp) = &
                              (FIVE*D(1)%p(ix,iy,1,comp)-D(2)%p(ix,iy,1,comp))/&
                              FOUR
                         Intgrl(2)%p(ix,iy,1,comp) = &
                              HALF*(D(1)%p(ix,iy,1,comp)+D(2)%p(ix,iy,1,comp))
                         Intgrl(3)%p(ix,iy,1,comp) = &
                              (THREE*D(1)%p(ix,iy,1,comp)+D(2)%p(ix,iy,1,comp))/&
                              FOUR
                      enddo
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 1, nspec
                         Intgrl(1)%p(ix,iy,1,comp) = &
                              (THREE*D(1)%p(ix,iy,iz,comp)+D(2)%p(ix,iy,iz,comp))/&
                              FOUR
                         Intgrl(2)%p(ix,iy,iz,comp) = &
                              HALF*(D(1)%p(ix,iy,iz,comp)+D(2)%p(ix,iy,iz,comp))
                         Intgrl(3)%p(ix,iy,1,comp) = sdc_dt_fac(1)* &
                              Intgrl(1)%p(ix,iy,iz,comp) + &
                              sdc_dt_fac(2)*Intgrl(2)%p(ix,iy,iz,comp)
                      enddo

                   end do
                end do
             end do

          end select
       end do
    end do   

  end subroutine mk_I_AD

  subroutine get_single_location_array(full_array,sngl_pt,lev,box,x,y,z,hi)
    !*********************************
    ! this subroutine takes the full array of adv (or diff) data 
    ! and creates an array that holds adv (diff) at all the various times 
    ! (i.e. adv at all the quadrature points) 
    ! associated with a single spatial location
    !*********************************

    type(multifab),   intent(in   ) :: full_array(:,0:)
    real(kind=dp_t),  intent(  out) :: sngl_pt(0:,:)
    integer,          intent(in   ) :: lev,box,hi
    integer,          intent(in   ) :: x,y,z

    integer             :: t,s
    real(kind=dp_t), pointer :: aop(:,:,:,:)


    do t = 0, hi
       aop    => dataptr(full_array(lev,t), box)
       do s = 1, nspec
          sngl_pt(t,s) = aop(x,y,z,s)
       end do
    end do

  end subroutine get_single_location_array

  subroutine get_interp_pts_array(full_array,times_array,lev,box,hi)
    !*********************************
    ! this subroutine takes the full array of adv (or diff) data 
    ! and creates an array that holds adv (diff) at all the various times 
    ! (i.e. adv at all the quadrature points) 
    ! associated with a single spatial location
    !*********************************

    type(multifab),  intent(in   ) :: full_array(:,:)
    type(dpt_pntr),  intent(  out) :: times_array(:)
    integer,         intent(in   ) :: lev,box,hi

    integer             :: t,s
    real(kind=dp_t), pointer :: aop(:,:,:,:)


    do t = 1, hi
       aop    => dataptr(full_array(lev,t), box)
       times_array(t)%p => aop
    end do

  end subroutine get_interp_pts_array

end module sdc_interpolation
