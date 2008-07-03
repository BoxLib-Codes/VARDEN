module rxns_integrator

  use bl_types
  use bl_constants_module
  use define_bc_module
  use ml_restriction_module
  use ml_layout_module
  use multifab_module
  use multifab_fill_ghost_module
  use multifab_physbc_module
  use probin_module,  only : nscal
  use reaction_fn

  implicit none

  public :: react,provisional,sdc

contains 

  ! used with strang splitting
  subroutine react(mla,the_bc_tower,s,dt,f,is_conserv)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(dp_t),      intent(in   ) :: dt
    interface
       subroutine f(soln,u,t,adv,diff)
         use bl_types
         real(dp_t),  intent(  out) :: soln(:)
         real(dp_t),  intent(in   ) :: u(:)
         real(dp_t),  intent(in   ) :: t
       end subroutine f
    end interface
    logical,         intent(in   )  :: is_converv

    ! local
    integer             :: n,i,dm,ng,comp,nlevs,ix,iy,iz
    integer             :: lo(s(1)%dim),hi(s(1)%dim)
    real(dp_t), pointer :: sop(:,:,:,:)
                               


    nlevs = mla%nlevel
    dm    = mla%dim
    
    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sop    => dataptr(s(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             do ix = lo(1),hi(1)
                do iy = lo(2), hi(2)                   
                   call vode(sop(ix,iy,1,:),dt,f,is_conserv)
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call vode(sop(ix,iy,iz,:),dt,f,is_converv)
                   end do
                end do
             end do

          end select
       end do

       ! fill ghost cells  
       call multifab_fill_boundary(s(n))       
       call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n))

    end do

    do n = nlevs,2, -1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1),&
                                      the_bc_tower%bc_tower_array(n),&
                                      1,dm+1,nscal)
    end do

  end subroutine react

  subroutine vode(u,dt,f,is_converv)
    real(kind=dp_t), intent(inout) ::  u(:)
    real(kind=dp_t), intent(in   ) ::  dt            
    interface
       subroutine f(soln,u,t)
         use bl_types
         real(kind=dp_t),  intent(  out)  :: soln(:)
         real(kind=dp_t),  intent(in   )  :: u(:)
         real(kind=dp_t),  intent(in   )  :: t
       end subroutine f
    end interface
    logical,          intent(in   ) :: is_converv

    integer i,s,ns
    real(kind=dp_t)  dtl, tl
!    real(kind=dp_t)  u_star(nspecies)
!    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
    real(kind=dp_t), allocatable ::  k1(:),k2(:),k3(:),k4(:)

    
    ! first component is density
    ns = nscal-1
    allocate(k1(ns), k2(ns), k3(ns),k4(ns))

    if (is_converv) then
       do i = 2,ns+1
          u(i) = u(i)/u(1)
       enddo
    endif

    dtl = dt/100.d0
    tl = zero

    do i = 1,100

       ! euler predictor-corrector method
!        rxn1 = f(u,tl,ns,c1,c2)
!        do s = 1,ns
!            u(s) = u(s) + dtl*(rxn1(s)) 
!            u_star(s) = u(s)+dtl*(rxn1(s))
!         enddo
!         rxn2 = f(u_star,tl,ns,c1,c2)
!         do s = 1,ns
!           u(s)=u(s)+dtl*(rxn1(s)+rxn2(s))/2.d0
!         enddo

       ! Runge-Kutta
       call f(k1,u,               tl           )
       call f(k2,u + half*dtl*k1, tl + half*dtl)
       call f(k3,u + half*dtl*k2, tl + half*dtl)
       call f(k4,u +      dtl*k3, tl +      dtl)
       
       do s = 1,ns
          u(s+1) = u(s+1) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/6.d0
       end do
       
       tl = tl + dtl
    enddo

    if (is_converv) then
       do i = 1,ns
          u(i+1) = s(i)*u(1)
       enddo
    endif
    
    return
    
  end subroutine vode


  ! used with SDC
  subroutine react(mla,the_bc_tower,s,dt,f,is_conserv,adv,diff)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(dp_t),      intent(in   ) :: dt
    interface
       subroutine f(soln,u,t,adv,diff)
         use bl_types
         real(dp_t),  intent(  out) :: soln(:)
         real(dp_t),  intent(in   ) :: u(:)
         real(dp_t),  intent(in   ) :: t
         real(dp_t),  intent(in   ) :: adv(:,:)
         real(dp_t),  intent(in   ) :: diff(:,:)
       end subroutine f
    end interface
    logical,         intent(in   )  :: is_converv
    type(multifab) , intent(in   )  :: adv(:,:), diff(:,:)

    ! local
    integer                 :: n,i,dm,ng,comp,nlevs
    integer                 :: ix,iy,iz,hi
    integer                 :: lo(s(1)%dim),hi(s(1)%dim)
    real(dp_t), pointer     :: sop(:,:,:,:)
    real(dp_t), allocatable :: A(:,:), D(:,:) 
                               


    nlevs = mla%nlevel
    dm    = mla%dim

    hi = size(diff,dim=1)-1
    allocate(A(0:hi,nscal-1),D(0:hi,nscal-1)
    
    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sop    => dataptr(s(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             do ix = lo(1),hi(1)
                do iy = lo(2), hi(2)          
                   call get_single_location_array(adv,A,n,i,ix,iy,1,hi)
                   call get_single_location_array(diff,D,n,i,ix,iy,1,hi)
                   call vode(sop(ix,iy,1,:),dt,f,is_conserv,A,D)
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call get_single_location_array(adv,A,n,i,ix,iy,iz,hi)
                      call get_single_location_array(diff,D,n,i,ix,iy,iz,hi)
                      call vode(sop(ix,iy,iz,:),dt,f,is_converv,A,D)
                   end do
                end do
             end do

          end select
       end do

       ! fill ghost cells  
       call multifab_fill_boundary(s(n))       
       call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n))

    end do

    do n = nlevs,2, -1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1),&
                                      the_bc_tower%bc_tower_array(n),&
                                      1,dm+1,nscal)
    end do


    deallocate(A)
    deallocate(D)

  end subroutine react

  subroutine vode(u,dt,f,is_converv,adv,diff)
    real(kind=dp_t), intent(inout) ::  u(:)
    real(kind=dp_t), intent(in   ) ::  dt            
    interface
       subroutine f(soln,u,t,adv,diff,dt)
         use bl_types
         real(kind=dp_t),  intent(  out)  :: soln(:)
         real(kind=dp_t),  intent(in   )  :: u(:)
         real(kind=dp_t),  intent(in   )  :: t,dt
         real(dp_t),       intent(in   ) :: adv(:,:)
         real(dp_t),       intent(in   ) :: diff(:,:)
       end subroutine f
    end interface
    logical,          intent(in   ) :: is_converv
    real(dp_t),       intent(in   ) :: adv(:,:)
    real(dp_t),       intent(in   ) :: diff(:,:)

    integer i,s,ns
    real(kind=dp_t)  dtl, tl
!    real(kind=dp_t)  u_star(nspecies)
!    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
    real(kind=dp_t), allocatable ::  k1(:),k2(:),k3(:),k4(:)

    
    ! first component is density
    ns = nscal-1
    allocate(k1(ns), k2(ns), k3(ns),k4(ns))

    if (is_converv) then
       do i = 2,ns+1
          u(i) = u(i)/u(1)
       enddo
    endif

    dtl = dt/100.d0
    tl = zero

    do i = 1,100

       ! euler predictor-corrector method
!        rxn1 = f(u,tl,ns,c1,c2)
!        do s = 1,ns
!            u(s) = u(s) + dtl*(rxn1(s)) 
!            u_star(s) = u(s)+dtl*(rxn1(s))
!         enddo
!         rxn2 = f(u_star,tl,ns,c1,c2)
!         do s = 1,ns
!           u(s)=u(s)+dtl*(rxn1(s)+rxn2(s))/2.d0
!         enddo

       ! Runge-Kutta
       call f(k1,u,               tl           ,adv,diff,dt)
       call f(k2,u + half*dtl*k1, tl + half*dtl,adv,diff,dt)
       call f(k3,u + half*dtl*k2, tl + half*dtl,adv,diff,dt)
       call f(k4,u +      dtl*k3, tl +      dtl,adv,diff,dt)
       
       do s = 1,ns
          u(s+1) = u(s+1) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/6.d0
       end do
       
       tl = tl + dtl
    enddo

    if (is_converv) then
       do i = 1,ns
          u(i+1) = s(i)*u(1)
       enddo
    endif

    deallocate(k1)
    deallocate(k2)
    deallocate(k3)
    deallocate(k4)
    
    return
    
  end subroutine vode

  ! Takes the full array of adv (or diff) data and creates an array that 
  ! holds adv (diff) at all the various times (i.e. adv at all the quadrature points) 
  ! associated with a single spatial location
  subroutine get_single_location_array(full_array,sngl_pt,lev,box,x,y,z,hi)
    type(mulitfab),   intent(in   ) :: full_array(:,:)
    real(kind=dp_t),  intent(  out) :: sngle_pt(:,:)
    integer,          intent(in   ) :: lev,box,hi
    real(kind=dp_t),  intent(in   ) :: x,y,z

    integer             :: t,s
    real(dp_t), pointer :: aop(:,:,:,:)


    do t = 0, hi
       aop    => dataptr(full_array(t,lev), box)
       do s = 1, nscal-1
          sngl_pt(t,s) = aop(x,y,z,s)
       end do
    end do

  end subroutine get_single_location_array


!--------------------------------------------------------------------
! make sure to adjust integral to 0 -> dt    


    ! for provisional soln in SDC
    subroutine provisional(soln,u,t,adv,diff,dt)
      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t,dt
      real(kind=dp_t),  intent(in   ) :: adv(:,:)
      real(kind=dp_t),  intent(in   ) :: diff(:,:)

      real(kind=dp_t),  intent(  out) :: rxns(:)


      ! compute the rxns
      call f_rxn(rxns,u,t)

      ! compute the adv-diff
      call AD_provisional(soln,adv,diff,t)

      soln = soln + rxns

    end subroutine provisional

    ! for computing SDC corrections
    subroutine sdc(soln,u,t,adv,diff,dt)
      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t,dt
      real(kind=dp_t),  intent(in   ) :: adv(:,:)
      real(kind=dp_t),  intent(in   ) :: diff(:,:)

      real(kind=dp_t),  intent(  out) :: rxns(:)
      integer,                        :: final,s


      final = size(A,dim=1)-1
      if(final+1 < 2) write(*,*)'ERROR: must have at least 2 values in ADV & DIFF'

      ! compute the rxns
      call f_rxn(rxns,u,t)

      ! compute the adv-diff
      call sdc_interpolant(soln,adv,diff,t,dt)

      do s = 1, nscal-1
         soln(s) = soln(s) + rxns(s) + dt*(adv(final,s) - adv(0,s) &
                                           + diff(final,s) - diff(0,s))

    end subroutine sdc
    
!-----------------------------------------------------------------------

  end module rxns_integrator
