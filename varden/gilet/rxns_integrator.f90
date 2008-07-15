module rxns_integrator

  use bl_types
  use bl_constants_module
  use define_bc_module
  use ml_restriction_module
  use ml_layout_module
  use multifab_module
  use multifab_fill_ghost_module
  use multifab_physbc_module
  use probin_module,  only : nscal, mass_fractions
  use sdc_interpolation
  use reaction_fn

  implicit none

  interface react
     module procedure react_strang
     module procedure react_sdc
  end interface

  public :: provisional,sdc

contains 

  ! used with strang splitting
  subroutine react_strang(mla,the_bc_tower,s,dx,dt,t,f)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t),      intent(in   ) :: dx(:,:)
    real(kind=dp_t),      intent(in   ) :: dt
    real(kind=dp_t),      intent(in   ) :: t
    interface
       subroutine f(soln,u,t)
         use bl_types
         real(kind=dp_t),  intent(  out) :: soln(:)
         real(kind=dp_t),  intent(in   ) :: u(:)
         real(kind=dp_t),  intent(in   ) :: t
       end subroutine f
    end interface

    ! local
    integer             :: n,i,dm,ng,nlevs,ix,iy,iz
    integer             :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
                               


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
                   call vode(sop(ix,iy,1,:),dt,f)
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call vode(sop(ix,iy,iz,:),dt,f)
                   end do
                end do
             end do

          end select
       end do

       ! fill ghost cells  
       call multifab_fill_boundary(s(n))       
       call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n),&
                            dx(n,:),t)

    end do

    do n = nlevs,2, -1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1),&
                                      the_bc_tower%bc_tower_array(n),&
                                      1,dm+1,nscal,dx(n-1:n,:),t)
    end do

  end subroutine react_strang

  subroutine vode(u,dt,f)
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

    integer i,s
    real(kind=dp_t)  dtl, tl
!    real(kind=dp_t)  u_star(nspecies)
!    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
    real(kind=dp_t), allocatable ::  k1(:),k2(:),k3(:),k4(:)

    
    allocate(k1(nscal), k2(nscal), k3(nscal),k4(nscal))
    k1(:)=zero
    k2(:)=zero
    k3(:)=zero
    k4(:)=zero


    ! either all the tracers are conserv or not convervative
    if (mass_fractions) then
       do i = 2,nscal
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
       
       do s = 2,nscal
          u(s) = u(s) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/6.d0
       end do
       
       tl = tl + dtl
    enddo

    ! either all the tracers are conserv or not convervative
    if (mass_fractions) then
       do i = 2,nscal
          u(i) = u(i)*u(1)
       enddo
    endif
    
    deallocate(k1,k2,k3,k4)

    return
    
  end subroutine vode


  ! used with SDC
  subroutine react_sdc(mla,the_bc_tower,s,dx,dt,t,f,adv,diff)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t),      intent(in   ) :: dx(:,:)
    real(kind=dp_t),      intent(in   ) :: dt
    real(kind=dp_t),      intent(in   ) :: t
    interface
       subroutine f(soln,u,t,adv,diff,dt)
         use bl_types
         real(kind=dp_t),  intent(  out) :: soln(:)
         real(kind=dp_t),  intent(in   ) :: u(:)
         real(kind=dp_t),  intent(in   ) :: t,dt
         real(kind=dp_t),  intent(in   ) :: adv(:,:)
         real(kind=dp_t),  intent(in   ) :: diff(:,:)
       end subroutine f
    end interface
    type(multifab) , intent(in   )  :: adv(:,:), diff(:,:)

    ! local
    integer                 :: n,i,dm,ng,nlevs
    integer                 :: ix,iy,iz,D_last,A_last
    integer                 :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer     :: sop(:,:,:,:)
    real(kind=dp_t), allocatable :: A(:,:), D(:,:) 
                               


    nlevs = mla%nlevel
    dm    = mla%dim

    A_last = size(adv,dim=1)-1
    allocate(A(0:A_last,nscal-1))
    D_last = size(diff,dim=1)-1
    allocate(D(0:D_last,nscal-1))
    
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
                   call get_single_location_array(adv,A,n,i,ix,iy,1,A_last)
                   call get_single_location_array(diff,D,n,i,ix,iy,1,D_last)
                   call vode_sdc(sop(ix,iy,1,:),dt,f,A,D)
                end do
             end do

          case (3)
             do ix = lo(1), hi(1)
                do iy = lo(2), hi(2)
                   do iz = lo(3), hi(3)
                      call get_single_location_array(adv,A,n,i,ix,iy,iz,A_last)
                      call get_single_location_array(diff,D,n,i,ix,iy,iz,D_last)
                      call vode_sdc(sop(ix,iy,iz,:),dt,f,A,D)
                   end do
                end do
             end do

          end select
       end do

       ! fill ghost cells  
       call multifab_fill_boundary(s(n))       
       call multifab_physbc(s(n),1,dm+1,nscal,the_bc_tower%bc_tower_array(n),&
                            dx(n,:),t)

    end do

    do n = nlevs,2, -1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
       call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                      the_bc_tower%bc_tower_array(n-1),&
                                      the_bc_tower%bc_tower_array(n),&
                                      1,dm+1,nscal,dx(n-1:n,:),t)
    end do


    deallocate(A)
    deallocate(D)

  end subroutine react_sdc

  subroutine vode_sdc(u,dt,f,adv,diff)
    real(kind=dp_t), intent(inout) ::  u(:)
    real(kind=dp_t), intent(in   ) ::  dt            
    interface
       subroutine f(soln,u,t,adv,diff,dt)
         use bl_types
         real(kind=dp_t),  intent(  out) :: soln(:)
         real(kind=dp_t),  intent(in   ) :: u(:)
         real(kind=dp_t),  intent(in   ) :: t,dt
         real(kind=dp_t),  intent(in   ) :: adv(:,:)
         real(kind=dp_t),  intent(in   ) :: diff(:,:)
       end subroutine f
    end interface
    real(kind=dp_t),       intent(in   ) :: adv(:,:)
    real(kind=dp_t),       intent(in   ) :: diff(:,:)

    integer i,s
    real(kind=dp_t)  dtl, tl
!    real(kind=dp_t)  u_star(nspecies)
!    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
    real(kind=dp_t), allocatable ::  k1(:),k2(:),k3(:),k4(:)

    
    allocate(k1(nscal), k2(nscal), k3(nscal), k4(nscal))
    k1(:) = zero
    k2(:) = zero
    k3(:) = zero
    k4(:) = zero

    ! either all the tracers are conserv or not convervative
    ! change to use a mix?
    if (mass_fractions) then
       do i = 2,nscal
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
       
       do s = 2,nscal
          u(s) = u(s) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/6.d0
       end do
       
       tl = tl + dtl
    enddo

    ! either all the tracers are conserv or not convervative
    ! change to use a mix?
    if (mass_fractions) then
       do i = 2,nscal
          u(i) = u(i)*u(1)
       enddo
    endif

    deallocate(k1)
    deallocate(k2)
    deallocate(k3)
    deallocate(k4)
    
    return
    
  end subroutine vode_sdc


!--------------------------------------------------------------------
! make sure to adjust integral to 0 -> dt    


    ! for provisional soln in SDC
    subroutine provisional(soln,u,t,adv,diff,dt)
      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t,dt
      real(kind=dp_t),  intent(in   ) :: adv(:,:)
      real(kind=dp_t),  intent(in   ) :: diff(:,:)

      real(kind=dp_t) :: rxns(nscal)

      rxns(:) = zero
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

      real(kind=dp_t) :: rxns(nscal)


      rxns(:) = zero
      ! compute the rxns
      call f_rxn(rxns,u,t)

      ! compute the adv-diff
      call sdc_interpolant(soln,adv,diff,t,dt)

      soln = soln + rxns 
      
    end subroutine sdc
    
!-----------------------------------------------------------------------

  end module rxns_integrator
