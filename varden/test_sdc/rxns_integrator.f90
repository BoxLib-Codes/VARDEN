module rxns_integrator

  use bl_types
  use bl_constants_module
  use define_bc_module
  use ml_restriction_module
  use ml_layout_module
  use multifab_module
  use multifab_fill_ghost_module
  use multifab_physbc_module
  use probin_module,  only : nscal, mass_fractions, n_rxn_steps
  use sdc_interpolation
!  use reaction_fn

  implicit none

  interface react
     module procedure react_strang
     module procedure react_sdc
  end interface


contains 

  ! used with strang splitting
  subroutine react_strang(mla,the_bc_tower,s,dx,dt,t)!,f)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t),      intent(in   ) :: dx(:,:)
    real(kind=dp_t),      intent(in   ) :: dt
    real(kind=dp_t),      intent(in   ) :: t
!     interface
!        subroutine f(soln,u,t)
!          use bl_types
!          real(kind=dp_t),  intent(  out) :: soln(:)
!          real(kind=dp_t),  intent(in   ) :: u(:)
!          real(kind=dp_t),  intent(in   ) :: t
!        end subroutine f
!     end interface

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
             do iy = lo(2), hi(2)                   
                do ix = lo(1),hi(1)
                   call vode(sop(ix,iy,1,:))!,f)
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
                      call vode(sop(ix,iy,iz,:))!,f)
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


  contains

    subroutine vode(u)!,f)
      real(kind=dp_t), intent(inout) ::  u(:)
   !    interface
!          subroutine f(soln,u,t)
!            use bl_types
!            real(kind=dp_t),  intent(  out)  :: soln(:)
!            real(kind=dp_t),  intent(in   )  :: u(:)
!            real(kind=dp_t),  intent(in   )  :: t
!          end subroutine f
!       end interface

      integer i,s
      real(kind=dp_t)  dtl, tl
      !    real(kind=dp_t)  u_star(nspecies)
      !    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
      real(kind=dp_t)  k1(nscal),k2(nscal),k3(nscal),k4(nscal)


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

      dtl = dt/dble(n_rxn_steps)
      tl = zero

      do i = 1,n_rxn_steps

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
         call f_rxn(k1,u,               tl           )
         call f_rxn(k2,u + half*dtl*k1, tl + half*dtl)
         call f_rxn(k3,u + half*dtl*k2, tl + half*dtl)
         call f_rxn(k4,u +      dtl*k3, tl +      dtl)

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

      return

    end subroutine vode

    subroutine f_rxn(soln,u,t)  
      use probin_module, only : k_rxn1, k_rxn2

      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t

      real(kind=dp_t) :: u2,u3,u4

      ! first component is density 
      u2 = merge(u(2), zero, u(2) > 0.d0)
      u3 = merge(u(3), zero, u(3) > 0.d0)
      u4 = merge(u(4), zero, u(4) > 0.d0)

      soln(2) =    -k_rxn1*u3*u2 + half*k_rxn2*u4
      soln(3) =    -k_rxn1*u3*u2 + half*k_rxn2*u4
      soln(4) = two*k_rxn1*u3*u2 -      k_rxn2*u4

    end subroutine f_rxn


  end subroutine react_strang



  ! used with SDC
  subroutine react_sdc(mla,the_bc_tower,s,dx,dt,t,adv,adv_rho,diff,sdc_flag)  

    type(ml_layout), intent(in   ) :: mla   
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t),      intent(in   ) :: dx(:,:)
    real(kind=dp_t),      intent(in   ) :: dt
    real(kind=dp_t),      intent(in   ) :: t
    type(multifab) , intent(in   )  :: adv(:,:), diff(:,:)
    type(multifab) , intent(in   )  :: adv_rho(:)
    integer        , intent(in   )  :: sdc_flag

    ! local
    integer                 :: n,i,dm,ng,nlevs
    integer                 :: ix,iy,iz,D_last,A_last
    integer                 :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer     :: sop(:,:,:,:),arop(:,:,:,:)
!    real(kind=dp_t) :: A(0:size(adv,dim=2)-1,nspec), D(0:size(diff,dim=2)-1,nspec)
    type(dpt_pntr) :: A(0:size(adv,dim=2)-1), D(0:size(diff,dim=2)-1)         


    nlevs = mla%nlevel
    dm    = mla%dim

    A_last = size(adv,dim=2)-1
    D_last = size(diff,dim=2)-1
   
    ng = s(1)%ng
    iz = 1

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          call get_interp_pts_array(adv,A,n,i,A_last)
          call get_interp_pts_array(diff,D,n,i,D_last)
          if(mass_fractions) arop => dataptr(adv_rho(n),i)
          sop    => dataptr(s(n), i)
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)    
                do ix = lo(1),hi(1)
 !                  call get_times_array(adv,A,n,i,A_last)
 !                  call get_times_array(diff,D,n,i,D_last)
                   call vode_sdc(sop(ix,iy,iz,:))
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
 !                     call get_times_array(adv,A,n,i,A_last)
 !                     call get_times_array(diff,D,n,i,D_last)
                      call vode_sdc(sop(ix,iy,iz,:))
!                     call vode_sdc(sop(ix,iy,A(:,ix,iy,iz,:),D(:,ix,iy,iz,:))
                      ! call vode_sdc(sop(ix,iy,iz,:),dt,provisional)!,A,D)
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
 
  contains
  
    subroutine vode_sdc(u)
    !subroutine vode_sdc(u,adv,diff)
      real(kind=dp_t), intent(inout) ::  u(:)
     ! real(kind=dp_t),       intent(in   ) :: adv(:,:)
      !real(kind=dp_t),       intent(in   ) :: diff(:,:)

      integer i,s
      real(kind=dp_t)  dtl, tl
      !    real(kind=dp_t)  u_star(nspecies)
      !    real(kind=dp_t)  rxn1(nspecies), rxn2(nspecies)
      real(kind=dp_t) ::  k1(nscal), k2(nscal), k3(nscal), k4(nscal)


!unneeded
!      k1(:) = zero
!      k2(:) = zero
!      k3(:) = zero
!      k4(:) = zero


      dtl = dt/dble(n_rxn_steps)
      tl = zero

      select case (sdc_flag)

      ! for provisional integration   
      case(1)    
         do i = 1,n_rxn_steps

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
            call provisional(k1,u              , tl           )
            call provisional(k2,u + half*dtl*k1, tl + half*dtl)
            call provisional(k3,u + half*dtl*k2, tl + half*dtl)
            call provisional(k4,u +      dtl*k3, tl +      dtl)

!            call provisional(k1,u              , tl           ,adv,diff)
!            call provisional(k2,u + half*dtl*k1, tl + half*dtl,adv,diff)
!            call provisional(k3,u + half*dtl*k2, tl + half*dtl,adv,diff)
!            call provisional(k4,u +      dtl*k3, tl +      dtl,adv,diff)
            
            do s = 1,nscal
               u(s) = u(s) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/six
            end do

            tl = tl + dtl
         enddo

      ! for SDC integration   
      case(2)    
         do i = 1,n_rxn_steps

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
            call sdc(k1,u              , tl           )
            call sdc(k2,u + half*dtl*k1, tl + half*dtl)
            call sdc(k3,u + half*dtl*k2, tl + half*dtl)
            call sdc(k4,u +      dtl*k3, tl +      dtl)

!            call sdc(k1,u              , tl           ,adv,diff)
!            call sdc(k2,u + half*dtl*k1, tl + half*dtl,adv,diff)
!            call sdc(k3,u + half*dtl*k2, tl + half*dtl,adv,diff)
!            call sdc(k4,u +      dtl*k3, tl +      dtl,adv,diff)
            
            do s = 1,nscal
               u(s) = u(s) + dtl*(k1(s) + 2.d0*k2(s) + 2.d0*k3(s) + k4(s))/six
            end do

            tl = tl + dtl
         enddo

      end select

      return
    end subroutine vode_sdc


    !--------------------------------------------------------------------
    ! The integrals
    ! make sure to adjust integral to 0 -> dt        
  
    ! for provisional soln in SDC
    subroutine provisional(soln,u,t)!,adv,diff)
      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t!,dt
!      real(kind=dp_t),  intent(in   ) :: adv(:,:)
!      real(kind=dp_t),  intent(in   ) :: diff(:,:)

      integer         :: i
      real(kind=dp_t) :: rxns(nscal)

!      rxns(:) = zero
      ! compute the rxns 
      call f_rxn(rxns,u,t)
      soln = rxns

      ! compute the adv-diff for the provisional solution
      ! depends only on choice of provisional method
      do i = 1, nspec
         soln(i+1) = soln(i+1)+ A(0)%p(ix,iy,iz,i) + D(1)%p(ix,iy,iz,i)
!         soln(i+1) = adv(0,i) + diff(1,i)
      end do

    end subroutine provisional


    ! for computing SDC corrections
    subroutine sdc(soln,u,t)!,adv,diff)
      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t!,dt
!      real(kind=dp_t),  intent(in   ) :: adv(:,:)
!      real(kind=dp_t),  intent(in   ) :: diff(:,:)

      integer         :: i
      real(kind=dp_t) :: rxns(nscal)


!      rxns(:) = zero

      ! compute the rxns
      call f_rxn(rxns,u,t)
      soln = rxns

      ! compute the adv-diff: 
      ! (interpolation of A+D) + (D_new-D_old) + (A_new-A_old)
      do i = 1, nspec
           soln(i+1) = soln(i+1) + A(0)%p(ix,iy,iz,i) + D(0)%p(ix,iy,iz,i)& 
                       - D(2)%p(ix,iy,iz,i) + D(1)%p(ix,iy,iz,i)& 
                       + t*(D(2)%p(ix,iy,iz,i) - D(0)%p(ix,iy,iz,i))/dt
      end do

    end subroutine sdc

    !-----------------------------------------------------------------------
 
    subroutine f_rxn(soln,u,t)  
     use probin_module, only : k_rxn1, k_rxn2

      real(kind=dp_t),  intent(  out) :: soln(:)
      real(kind=dp_t),  intent(in   ) :: u(:)
      real(kind=dp_t),  intent(in   ) :: t

      real(kind=dp_t) :: u2,u3,u4,k1,k2,k3
      real(kind=dp_t) :: rho

      ! first component is density 
       u2 = merge(u(2), zero, u(2) > zero)
       u3 = merge(u(3), zero, u(3) > zero)
       u4 = merge(u(4), zero, u(4) > zero)

!       k1 = merge(k_rxn2,zero,u(2)<= one)
!       k2 = merge(k_rxn2,zero,u(3)<= one)
!       k3 = merge(k_rxn1,zero,u(4)<= one)

      u2 = u(2)
      u3 = u(3)
      u4 = u(4)

      k1 = k_rxn2
      k2 = k_rxn2
      k3 = k_rxn1

      if(mass_fractions) then
         rho = u(1)
         soln(1) = arop(ix,iy,iz,1)
         soln(2) = -k_rxn1*u3*u2/rho + half*k1*u4
         soln(3) = -k_rxn1*u3*u2/rho + half*k2*u4
         soln(4) =  two*k3*u3*u2/rho - k_rxn2*u4
      else
         soln(1) = zero
         soln(2) = -k_rxn1*u3*u2 + half*k1*u4
         soln(3) = -k_rxn1*u3*u2 + half*k2*u4
         soln(4) =  two*k3*u3*u2 - k_rxn2*u4
      end if

    end subroutine f_rxn

  end subroutine react_sdc

end module rxns_integrator
