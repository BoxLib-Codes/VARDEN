module update_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module
  use probin_module, only : mass_fractions, sdc_iters, nscal

  implicit none

  private

  public :: update

contains

!   subroutine mk_adv(mla,umac,sedge,flux,adv_s,dx,&
!        dt,t,is_vel,is_cons,the_bc_level,adv_rho)

!     use bl_constants_module
!     use multifab_physbc_module
!     use ml_restriction_module, only: ml_cc_restriction
!     use multifab_fill_ghost_module

!     type(ml_layout)   , intent(in   ) :: mla
!     type(multifab)    , intent(in   ) :: umac(:,:)
!     type(multifab)    , intent(in   ) :: sedge(:,:)
!     type(multifab)    , intent(in   ) :: flux(:,:)
!     type(multifab)    , intent(inout) :: adv_s(:)
!     real(kind = dp_t) , intent(in   ) :: dx(:,:),dt,t
!     logical           , intent(in   ) :: is_vel,is_cons(:)
!     type(bc_level)    , intent(in   ) :: the_bc_level(:)
!     type(multifab)    , intent(inout), optional :: adv_rho(:)

!     ! local
!     real(kind=dp_t), pointer :: ap(:,:,:,:)
!     real(kind=dp_t), pointer :: arhop(:,:,:,:)
!     real(kind=dp_t), pointer :: ump(:,:,:,:)
!     real(kind=dp_t), pointer :: vmp(:,:,:,:)
!     real(kind=dp_t), pointer :: wmp(:,:,:,:)
!     real(kind=dp_t), pointer :: sepx(:,:,:,:)
!     real(kind=dp_t), pointer :: sepy(:,:,:,:)
!     real(kind=dp_t), pointer :: sepz(:,:,:,:)
!     real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
!     real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
!     real(kind=dp_t), pointer :: fluxpz(:,:,:,:)

!     integer :: lo(adv_s(1)%dim),hi(adv_s(1)%dim)
!     integer :: i,ng,dm,n,nlevs

!     nlevs = mla%nlevel
!     dm    = mla%dim

!     ng = adv_s(1)%ng

!     do n=1,nlevs

!        do i = 1, adv_s(n)%nboxes
!           if ( multifab_remote(adv_s(n),i) ) cycle
!           ap     => dataptr(adv_s(n),i)
!           ump    => dataptr(umac(n,1),i)
!           vmp    => dataptr(umac(n,2),i)
!           sepx   => dataptr(sedge(n,1),i)
!           sepy   => dataptr(sedge(n,2),i)
!           fluxpx => dataptr(flux(n,1),i)
!           fluxpy => dataptr(flux(n,2),i)
!           lo = lwb(get_box(adv_s(n),i))
!           hi = upb(get_box(adv_s(n),i))
!           if (mass_fractions .AND. (.not. is_vel) .AND. sdc_iters >= 0) then
!              arhop => dataptr(adv_rho(n),i)
!              select case (dm)
!              case (2)
!                 call mkadv_2d( ump(:,:,1,1), vmp(:,:,1,1), &
!                      sepx(:,:,1,:), sepy(:,:,1,:), &
!                      fluxpx(:,:,1,:), fluxpy(:,:,1,:), ap(:,:,1,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons, arhop(:,:,1,:))
!              case (3)
!                 wmp    => dataptr( umac(n,3),i)
!                 sepz   => dataptr(sedge(n,3),i)
!                 fluxpz => dataptr( flux(n,3),i)
!                 call mkadv_3d(sop(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
!                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
!                      fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
!                      ap(:,:,:,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons,arhop(:,:,:,:))
!              end select
!           else !adv_rho not needed without mass-fractions
!              select case (dm)
!              case (2)
!                 call mkadv_2d(ump(:,:,1,1), vmp(:,:,1,1), &
!                      sepx(:,:,1,:), sepy(:,:,1,:), &
!                      fluxpx(:,:,1,:), fluxpy(:,:,1,:), ap(:,:,1,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
!              case (3)
!                 wmp    => dataptr( umac(n,3),i)
!                 sepz   => dataptr(sedge(n,3),i)
!                 fluxpz => dataptr( flux(n,3),i)
!                 call mkadv_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
!                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
!                      fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
!                      ap(:,:,:,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
!              end select
!           end if
!        end do

!     enddo ! end loop over levels

! contains

!   subroutine mkadv_2d(umac,vmac,sedgex,sedgey,fluxx,fluxy,&
!        adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

!     use bl_constants_module

!     integer           , intent(in   ) :: lo(:), hi(:), ng
!     real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,:)  
!     real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
!     real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
!     real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
!     real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,:) 
!     real (kind = dp_t), intent(in   ) :: dx(:)
!     real (kind = dp_t), intent(in   ) :: dt
!     logical           , intent(in   ) :: is_vel
!     logical           , intent(in   ) :: is_cons(:)
!     real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1):,lo(2):,:)  

!     integer :: i, j, comp
!     real (kind = dp_t) ubar,vbar
!     real (kind = dp_t) ugradu,ugradv,ugrads
!     real (kind = dp_t) :: divsu

!     if (.not. is_vel) then

!        do comp = 1,nscal
!           if (is_cons(comp)) then
!              do j = lo(2), hi(2)
!                 do i = lo(1), hi(1)
!                    divsu = (fluxx(i+1,j,comp)-fluxx(i,j,comp))/dx(1) &
!                          + (fluxy(i,j+1,comp)-fluxy(i,j,comp))/dx(2)
!                    if (comp > 1) then; adv_s(i,j,comp-1) = -divsu
!                    else
!                       if (mass_fractions .AND.sdc_iters >= 0) adv_rho(i,j,1) = -divsu
!                    endif
!                 enddo
!              enddo
!           else
!              do j = lo(2), hi(2)
!                 do i = lo(1), hi(1)
!                    ubar = HALF*(umac(i,j) + umac(i+1,j))
!                    vbar = HALF*(vmac(i,j) + vmac(i,j+1))
!                    ugrads = ubar*(sedgex(i+1,j,comp) - sedgex(i,j,comp))/dx(1) + &
!                             vbar*(sedgey(i,j+1,comp) - sedgey(i,j,comp))/dx(2)
!                    if (comp > 1) adv_s(i,j,comp-1) = -ugrads
! ! mass_fractions require is_cons = TRUE
! !                   else; if (mass_fractions) adv_rho(i,j,1) = -ugrads
! !                   endif
!                 enddo
!              enddo
!           end if
!        end do

!     else if (is_vel) then 

!        do j = lo(2), hi(2)
!           do i = lo(1), hi(1)

!              ubar = HALF*(umac(i,j) + umac(i+1,j))
!              vbar = HALF*(vmac(i,j) + vmac(i,j+1))

!              ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
!                       vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

!              ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
!                       vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

!              adv_s(i,j,1) = -ugradu
!              adv_s(i,j,2) = -ugradv
             
!           enddo
!        enddo
!     end if

!   end subroutine mkadv_2d

!   subroutine mkadv_3d(umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
!                        adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

!     use bl_constants_module

!     integer           , intent(in   ) :: lo(:), hi(:), ng
!     real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxz(lo(1)   :,lo(2)   :,lo(3)   :,:) 
!     real (kind = dp_t), intent(in   ) :: dx(:)
!     real (kind = dp_t), intent(in   ) :: dt
!     logical           , intent(in   ) :: is_vel
!     logical           , intent(in   ) :: is_cons(:)
!     real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1)   :,lo(2)   :,lo(3)   :,:)  

!     !     Local variables
!     integer :: i, j, k, comp
!     real (kind = dp_t) ubar,vbar,wbar
!     real (kind = dp_t) :: ugradu,ugradv,ugradw,ugrads
!     real (kind = dp_t) :: divsu

!     if (.not. is_vel) then

!        do comp = 1,size(sold,dim=4)
!           if (is_cons(comp)) then
!              do k = lo(3), hi(3)
!                 do j = lo(2), hi(2)
!                    do i = lo(1), hi(1)
!                       divsu = (fluxx(i+1,j,k,comp)-fluxx(i,j,k,comp))/dx(1) &
!                             + (fluxy(i,j+1,k,comp)-fluxy(i,j,k,comp))/dx(2) &
!                             + (fluxz(i,j,k+1,comp)-fluxz(i,j,k,comp))/dx(3)
!                       if (comp > 1) then; adv_s(i,j,k,comp-1) = -divsu
!                       else; if (sdc_iters >= 0) adv_rho(i,j,k,1) = -divsu
!                       endif
!                    enddo
!                 enddo
!              enddo
!           else 

!              do k = lo(3), hi(3)
!                 do j = lo(2), hi(2)
!                    do i = lo(1), hi(1)
!                       ubar = half*(umac(i,j,k) + umac(i+1,j,k))
!                       vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
!                       wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))
!                       ugrads = ubar*(sedgex(i+1,j,k,comp) - sedgex(i,j,k,comp))/dx(1) + &
!                            vbar*(sedgey(i,j+1,k,comp) - sedgey(i,j,k,comp))/dx(2) + &
!                            wbar*(sedgez(i,j,k+1,comp) - sedgez(i,j,k,comp))/dx(3)
!                       if (comp > 1) adv_s(i,j,k,comp) = -ugrads
!                    enddo
!                 enddo
!              enddo
!           end if
!        enddo

!     else if (is_vel) then

!        do k = lo(3), hi(3)
!           do j = lo(2), hi(2)
!              do i = lo(1), hi(1)
!                 ubar = half*(umac(i,j,k) + umac(i+1,j,k))
!                 vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
!                 wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))

!                 ugradu = ubar*(sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)

!                 ugradv = ubar*(sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)

!                 ugradw = ubar*(sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)


!                 adv_s(i,j,k,1) = -ugradu
!                 adv_s(i,j,k,2) = -ugradv
!                 adv_s(i,j,k,3) = -ugradw
!              enddo
!           enddo
!        enddo

!     end if

!   end subroutine mkadv_3d
! end subroutine mk_adv

! need to finish changing this subroutine
!   subroutine adv_update(mla,sold,force,snew,adv_s,dx,&
!        dt,t,is_vel,the_bc_level,adv_rho)

!     use bl_constants_module
!     use multifab_physbc_module
!     use ml_restriction_module, only: ml_cc_restriction
!     use multifab_fill_ghost_module

!     type(ml_layout)   , intent(in   ) :: mla
!     type(multifab)    , intent(in   ) :: sold(:)
!     type(multifab)    , intent(in   ) :: force(:)
!     type(multifab)    , intent(inout) :: snew(:)
!     type(multifab)    , intent(inout) :: adv_s(:)
!     real(kind = dp_t) , intent(in   ) :: dx(:,:),dt,t
!     logical           , intent(in   ) :: is_vel
!     type(bc_level)    , intent(in   ) :: the_bc_level(:)
!     type(multifab)    , intent(inout), optional :: adv_rho(:)

!     ! local
!     real(kind=dp_t), pointer :: sop(:,:,:,:)
!     real(kind=dp_t), pointer :: snp(:,:,:,:)
!     real(kind=dp_t), pointer :: ap(:,:,:,:)
!     real(kind=dp_t), pointer :: arhop(:,:,:,:)
!     real(kind=dp_t), pointer :: fp(:,:,:,:)

!     integer :: lo(sold(1)%dim),hi(sold(1)%dim)
!     integer :: i,ng,dm,nscal,n,nlevs

!     nlevs = mla%nlevel
!     dm    = mla%dim

!     ng = sold(1)%ng
!     nscal = multifab_ncomp(sold(1))

!     do n=1,nlevs

!        do i = 1, sold(n)%nboxes
!           if ( multifab_remote(sold(n),i) ) cycle
!           sop    => dataptr(sold(n),i)
!           snp    => dataptr(snew(n),i)
!           ap     => dataptr(adv_s(n),i)
!           fp     => dataptr(force(n),i)
!           lo = lwb(get_box(sold(n),i))
!           hi = upb(get_box(sold(n),i))
!           if (mass_fractions .AND. (.not. is_vel) .AND. sdc_iters >= 0) then
!              arhop => dataptr(adv_rho(n),i)
!              select case (dm)
!              case (2)
!                 call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
!                      sepx(:,:,1,:), sepy(:,:,1,:), &
!                      fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
!                      fp(:,:,1,:) , snp(:,:,1,:), ap(:,:,1,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons, arhop(:,:,1,:))
!              case (3)
!                 wmp    => dataptr( umac(n,3),i)
!                 sepz   => dataptr(sedge(n,3),i)
!                 fluxpz => dataptr( flux(n,3),i)
!                 call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
!                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
!                      fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
!                      fp(:,:,:,:) , snp(:,:,:,:), ap(:,:,:,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons,arhop(:,:,:,:))
!              end select
!           else !adv_rho not needed without mass-fractions
!              select case (dm)
!              case (2)
!                 call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
!                      sepx(:,:,1,:), sepy(:,:,1,:), &
!                      fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
!                      fp(:,:,1,:) , snp(:,:,1,:), ap(:,:,1,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
!              case (3)
!                 wmp    => dataptr( umac(n,3),i)
!                 sepz   => dataptr(sedge(n,3),i)
!                 fluxpz => dataptr( flux(n,3),i)
!                 call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
!                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
!                      fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
!                      fp(:,:,:,:) , snp(:,:,:,:), ap(:,:,:,:), &
!                      lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
!              end select
!           end if
!        end do

!        if (.not. is_vel) then

!           call multifab_fill_boundary(snew(n))
!           call multifab_physbc(snew(n),1,dm+1,nscal,the_bc_level(n),dx(n,:),t)

!        else if (is_vel) then

!           call multifab_fill_boundary(snew(n))
!           call multifab_physbc(snew(n),1,1,dm,the_bc_level(n),dx(n,:),t)

!        end if

!     enddo ! end loop over levels

!     do n=nlevs,2,-1

!        if (.not. is_vel) then

!           call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

!           call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
!                                          the_bc_level(n-1),the_bc_level(n), &
!                                          1,dm+1,nscal,dx(n-1:n,:),t)

!        else if (is_vel) then

!           call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

!           call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
!                                          the_bc_level(n-1),the_bc_level(n), &
!                                          1,1,dm,dx(n-1:n,:),t)

!        end if

!     enddo ! end loop over levels

! contains

!   subroutine adv_update_2d(sold,force,snew,&
!        adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

!     use bl_constants_module

!     integer           , intent(in   ) :: lo(:), hi(:), ng
!     real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)  
!     real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)  
!     real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
!     real (kind = dp_t), intent(in   ) :: dx(:)
!     real (kind = dp_t), intent(in   ) :: dt
!     logical           , intent(in   ) :: is_vel
!     real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1):,lo(2):,:)  

!     integer :: i, j, comp
!     real (kind = dp_t) ubar,vbar
!     real (kind = dp_t) ugradu,ugradv,ugrads
!     real (kind = dp_t) :: divsu

!     if (.not. is_vel) then

!        do comp = 1,size(sold,dim=3)
!           if (is_cons(comp)) then
!              do j = lo(2), hi(2)
!                 do i = lo(1), hi(1)
!                    divsu = (fluxx(i+1,j,comp)-fluxx(i,j,comp))/dx(1) &
!                          + (fluxy(i,j+1,comp)-fluxy(i,j,comp))/dx(2)
!                    if (comp > 1) then; adv_s(i,j,comp-1) = -divsu
!                    else
!                       if (mass_fractions .AND.sdc_iters >= 0) adv_rho(i,j,1) = -divsu
!                    endif
!                    snew(i,j,comp) = sold(i,j,comp) - dt*divsu + dt*force(i,j,comp)
!                 enddo
!              enddo
!           else
!              do j = lo(2), hi(2)
!                 do i = lo(1), hi(1)
!                    ubar = HALF*(umac(i,j) + umac(i+1,j))
!                    vbar = HALF*(vmac(i,j) + vmac(i,j+1))
!                    ugrads = ubar*(sedgex(i+1,j,comp) - sedgex(i,j,comp))/dx(1) + &
!                             vbar*(sedgey(i,j+1,comp) - sedgey(i,j,comp))/dx(2)
!                    if (comp > 1) adv_s(i,j,comp-1) = -ugrads
! ! mass_fractions require is_cons = TRUE
! !                   else; if (mass_fractions) adv_rho(i,j,1) = -ugrads
! !                   endif
!                    if (comp > 1) adv_s(i,j,comp-1) = -ugrads
!                    snew(i,j,comp) = sold(i,j,comp) - dt * ugrads + dt * force(i,j,comp)
!                 enddo
!              enddo
!           end if
!        end do

!     else if (is_vel) then 

!        do j = lo(2), hi(2)
!           do i = lo(1), hi(1)

!              ubar = HALF*(umac(i,j) + umac(i+1,j))
!              vbar = HALF*(vmac(i,j) + vmac(i,j+1))

!              ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
!                       vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

!              ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
!                       vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

!              snew(i,j,1) = sold(i,j,1) - dt * ugradu + dt * force(i,j,1)
!              snew(i,j,2) = sold(i,j,2) - dt * ugradv + dt * force(i,j,2)

!              adv_s(i,j,1) = -ugradu
!              adv_s(i,j,2) = -ugradv
             
!           enddo
!        enddo
!     end if

!   end subroutine adv_update_2d

!   subroutine adv_update_3d(sold,umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
!                        force,snew,adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

!     use bl_constants_module

!     integer           , intent(in   ) :: lo(:), hi(:), ng
!     real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
!     real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
!     real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
!     real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)  
!     real (kind = dp_t), intent(in   ) ::   fluxz(lo(1)   :,lo(2)   :,lo(3)   :,:) 
!     real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
!     real (kind = dp_t), intent(in   ) :: dx(:)
!     real (kind = dp_t), intent(in   ) :: dt
!     logical           , intent(in   ) :: is_vel
!     logical           , intent(in   ) :: is_cons(:)
!     real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1)   :,lo(2)   :,lo(3)   :,:)  

!     !     Local variables
!     integer :: i, j, k, comp
!     real (kind = dp_t) ubar,vbar,wbar
!     real (kind = dp_t) :: ugradu,ugradv,ugradw,ugrads
!     real (kind = dp_t) :: divsu

!     if (.not. is_vel) then

!        do comp = 1,size(sold,dim=4)
!           if (is_cons(comp)) then
!              do k = lo(3), hi(3)
!                 do j = lo(2), hi(2)
!                    do i = lo(1), hi(1)
!                       divsu = (fluxx(i+1,j,k,comp)-fluxx(i,j,k,comp))/dx(1) &
!                             + (fluxy(i,j+1,k,comp)-fluxy(i,j,k,comp))/dx(2) &
!                             + (fluxz(i,j,k+1,comp)-fluxz(i,j,k,comp))/dx(3)
!                       if (comp > 1) then; adv_s(i,j,k,comp-1) = -divsu
!                       else; if (sdc_iters >= 0) adv_rho(i,j,k,1) = -divsu
!                       endif
!                       snew(i,j,k,comp) = sold(i,j,k,comp) - dt*divsu + dt*force(i,j,k,comp)
!                    enddo
!                 enddo
!              enddo
!           else 

!              do k = lo(3), hi(3)
!                 do j = lo(2), hi(2)
!                    do i = lo(1), hi(1)
!                       ubar = half*(umac(i,j,k) + umac(i+1,j,k))
!                       vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
!                       wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))
!                       ugrads = ubar*(sedgex(i+1,j,k,comp) - sedgex(i,j,k,comp))/dx(1) + &
!                            vbar*(sedgey(i,j+1,k,comp) - sedgey(i,j,k,comp))/dx(2) + &
!                            wbar*(sedgez(i,j,k+1,comp) - sedgez(i,j,k,comp))/dx(3)
!                       if (comp > 1) adv_s(i,j,k,comp) = -ugrads
!                       snew(i,j,k,comp) = sold(i,j,k,comp) - dt*ugrads + dt*force(i,j,k,comp)
!                    enddo
!                 enddo
!              enddo
!           end if
!        enddo

!     else if (is_vel) then

!        do k = lo(3), hi(3)
!           do j = lo(2), hi(2)
!              do i = lo(1), hi(1)
!                 ubar = half*(umac(i,j,k) + umac(i+1,j,k))
!                 vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
!                 wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))

!                 ugradu = ubar*(sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)

!                 ugradv = ubar*(sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)

!                 ugradw = ubar*(sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1) + &
!                      vbar*(sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2) + &
!                      wbar*(sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)

!                 snew(i,j,k,1) = sold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
!                 snew(i,j,k,2) = sold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
!                 snew(i,j,k,3) = sold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

!                 adv_s(i,j,k,1) = -ugradu
!                 adv_s(i,j,k,2) = -ugradv
!                 adv_s(i,j,k,3) = -ugradw
!              enddo
!           enddo
!        enddo

!     end if

!   end subroutine adv_update_3d
! end subroutine adv_update

  subroutine update(mla,sold,umac,sedge,flux,force,snew,adv_s,dx,&
       dt,t,is_vel,is_cons,the_bc_level,adv_rho)

    use bl_constants_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    type(multifab)    , intent(in   ) :: sedge(:,:)
    type(multifab)    , intent(in   ) :: flux(:,:)
    type(multifab)    , intent(in   ) :: force(:)
    type(multifab)    , intent(inout) :: snew(:)
    type(multifab)    , intent(inout) :: adv_s(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt,t
    logical           , intent(in   ) :: is_vel,is_cons(:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)
    type(multifab)    , intent(inout), optional :: adv_rho(:)

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: ap(:,:,:,:)
    real(kind=dp_t), pointer :: arhop(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    integer :: lo(sold(1)%dim),hi(sold(1)%dim)
    integer :: i,ng,dm,nscal,n,nlevs

    nlevs = mla%nlevel
    dm    = mla%dim

    ng = sold(1)%ng
    nscal = multifab_ncomp(sold(1))

    do n=1,nlevs

       do i = 1, sold(n)%nboxes
          if ( multifab_remote(sold(n),i) ) cycle
          sop    => dataptr(sold(n),i)
          snp    => dataptr(snew(n),i)
          ap     => dataptr(adv_s(n),i)
          ump    => dataptr(umac(n,1),i)
          vmp    => dataptr(umac(n,2),i)
          sepx   => dataptr(sedge(n,1),i)
          sepy   => dataptr(sedge(n,2),i)
          fluxpx => dataptr(flux(n,1),i)
          fluxpy => dataptr(flux(n,2),i)
          fp     => dataptr(force(n),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          if (mass_fractions .AND. (.not. is_vel) .AND. sdc_iters >= 0) then
             arhop => dataptr(adv_rho(n),i)
             select case (dm)
             case (2)
                call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                     sepx(:,:,1,:), sepy(:,:,1,:), &
                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                     fp(:,:,1,:) , snp(:,:,1,:), ap(:,:,1,:), &
                     lo, hi, ng, dx(n,:), dt, is_vel, is_cons, arhop(:,:,1,:))
             case (3)
                wmp    => dataptr( umac(n,3),i)
                sepz   => dataptr(sedge(n,3),i)
                fluxpz => dataptr( flux(n,3),i)
                call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                     sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                     fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                     fp(:,:,:,:) , snp(:,:,:,:), ap(:,:,:,:), &
                     lo, hi, ng, dx(n,:), dt, is_vel, is_cons,arhop(:,:,:,:))
             end select
          else !adv_rho not needed without mass-fractions
             select case (dm)
             case (2)
                call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                     sepx(:,:,1,:), sepy(:,:,1,:), &
                     fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                     fp(:,:,1,:) , snp(:,:,1,:), ap(:,:,1,:), &
                     lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
             case (3)
                wmp    => dataptr( umac(n,3),i)
                sepz   => dataptr(sedge(n,3),i)
                fluxpz => dataptr( flux(n,3),i)
                call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                     sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                     fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                     fp(:,:,:,:) , snp(:,:,:,:), ap(:,:,:,:), &
                     lo, hi, ng, dx(n,:), dt, is_vel, is_cons)
             end select
          end if
       end do

       if (.not. is_vel) then

          call multifab_fill_boundary(snew(n))
          call multifab_physbc(snew(n),1,dm+1,nscal,the_bc_level(n),dx(n,:),t)

       else if (is_vel) then

          call multifab_fill_boundary(snew(n))
          call multifab_physbc(snew(n),1,1,dm,the_bc_level(n),dx(n,:),t)

       end if

    enddo ! end loop over levels

    do n=nlevs,2,-1

       if (.not. is_vel) then

          call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         1,dm+1,nscal,dx(n-1:n,:),t)

       else if (is_vel) then

          call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

          call multifab_fill_ghost_cells(snew(n),snew(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         1,1,dm,dx(n-1:n,:),t)

       end if

    enddo ! end loop over levels

contains

  subroutine update_2d(sold,umac,vmac,sedgex,sedgey,fluxx,fluxy,force,snew,&
       adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
    real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,:) 
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: is_vel
    logical           , intent(in   ) :: is_cons(:)
    real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1):,lo(2):,:)  

    integer :: i, j, comp
    real (kind = dp_t) ubar,vbar
    real (kind = dp_t) ugradu,ugradv,ugrads
    real (kind = dp_t) :: divsu

    if (.not. is_vel) then

       do comp = 1,size(sold,dim=3)
          if (is_cons(comp)) then
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   divsu = (fluxx(i+1,j,comp)-fluxx(i,j,comp))/dx(1) &
                         + (fluxy(i,j+1,comp)-fluxy(i,j,comp))/dx(2)
                   if (comp > 1) then; adv_s(i,j,comp-1) = -divsu
                   else
                      if (mass_fractions .AND.sdc_iters >= 0) adv_rho(i,j,1) = -divsu
                   endif
                   snew(i,j,comp) = sold(i,j,comp) - dt*divsu + dt*force(i,j,comp)
                enddo
             enddo
          else
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   ubar = HALF*(umac(i,j) + umac(i+1,j))
                   vbar = HALF*(vmac(i,j) + vmac(i,j+1))
                   ugrads = ubar*(sedgex(i+1,j,comp) - sedgex(i,j,comp))/dx(1) + &
                            vbar*(sedgey(i,j+1,comp) - sedgey(i,j,comp))/dx(2)
                   if (comp > 1) adv_s(i,j,comp-1) = -ugrads
! mass_fractions require is_cons = TRUE
!                   else; if (mass_fractions) adv_rho(i,j,1) = -ugrads
!                   endif
                   if (comp > 1) adv_s(i,j,comp-1) = -ugrads
                   snew(i,j,comp) = sold(i,j,comp) - dt * ugrads + dt * force(i,j,comp)
                enddo
             enddo
          end if
       end do

    else if (is_vel) then 

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ubar = HALF*(umac(i,j) + umac(i+1,j))
             vbar = HALF*(vmac(i,j) + vmac(i,j+1))

             ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
                      vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

             ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
                      vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

             snew(i,j,1) = sold(i,j,1) - dt * ugradu + dt * force(i,j,1)
             snew(i,j,2) = sold(i,j,2) - dt * ugradv + dt * force(i,j,2)

             adv_s(i,j,1) = -ugradu
             adv_s(i,j,2) = -ugradv
             
          enddo
       enddo
    end if

  end subroutine update_2d

  subroutine update_3d(sold,umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
                       force,snew,adv_s,lo,hi,ng,dx,dt,is_vel,is_cons,adv_rho)

    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(  out) ::    adv_s(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
    real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
    real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
    real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) ::   fluxz(lo(1)   :,lo(2)   :,lo(3)   :,:) 
    real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: is_vel
    logical           , intent(in   ) :: is_cons(:)
    real (kind = dp_t), intent(  out), optional :: adv_rho(lo(1)   :,lo(2)   :,lo(3)   :,:)  

    !     Local variables
    integer :: i, j, k, comp
    real (kind = dp_t) ubar,vbar,wbar
    real (kind = dp_t) :: ugradu,ugradv,ugradw,ugrads
    real (kind = dp_t) :: divsu

    if (.not. is_vel) then

       do comp = 1,size(sold,dim=4)
          if (is_cons(comp)) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      divsu = (fluxx(i+1,j,k,comp)-fluxx(i,j,k,comp))/dx(1) &
                            + (fluxy(i,j+1,k,comp)-fluxy(i,j,k,comp))/dx(2) &
                            + (fluxz(i,j,k+1,comp)-fluxz(i,j,k,comp))/dx(3)
                      if (comp > 1) then; adv_s(i,j,k,comp-1) = -divsu
                      else; if (sdc_iters >= 0) adv_rho(i,j,k,1) = -divsu
                      endif
                      snew(i,j,k,comp) = sold(i,j,k,comp) - dt*divsu + dt*force(i,j,k,comp)
                   enddo
                enddo
             enddo
          else 

             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = lo(1), hi(1)
                      ubar = half*(umac(i,j,k) + umac(i+1,j,k))
                      vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
                      wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))
                      ugrads = ubar*(sedgex(i+1,j,k,comp) - sedgex(i,j,k,comp))/dx(1) + &
                           vbar*(sedgey(i,j+1,k,comp) - sedgey(i,j,k,comp))/dx(2) + &
                           wbar*(sedgez(i,j,k+1,comp) - sedgez(i,j,k,comp))/dx(3)
                      if (comp > 1) adv_s(i,j,k,comp) = -ugrads
                      snew(i,j,k,comp) = sold(i,j,k,comp) - dt*ugrads + dt*force(i,j,k,comp)
                   enddo
                enddo
             enddo
          end if
       enddo

    else if (is_vel) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                ubar = half*(umac(i,j,k) + umac(i+1,j,k))
                vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
                wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))

                ugradu = ubar*(sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1) + &
                     vbar*(sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2) + &
                     wbar*(sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)

                ugradv = ubar*(sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1) + &
                     vbar*(sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2) + &
                     wbar*(sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)

                ugradw = ubar*(sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1) + &
                     vbar*(sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2) + &
                     wbar*(sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)

                snew(i,j,k,1) = sold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
                snew(i,j,k,2) = sold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
                snew(i,j,k,3) = sold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

                adv_s(i,j,k,1) = -ugradu
                adv_s(i,j,k,2) = -ugradv
                adv_s(i,j,k,3) = -ugradw
             enddo
          enddo
       enddo

    end if

  end subroutine update_3d

end subroutine update

end module update_module
