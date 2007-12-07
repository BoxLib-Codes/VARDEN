module update_module

  use bl_types
  use multifab_module
  use multifab_physbc_module
  use define_bc_module

  implicit none

  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  private
  public :: update

  contains

   subroutine update(sold,umac,sedge,flux,force,snew,rhohalf, &
                     dx,dt,is_vel,is_cons,the_bc_level)

    implicit none

    type(multifab)    , intent(in   ) :: sold
    type(multifab)    , intent(in   ) :: umac(:)
    type(multifab)    , intent(in   ) :: sedge(:)
    type(multifab)    , intent(in   ) :: flux(:)
    type(multifab)    , intent(in   ) :: force
    type(multifab)    , intent(inout) :: snew
    type(multifab)    , intent(inout) :: rhohalf
    real(kind = dp_t) , intent(in   ) :: dx(:),dt
    logical           , intent(in   ) :: is_vel, is_cons(:)
    type(bc_level)    , intent(in   ) :: the_bc_level
    
    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    integer :: lo(sold%dim),hi(sold%dim)
    integer :: i,ng,dm,bc_comp,nscal

    dm = sold%dim
    ng = sold%ng

    do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         snp  => dataptr(snew, i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         fluxpx => dataptr(flux(1), i)
         fluxpy => dataptr(flux(2), i)
          rp  => dataptr(rhohalf, i)
          fp  => dataptr(force, i)
         lo = lwb(get_box(sold, i))
         hi = upb(get_box(sold, i))
         select case (dm)
         case (2)
            call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                           sepx(:,:,1,:), sepy(:,:,1,:), &
                           fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                           fp(:,:,1,:) , snp(:,:,1,:), &
                           rp(:,:,1,1) , &
                           lo, hi, ng, dx,dt,is_vel,is_cons)
         case (3)
            wmp    => dataptr( umac(3), i)
            sepz   => dataptr(sedge(3), i)
            fluxpz => dataptr( flux(3), i)
            call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                           sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                           fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                           fp(:,:,:,:) , snp(:,:,:,:), &
                           rp(:,:,:,1) , &
                           lo, hi, ng, dx,dt,is_vel,is_cons)
         end select
      end do

      if (.not. is_vel) then

         call multifab_fill_boundary(snew)
         call multifab_fill_boundary(rhohalf)

         bc_comp = dm+1
         nscal = multifab_ncomp(snew)
         call multifab_physbc(snew   ,1,bc_comp,nscal,dx,the_bc_level)
         call multifab_physbc(rhohalf,1,bc_comp,    1,dx,the_bc_level)

      else if (is_vel) then

         call multifab_fill_boundary(snew)
         call multifab_physbc(snew,1,1,dm,dx,the_bc_level)

      end if

   end subroutine update

   subroutine update_2d (sold,umac,vmac,sedgex,sedgey,fluxx,fluxy,force,snew,rhohalf, &
                         lo,hi,ng,dx,dt,is_vel,is_cons)

      implicit none

      integer           , intent(in   ) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   fluxx(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   fluxy(lo(1)   :,lo(2)   :,:) 
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt
      logical           , intent(in   ) :: is_vel
      logical           , intent(in   ) :: is_cons(:)

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
             snew(i,j,comp) = sold(i,j,comp) - dt * divsu + dt * force(i,j,comp)
             if (comp.eq.1) rhohalf(i,j) = HALF * (sold(i,j,1) + snew(i,j,1))
           enddo
           enddo
        else
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             ubar = HALF*(umac(i,j) + umac(i+1,j))
             vbar = HALF*(vmac(i,j) + vmac(i,j+1))
             ugrads = ubar*(sedgex(i+1,j,comp) - sedgex(i,j,comp))/dx(1) + &
                      vbar*(sedgey(i,j+1,comp) - sedgey(i,j,comp))/dx(2)
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

           enddo
         enddo
    end if

   end subroutine update_2d

   subroutine update_3d (sold,umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
                         force,snew,rhohalf,lo,hi,ng,dx,dt,is_vel,is_cons)

      implicit none

      integer           , intent(in   ) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
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
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt
      logical           , intent(in   ) :: is_vel
      logical           , intent(in   ) :: is_cons(:)

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
             snew(i,j,k,comp) = sold(i,j,k,comp) - dt * divsu + dt * force(i,j,k,comp)
             if (comp.eq.1)  rhohalf(i,j,k) = HALF * (sold(i,j,k,1) + snew(i,j,k,1))
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
             snew(i,j,k,comp) = sold(i,j,k,comp) - dt * ugrads + dt * force(i,j,k,comp)
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
         enddo
         enddo
         enddo

      end if

   end subroutine update_3d

end module update_module
