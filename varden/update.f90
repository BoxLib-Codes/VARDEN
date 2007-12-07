module update_module

  use bl_types
  use multifab_module

  implicit none

  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  contains

   subroutine update_2d (lev,sold,umac,vmac,sedgex,sedgey,fluxx,fluxy,force,snew,rhohalf, &
                         lo,hi,ng,dx,dt,is_vel,is_cons)

      implicit none

      integer           , intent(in   ) :: lev, lo(:), hi(:), ng
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

   subroutine update_3d (lev,sold,umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
                         force,snew,rhohalf,lo,hi,ng,dx,dt,is_vel,is_cons)

      implicit none

      integer           , intent(in   ) :: lev, lo(:), hi(:), ng
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
