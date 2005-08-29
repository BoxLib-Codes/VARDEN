module update_module

  use bl_types
  use multifab_module

  implicit none

  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  contains

   subroutine update_2d (sold,umac,vmac,sedgex,sedgey,force,snew,rhohalf, &
                         lo,hi,ng,dx,time,dt,is_vel,is_cons)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt
      logical           , intent(in   ) :: is_vel
      logical           , intent(in   ) :: is_cons(:)

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax

      if (.not. is_vel) then

        do n = 1,size(sold,dim=3)
         if (is_cons(n)) then
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             divsu = (umac(i+1,j) * sedgex(i+1,j,n) &
                     -umac(i  ,j) * sedgex(i  ,j,n) ) / dx(1) + &
                     (vmac(i,j+1) * sedgey(i,j+1,n) &
                     -vmac(i,j  ) * sedgey(i,j  ,n) ) / dx(2)

             snew(i,j,n) = sold(i,j,n) - dt * divsu + dt * force(i,j,n)
             if (n.eq.1) rhohalf(i,j) = HALF * (sold(i,j,1) + snew(i,j,1))
           enddo
           enddo
           smax = snew(lo(1),lo(2),n) 
           smin = snew(lo(1),lo(2),n) 
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             smax = max(smax,snew(i,j,n))
             smin = min(smin,snew(i,j,n))
           enddo
           enddo
           print *,'MIN/MAX OF RHO ',smin,smax
        else
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             ubar = HALF*(umac(i,j) + umac(i+1,j))
             vbar = HALF*(vmac(i,j) + vmac(i,j+1))
             ugrads = ubar*(sedgex(i+1,j,n) - sedgex(i,j,n))/dx(1) + &
                      vbar*(sedgey(i,j+1,n) - sedgey(i,j,n))/dx(2)
             snew(i,j,n) = sold(i,j,n) - dt * ugrads + dt * force(i,j,n)
           enddo
           enddo
           smax = snew(lo(1),lo(2),n) 
           smin = snew(lo(1),lo(2),n) 
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             smax = max(smax,snew(i,j,n))
             smin = min(smin,snew(i,j,n))
           enddo
           enddo
           print *,'MIN/MAX OF TEMP ',smin,smax
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
         umax = snew(lo(1),lo(2),1) 
         umin = snew(lo(1),lo(2),1) 
         vmax = snew(lo(1),lo(2),2) 
         vmin = snew(lo(1),lo(2),2) 
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             umax = max(umax,snew(i,j,1))
             umin = min(umin,snew(i,j,1))
             vmax = max(vmax,snew(i,j,2))
             vmin = min(vmin,snew(i,j,2))
           enddo
         enddo
         print *,'MIN/MAX OF UNEW ',umin,umax
         print *,'MIN/MAX OF VNEW ',vmin,vmax

    end if

   end subroutine update_2d

   subroutine update_3d (sold,umac,vmac,wmac,sedgex,sedgey,sedgez,force,snew,rhohalf, &
                         lo,hi,ng,dx,time,dt,is_vel,is_cons)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt
      logical           , intent(in   ) :: is_vel
      logical           , intent(in   ) :: is_cons(:)

!     Local variables
      integer :: i, j, k, n
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) :: ugradu,ugradv,ugradw,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax,wmin,wmax

      if (.not. is_vel) then

        do n = 1,size(sold,dim=4)
         if (is_cons(n)) then
           do k = lo(3), hi(3)
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)

             divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,n) &
                     -umac(i  ,j,k) * sedgex(i  ,j,k,n) ) / dx(1) &
                    +(vmac(i,j+1,k) * sedgey(i,j+1,k,n) &
                     -vmac(i,j  ,k) * sedgey(i,j  ,k,n) ) / dx(2) &
                    +(wmac(i,j,k+1) * sedgez(i,j,k+1,n) &
                     -wmac(i,j,k  ) * sedgez(i,j,k  ,n) ) / dx(3)

             snew(i,j,k,n) = sold(i,j,k,n) - dt * divsu + dt * force(i,j,k,n)
             if (n.eq.1)  rhohalf(i,j,k) = HALF * (sold(i,j,k,1) + snew(i,j,k,1))
           enddo
           enddo
           enddo
           smin = snew(lo(1),lo(2),lo(3),n) 
           smax = snew(lo(1),lo(2),lo(3),n) 
           do k = lo(3), hi(3)
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             smax = max(smax,snew(i,j,k,n))
             smin = min(smin,snew(i,j,k,n))
           enddo
           enddo
           enddo
           print *,'MIN/MAX OF RHO ',smin,smax
         else 

           do k = lo(3), hi(3)
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             ubar = half*(umac(i,j,k) + umac(i+1,j,k))
             vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
             wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))
             ugrads = ubar*(sedgex(i+1,j,k,n) - sedgex(i,j,k,n))/dx(1) + &
                      vbar*(sedgey(i,j+1,k,n) - sedgey(i,j,k,n))/dx(2) + &
                      wbar*(sedgez(i,j,k+1,n) - sedgez(i,j,k,n))/dx(3)
             snew(i,j,k,n) = sold(i,j,k,n) - dt * ugrads + dt * force(i,j,k,n)
           enddo
           enddo
           enddo
           smin = snew(lo(1),lo(2),lo(3),n) 
           smax = snew(lo(1),lo(2),lo(3),n) 
           do k = lo(3), hi(3)
           do j = lo(2), hi(2)
           do i = lo(1), hi(1)
             smax = max(smax,snew(i,j,k,n))
             smin = min(smin,snew(i,j,k,n))
           enddo
           enddo
           enddo
           print *,'MIN/MAX OF TEMP ',smin,smax
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
         umax = snew(lo(1),lo(2),lo(3),1) 
         umin = snew(lo(1),lo(2),lo(3),1) 
         vmax = snew(lo(1),lo(2),lo(3),2) 
         vmin = snew(lo(1),lo(2),lo(3),2) 
         wmax = snew(lo(1),lo(2),lo(3),3) 
         wmin = snew(lo(1),lo(2),lo(3),3) 
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
         do i = lo(1), hi(1)
             umax = max(umax,snew(i,j,k,1))
             umin = min(umin,snew(i,j,k,1))
             vmax = max(umax,snew(i,j,k,2))
             vmin = min(umin,snew(i,j,k,2))
             wmax = max(wmax,snew(i,j,k,3))
             wmin = min(wmin,snew(i,j,k,3))
         enddo
         enddo
         enddo
         print *,'MIN/MAX OF UNEW ',umin,umax
         print *,'MIN/MAX OF VNEW ',vmin,vmax
         print *,'MIN/MAX OF WNEW ',wmin,wmax

      end if

   end subroutine update_3d

end module update_module
