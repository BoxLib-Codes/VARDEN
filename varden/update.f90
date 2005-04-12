module update_module

  use bl_types
  use multifab_module

  implicit none

  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  contains

   subroutine update_2d (sold,umac,sedgex,sedgey,force,snew,rhohalf, &
                         lo,hi,ng_cell,ng_edge,dx,time,dt,is_cons)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng_cell,ng_edge
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_cell:,lo(2)-ng_cell:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_cell:,lo(2)-ng_cell:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)-ng_edge:,lo(2)-ng_edge:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt
      logical :: is_cons

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: u_symm, v_symm, u_symm_max, v_symm_max

      if (is_cons) then
         do n = 1,size(sold,dim=3)
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)

             divsu = (umac(i+1,j,1) * sedgex(i+1,j,n) &
                     -umac(i  ,j,1) * sedgex(i  ,j,n) ) / dx(1) + &
                     (umac(i,j+1,2) * sedgey(i,j+1,n) &
                     -umac(i,j  ,2) * sedgey(i,j  ,n) ) / dx(2)

             snew(i,j,n) = sold(i,j,n) - dt * divsu + dt * force(i,j,n)

           enddo
         enddo
         enddo

         do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rhohalf(i,j) = HALF * (sold(i,j,1) + snew(i,j,1))
           enddo
         enddo

      else
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)

             ubar = HALF*(umac(i,j,1) + umac(i+1,j,1))
             vbar = HALF*(umac(i,j,2) + umac(i,j+1,2))

             ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
                      vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

             ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
                      vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

             snew(i,j,1) = sold(i,j,1) - dt * ugradu + dt * force(i,j,1)
             snew(i,j,2) = sold(i,j,2) - dt * ugradv + dt * force(i,j,2)

           enddo
         enddo
         u_symm_max = -1.e20
         v_symm_max = -1.e20
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)/2
             u_symm = snew(i,j,1) + snew(hi(1)-i,j,1)
             u_symm_max = max(u_symm_max, abs(u_symm))
             v_symm = snew(i,j,2) - snew(hi(1)-i,j,2)
             v_symm_max = max(v_symm_max, abs(v_symm))
           end do
         end do
         print *,'USYMM MAX ',u_symm_max
         print *,'VSYMM MAX ',v_symm_max
      end if

   end subroutine update_2d

   subroutine update_3d (sold,umac,sedgex,sedgey,sedgez,force,snew,rhohalf, &
                         lo,hi,ng_cell,ng_edge,dx,time,dt,is_cons)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng_cell, ng_edge
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng_cell:,lo(2)-ng_cell:,lo(3)-ng_cell:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng_cell:,lo(2)-ng_cell:,lo(3)-ng_cell:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)-ng_edge:,lo(2)-ng_edge:,lo(3)-ng_edge:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent( in) :: time,dt
      logical :: is_cons

!     Local variables
      integer :: i, j, k, n
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) :: ugradu,ugradv,ugradw
      real (kind = dp_t) :: divsu

      if (is_cons) then
         do n = 1, size(sold,dim=4)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)

             divsu = (umac(i+1,j,k,1) * sedgex(i+1,j,k,n) &
                     -umac(i  ,j,k,1) * sedgex(i  ,j,k,n) ) / dx(1) &
                    +(umac(i,j+1,k,2) * sedgey(i,j+1,k,n) &
                     -umac(i,j  ,k,2) * sedgey(i,j  ,k,n) ) / dx(2) &
                    +(umac(i,j,k+1,3) * sedgez(i,j,k+1,n) &
                     -umac(i,j,k  ,3) * sedgez(i,j,k  ,n) ) / dx(3)

             snew(i,j,k,n) = sold(i,j,k,n) - dt * divsu + dt * force(i,j,k,n)

           enddo
         enddo
         enddo
         enddo
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              rhohalf(i,j,k) = HALF * (sold(i,j,k,1) + snew(i,j,k,1))
           enddo
         enddo
         enddo
      else 
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
           do i = lo(1), hi(1)

             ubar = half*(umac(i,j,k,1) + umac(i+1,j,k,1))
             vbar = half*(umac(i,j,k,2) + umac(i,j+1,k,2))
             wbar = half*(umac(i,j,k,3) + umac(i,j,k+1,3))

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
