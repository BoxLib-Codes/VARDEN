module advance_module

  use bl_types
  use multifab_module
  use macproject_module
  use viscous_module
  use mkflux_module
  use mkutrans_module
  use mkforce_module
  use hgproject_module
  use setbc_module

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

contains

   subroutine advance(uold,unew,sold,snew,rhohalf,&
                      umac,uedgex,uedgey,uedgez,&
                      sedgex,sedgey,sedgez, &
                      utrans,gp,ext_force,ext_scal_force,dx,time,dt, &
                      phys_bc,norm_vel_bc,tang_vel_bc, &
                      scal_bc,press_bc,visc_coef,diff_coef)

      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: umac
      type(multifab) , intent(inout) :: unew
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: uedgex,uedgey,uedgez,utrans
      type(multifab) , intent(inout) :: sedgex,sedgey,sedgez
      type(multifab) , intent(inout) :: gp
      type(multifab) , intent(inout) :: ext_force
      type(multifab) , intent(inout) :: ext_scal_force
      real(kind=dp_t), intent(inout) :: dx(:),time,dt
      integer        , intent(in   ) ::   phys_bc(:,:)
      integer        , intent(in   ) :: norm_vel_bc(:,:)
      integer        , intent(in   ) :: tang_vel_bc(:,:)
      integer        , intent(in   ) :: scal_bc(:,:)
      integer        , intent(in   ) :: press_bc(:,:)
      real(kind=dp_t), intent(in   ) :: visc_coef
      real(kind=dp_t), intent(in   ) :: diff_coef

      type(multifab) :: force,scal_force

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: unp(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer:: uepx(:,:,:,:)
      real(kind=dp_t), pointer:: uepy(:,:,:,:)
      real(kind=dp_t), pointer:: uepz(:,:,:,:)

      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:)
      real(kind=dp_t), pointer::  up(:,:,:,:)
      real(kind=dp_t), pointer::  sp(:,:,:,:)
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)

      integer :: irz
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: velpred
      integer :: dm,ng_u,ng_rho
      integer :: i
      integer :: comp
      logical :: is_conservative
      logical :: is_vel
      real(kind=dp_t) :: visc_fac,mu

      ng_u   = uold%ng
      ng_rho = rhohalf%ng
      dm      = uold%dim
   
      irz = 0

      call multifab_build(     force,ext_force%la,dm,1)
      call multifab_build(scal_force,ext_force%la, 2,1)

!     Impose boundary conditions on sold and uold.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold , i)
         sop => dataptr(sold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call setvelbc_2d (uop(:,:,1,:), lo, ng_u, phys_bc, visc_coef, irz)
              call setscalbc_2d(sop(:,:,1,:), lo, ng_u, phys_bc)
            case (3)
              call setvelbc_3d (uop(:,:,:,:), lo, ng_u, phys_bc, visc_coef)
              call setscalbc_3d(sop(:,:,:,:), lo, ng_u, phys_bc)
         end select
      end do
      call multifab_fill_boundary(uold)
      call multifab_fill_boundary(sold)
  
!     Create force.
      visc_fac = ONE
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(force, i)
          ep => dataptr(ext_force, i)
         gpp => dataptr(gp   , i)
          rp => dataptr(sold , i)
          up => dataptr(uold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkforce_2d(fp(:,:,1,:), ep(:,:,1,:), gpp(:,:,1,:), rp(:,:,1,1), up(:,:,1,:), &
                              ng_u, ng_u, &
                              dx, norm_vel_bc, tang_vel_bc, visc_coef, visc_fac)
            case (3)
              call mkforce_3d(fp(:,:,:,:), ep(:,:,:,:), gpp(:,:,:,:), rp(:,:,:,1), up(:,:,:,:), &
                              ng_u, ng_u, &
                              dx, norm_vel_bc, tang_vel_bc, visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(force)

!     Create utrans.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold  , i)
         unp => dataptr(unew  , i)
         utp => dataptr(utrans, i)
          fp => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkutrans_2d(uop(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                               lo,dx,dt,ng_u,phys_bc,irz)
            case (3)
              call mkutrans_3d(uop(:,:,:,:), utp(:,:,:,:), fp(:,:,:,:), &
                               lo,dx,dt,ng_u,phys_bc)
         end select
      end do

!     Create the edge states to be used for the MAC velocity 
      velpred = 1
      is_vel = .true.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold  , i)
         uepx => dataptr(uedgex, i)
         uepy => dataptr(uedgey, i)
         ump  => dataptr(umac  , i)
         utp  => dataptr(utrans, i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             ump(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, irz, phys_bc, velpred, ng_u)
            case (3)
              uepz => dataptr(uedgez, i)
              call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                             uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                             ump(:,:,:,:), utp(:,:,:,:), fp(:,:,:,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, phys_bc, velpred, ng_u)
         end select
      end do

      call macproject(umac,sold,dx,press_bc)
  
!     Create scalar force at time n.
      visc_fac = ONE
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(scal_force, i)
          ep => dataptr(ext_scal_force, i)
          sp => dataptr(sold , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sp(:,:,1,:), &
                                  ng_u, dx, scal_bc, diff_coef, visc_fac)
            case (3)
              call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sp(:,:,:,:), &
                                  ng_u, dx, scal_bc, diff_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force)

!     Create the edge states of scalar using the MAC velocity 
      velpred = 0
      is_vel = .false.
      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         uop  => dataptr(uold, i)
         sepx => dataptr(sedgex, i)
         sepy => dataptr(sedgey, i)
         ump  => dataptr(umac, i)
         utp  => dataptr(utrans, i)
          fp  => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             ump(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, irz, phys_bc, velpred, ng_u)
            case (3)
              sepz => dataptr(sedgez, i)
              call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                             sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                             ump(:,:,:,:), utp(:,:,:,:), fp(:,:,:,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, phys_bc, velpred, ng_u)
         end select
      end do
  
!     Create scalar force at time n+1/2.
      visc_fac = HALF
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(scal_force, i)
          ep => dataptr(ext_scal_force, i)
          sp => dataptr(sold , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call mkscalforce_2d(fp(:,:,1,:), ep(:,:,1,:), sp(:,:,1,:), &
                                  ng_u, dx, scal_bc, diff_coef, visc_fac)
            case (3)
              call mkscalforce_3d(fp(:,:,:,:), ep(:,:,:,:), sp(:,:,:,:), &
                                  ng_u, dx, scal_bc, diff_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(scal_force)
  
!     Update the scalar with conservative differencing
      is_conservative = .true.
      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         sop => dataptr(sold, i)
         ump => dataptr(umac, i)
         sepx => dataptr(sedgex, i)
         sepy => dataptr(sedgey, i)
         snp => dataptr(snew, i)
          rp => dataptr(rhohalf, i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_2d(sop(:,:,1,:), ump(:,:,1,:), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             fp(:,:,1,:) , snp(:,:,1,:), &
                             rp(:,:,1,1) , &
                             lo, hi, ng_u, dx, time, dt, is_conservative)
              call setscalbc_2d(rp(:,:,1,:), lo, ng_rho, phys_bc)
            case (3)
              uepz => dataptr(uedgez, i)
              call update_3d(sop(:,:,:,:) , ump(:,:,:,:),  &
                             sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                             fp(:,:,:,:) , snp(:,:,:,:), &
                             rp(:,:,:,1) , &
                             lo, hi, ng_u, dx, time, dt, is_conservative)
              call setscalbc_3d(rp(:,:,:,:), lo, ng_rho, phys_bc)
         end select
      end do
      call multifab_fill_boundary(rhohalf)

!     Do the diffusive solve for Crank-Nicolson discretization of tracer.
      if (diff_coef > ZERO) then
        comp = 2
        mu = HALF*dt*diff_coef
        call diff_scalar_solve(snew,dx,mu,scal_bc,comp)
      end if

!     Create the edge states of velocity using the MAC velocity 
      velpred = 0
      is_vel = .true.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop  => dataptr(uold, i)
         uepx => dataptr(uedgex, i)
         uepy => dataptr(uedgey, i)
         ump  => dataptr(umac, i)
         utp  => dataptr(utrans, i)
          fp  => dataptr(force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(uop(:,:,1,:), uop(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             ump(:,:,1,:), utp(:,:,1,:), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, irz, phys_bc, velpred, ng_u)
            case (3)
              call mkflux_3d(uop(:,:,:,:), uop(:,:,:,:), &
                             uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                             ump(:,:,:,:), utp(:,:,:,:), fp(:,:,:,:), &
                             lo, dx, dt, is_vel, &
                             visc_coef, phys_bc, velpred, ng_u)
         end select
      end do
  
!     Create the force at half-time using rhohalf and half the viscous term.
      visc_fac = HALF
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(force, i)
          ep => dataptr(ext_force, i)
         gpp => dataptr(gp   , i)
          rp => dataptr(rhohalf , i)
          up => dataptr(uold , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkforce_2d(fp(:,:,1,:), ep(:,:,1,:), gpp(:,:,1,:), rp(:,:,1,1), up(:,:,1,:), &
                              ng_u, ng_rho, &
                              dx, norm_vel_bc, tang_vel_bc, visc_coef, visc_fac)
            case (3)
              call mkforce_3d(fp(:,:,:,:), ep(:,:,:,:), gpp(:,:,:,:), rp(:,:,:,1), up(:,:,:,:), &
                              ng_u, ng_rho, &
                              dx, norm_vel_bc, tang_vel_bc, visc_coef, visc_fac)
         end select
      end do
      call multifab_fill_boundary(force)

!     Update the velocity with convective differencing
      is_conservative = .false.
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         uop => dataptr(uold, i)
         ump => dataptr(umac, i)
         uepx => dataptr(uedgex, i)
         uepy => dataptr(uedgey, i)
         unp => dataptr(unew, i)
          fp => dataptr(force, i)
          rp => dataptr(rhohalf, i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_2d(uop(:,:,1,:), ump(:,:,1,:), &
                             uepx(:,:,1,:), uepy(:,:,1,:), &
                             fp(:,:,1,:), unp(:,:,1,:), &
                             rp(:,:,1,1), &
                             lo, hi, ng_u, dx, time, dt, is_conservative)
            case (3)
              uepz => dataptr(uedgez, i)
              call update_3d(uop(:,:,:,:), ump(:,:,:,:),  &
                             uepx(:,:,:,:), uepy(:,:,:,:), uepz(:,:,:,:), &
                             fp(:,:,:,:), unp(:,:,:,:), &
                             rp(:,:,:,1), &
                             lo, hi, ng_u, dx, time, dt, is_conservative)
         end select
      end do

!     Do the viscous solve for Crank-Nicolson discretization.
      if (visc_coef > ZERO) then
        mu = HALF*dt*visc_coef
        call visc_solve(unew,rhohalf,dx,mu,norm_vel_bc)
      end if

!     Project the new velocity field.
      call hgproject(unew,rhohalf,gp,dx,dt,press_bc)

      call multifab_destroy(force)

   end subroutine advance

   subroutine update_2d (sold,umac,sedgex,sedgey,force,snew,rhohalf, &
                         lo,hi,ng,dx,time,dt,is_cons)

      implicit none

      integer, intent(in) :: lo(2), hi(2), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) :: dx(2)
      real (kind = dp_t), intent(in   ) :: time,dt
      logical :: is_cons

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv
      real (kind = dp_t) :: divsu

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
      end if

   end subroutine update_2d

   subroutine update_3d (sold,umac,sedgex,sedgey,sedgez,force,snew,rhohalf, &
                         lo,hi,ng,dx,time,dt,is_cons)

      implicit none

      integer, intent(in) :: lo(3), hi(3), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(inout) :: rhohalf(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) :: dx(3)
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

end module advance_module
