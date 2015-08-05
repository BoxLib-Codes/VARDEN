module hgproject_module

  use bl_types
  use bl_error_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: hgproject

contains 

  subroutine hgproject(proj_type,mla,unew,uold,rhohalf,p,gp,dx,dt,the_bc_tower, &
                       press_comp,divu_rhs)

    use bl_constants_module
    use bc_module
    use proj_parameters
    use nodal_stencil_module
    use ml_restrict_fill_module
    use hg_multigrid_module        , only : hg_multigrid
    use hg_hypre_module            , only : hg_hypre
    use probin_module              , only : verbose, use_hypre
    use stencil_types_module

    integer        , intent(in   ) :: proj_type
    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: uold(:)
    type(multifab ), intent(inout) :: rhohalf(:)
    type(multifab ), intent(inout) :: p(:)
    type(multifab ), intent(inout) :: gp(:)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp

    type(multifab ), intent(inout), optional :: divu_rhs(:)

    ! Local  
    type(multifab) :: phi(mla%nlevel),gphi(mla%nlevel),rh(mla%nlevel)
    logical        :: nodal(mla%dim)
    integer        :: n,nlevs,dm,ng,stencil_type
    real(dp_t)     :: umin,umax,vmin,vmax,wmin,wmax

    real(dp_t)      :: rel_solver_eps
    real(dp_t)      :: abs_solver_eps

    stencil_type = ND_DENSE_STENCIL
!    stencil_type = ND_CROSS_STENCIL

    nlevs = mla%nlevel
    dm = mla%dim
    ng = nghost(unew(nlevs))

    nodal = .true.

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print *,'PROJ_TYPE IN HGPROJECT:',proj_type
    endif

    do n = 1, nlevs
       call multifab_build(  rh(n),mla%la(n),1,1,nodal)
       call multifab_build( phi(n), mla%la(n), 1, 1, nodal)
       call multifab_build(gphi(n), mla%la(n), dm, 0) 
       call setval( rh(n),ZERO,all=.true.)
       call setval(phi(n),ZERO,all=.true.)
    end do

    if (verbose .ge. 1) then
       umin = 1.d30
       vmin = 1.d30
       wmin = 1.d30
       umax = -1.d30
       vmax = -1.d30
       wmax = -1.d30
       do n = 1, nlevs
          umin = min(umin,multifab_min_c(unew(n),1))
          umax = max(umax,multifab_max_c(unew(n),1))
          vmin = min(vmin,multifab_min_c(unew(n),2))
          vmax = max(vmax,multifab_max_c(unew(n),2))
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,1001) umin,umax
          write(6,1002) vmin,vmax
          if (dm .eq. 3) write(6,1003) wmin,wmax
          write(6,1004)
       end if
    end if

1001 format('... x-velocity before projection ',e17.10,2x,e17.10)
1002 format('... y-velocity before projection ',e17.10,2x,e17.10)
1003 format('... z-velocity before projection ',e17.10,2x,e17.10)
1004 format(' ')

    call create_uvec_for_projection(nlevs,unew,uold,rhohalf,gp,dt,the_bc_tower,proj_type)

    if (nlevs .eq. 1) then
       rel_solver_eps = 1.d-12
    else if (nlevs .eq. 2) then
       rel_solver_eps = 1.d-11
    else
       rel_solver_eps = 1.d-10
    endif

    abs_solver_eps = -1.d0 

    if (use_hypre .eq. 1) then
       call hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                     press_comp,stencil_type, &
                     rel_solver_eps,abs_solver_eps,divu_rhs)
    else
       call hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                         press_comp,stencil_type, &
                         rel_solver_eps,abs_solver_eps,divu_rhs)
    end if

    call mkgphi(nlevs,gphi,phi,dx)

    call hg_update(mla,proj_type,unew,uold,gp,gphi,rhohalf,  &
                   p,phi,dt,the_bc_tower%bc_tower_array)

    if (verbose .ge. 1) then
       umin = 1.d30
       vmin = 1.d30
       wmin = 1.d30
       umax = -1.d30
       vmax = -1.d30
       wmax = -1.d30
       do n = 1, nlevs
          umin = min(umin,multifab_min_c(unew(n),1))
          umax = max(umax,multifab_max_c(unew(n),1))
          vmin = min(vmin,multifab_min_c(unew(n),2))
          vmax = max(vmax,multifab_max_c(unew(n),2))
          if (dm .eq. 3) then
             wmin = min(wmin,multifab_min_c(unew(n),3))
             wmax = max(wmax,multifab_max_c(unew(n),3))
          end if
       end do
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,1101) umin,umax
          write(6,1102) vmin,vmax
          if (dm .eq. 3) write(6,1103) wmin,wmax
          write(6,1104)
       end if
    end if

1101 format('... x-velocity  after projection ',e17.10,2x,e17.10)
1102 format('... y-velocity  after projection ',e17.10,2x,e17.10)
1103 format('... z-velocity  after projection ',e17.10,2x,e17.10)
1104 format(' ')

    do n = 1,nlevs
       call multifab_destroy(rh(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(gphi(n))
    end do

  contains

    subroutine create_uvec_for_projection(nlevs,unew,uold,rhohalf,gp,dt,the_bc_tower,proj_type)

      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: gp(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_tower) , intent(in   ) :: the_bc_tower
      integer        , intent(in   ) :: proj_type
  
      type(bc_level) :: bc
  
      real(kind=dp_t), pointer :: unp(:,:,:,:)
      real(kind=dp_t), pointer :: uop(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:)
  
      integer :: i,n,dm,ng_u,ng_g,ng_r
      integer :: lo(get_dim(uold(1))),hi(get_dim(uold(1)))
  
      type(bl_prof_timer), save :: bpt

      call build(bpt,"create_uvec_for_projection")

      dm = get_dim(unew(nlevs))
      ng_u = nghost(unew(nlevs))
      ng_g = nghost(gp(nlevs))
      ng_r = nghost(rhohalf(nlevs))
  
      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nfabs(unew(n))
            unp => dataptr(unew(n)   , i) 
            uop => dataptr(uold(n)   , i) 
            gpp => dataptr(gp(n)     , i)
             rp => dataptr(rhohalf(n), i)
             lo =  lwb(get_box(unew(n),i))
             hi =  upb(get_box(unew(n),i))
            select case (dm)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), uop(:,:,1,:), ng_u, rp(:,:,1,1), ng_r, &
                                     gpp(:,:,1,:), ng_g, dt, &
                                     bc%phys_bc_level_array(i,:,:), proj_type, lo, hi)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), uop(:,:,:,:), ng_u, rp(:,:,:,1), ng_r, &
                                     gpp(:,:,:,:), ng_g, dt, &
                                     bc%phys_bc_level_array(i,:,:), proj_type, lo, hi)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do 

      call destroy(bpt)

    end subroutine create_uvec_for_projection


    !   ******************************************************************************* !

    subroutine mkgphi(nlevs,gphi,phi,dx)

      integer       , intent(in   ) :: nlevs
      type(multifab), intent(inout) :: gphi(:)
      type(multifab), intent(in   ) :: phi(:)
      real(dp_t) :: dx(:,:)

      integer :: i,dm,n,lo(get_dim(phi(1))),hi(get_dim(phi(1)))
      integer :: ng_g, ng_p

      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt,"mkgphi")

      dm = get_dim(phi(1))

      ng_g = gphi(1)%ng
      ng_p = phi(1)%ng

      do n = 1, nlevs

         do i = 1, nfabs(phi(n))
            gph => dataptr(gphi(n),i)
            pp  => dataptr(phi(n),i)
            lo =  lwb(get_box(gphi(n),i))
            hi =  upb(get_box(gphi(n),i))
            select case (dm)
            case (2)
               call mkgphi_2d(lo, hi, gph(:,:,1,:), ng_g, pp(:,:,1,1), ng_p, dx(n,:))
            case (3)
               call mkgphi_3d(lo, hi, gph(:,:,:,:), ng_g, pp(:,:,:,1), ng_p, dx(n,:))
            end select
         end do

         call multifab_fill_boundary(gphi(n))

      end do

      call destroy(bpt)

    end subroutine mkgphi

    !   ********************************************************************************** !

    subroutine hg_update(mla,proj_type,unew,uold,gp,gphi,rhohalf,p,phi,dt,the_bc_level)

      type(ml_layout), intent(in   ) :: mla
      integer        , intent(in   ) :: proj_type
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(inout) :: gp(:)
      type(multifab) , intent(in   ) :: gphi(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: p(:)
      type(multifab) , intent(in   ) :: phi(:)
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_level) , intent(in   ) :: the_bc_level(:)

      ! local
      integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
      integer :: ng_u,ng_g,ng_h,ng_r,ng_p,ng_i

      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: uon(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 
      real(kind=dp_t), pointer ::  ph(:,:,:,:) 
      real(kind=dp_t), pointer ::  pp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt,"hg_update")

      nlevs = mla%nlevel
      dm    = mla%dim

      ng_u = unew(1)%ng
      ng_g = gp(1)%ng
      ng_h = gphi(1)%ng
      ng_r = rhohalf(1)%ng
      ng_p = p(1)%ng
      ng_i = phi(1)%ng

      do n = 1, nlevs

         do i = 1, nfabs(unew(n))
            upn => dataptr(unew(n),i)
            uon => dataptr(uold(n),i)
            gpp => dataptr(gp(n),i)
            gph => dataptr(gphi(n),i)
            rp  => dataptr(rhohalf(n),i)
            pp  => dataptr(p(n),i)
            ph  => dataptr(phi(n),i)
            lo =  lwb(get_box(phi(n),i))
            hi =  upb(get_box(phi(n),i))
            select case (dm)
            case (2)
               call hg_update_2d(proj_type, upn(:,:,1,:), uon(:,:,1,:), ng_u, &
                                 gpp(:,:,1,:), ng_g, gph(:,:,1,:), ng_h, rp(:,:,1,1), ng_r, &
                                 pp(:,:,1,1), ng_p, ph(:,:,1,1), ng_i, dt, lo, hi)
            case (3)
               call hg_update_3d(proj_type, upn(:,:,:,:), uon(:,:,:,:), ng_u, &
                                 gpp(:,:,:,:), ng_g, gph(:,:,:,:), ng_h, rp(:,:,:,1), ng_r, &
                                 pp(:,:,:,1), ng_p, ph(:,:,:,1), ng_i, dt, lo, hi)
            end select
         end do

      end do

      do n = nlevs, 2, -1
         call ml_cc_restriction(  gp(n-1),  gp(n),mla%mba%rr(n-1,:))
      end do

      do n = 1, nlevs
         call multifab_fill_boundary(gp(n))
         call multifab_fill_boundary(p(n))
      end do

      if (nlevs > 1) then
         call ml_restrict_and_fill(nlevs, unew, mla%mba%rr, the_bc_level)
      end if

      call destroy(bpt)

    end subroutine hg_update

    !   ******************************************************************************** !

    subroutine create_uvec_2d(unew,uold,ng_u,rhohalf,ng_r,gp,ng_g,dt,phys_bc,proj_type,lo,hi)

      use proj_parameters

      integer        , intent(in   ) :: ng_u,ng_r,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_u:,lo(2)-ng_u:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_u:,lo(2)-ng_u:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_r:,lo(2)-ng_r:)
      real(kind=dp_t), intent(inout) ::      gp(lo(1)-ng_g:,lo(2)-ng_g:,:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      real(kind=dp_t) :: dtinv

      dtinv = 1.d0 / dt

      if (phys_bc(1,1) .eq. INLET) gp(lo(1)-1,lo(2)-1:hi(2)+1,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(hi(1)+1,lo(2)-1:hi(2)+1,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(lo(1)-1:hi(1)+1,lo(2)-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(lo(1)-1:hi(1)+1,hi(2)+1,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) = ( &
             unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) - uold(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) ) * dtinv
         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) = ( &
             unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) - uold(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) ) * dtinv
     
      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) = unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1) + &
            dt*gp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)/rhohalf(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
         unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) = unew(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2) + &
            dt*gp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)/rhohalf(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

       else
     
          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(lo(1)-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(hi(1)+1,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,lo(2)-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,hi(2)+1,:) = ZERO


    end subroutine create_uvec_2d

!   *************************************************************************************** !

    subroutine create_uvec_3d(unew,uold,ng_u,rhohalf,ng_r,gp,ng_g,dt,phys_bc,proj_type,lo,hi)

      use proj_parameters

      integer        , intent(in   ) :: ng_u,ng_r,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
      real(kind=dp_t), intent(inout) ::      gp(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      ! local
      integer :: i,j,k,m
      real(kind=dp_t) :: dtinv

      dtinv = 1.d0 / dt

      if (phys_bc(1,1) .eq. INLET) gp(lo(1)-1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(lo(1)-1:hi(1)+1,lo(2)-1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(lo(1)-1:hi(1)+1,hi(2)+1,lo(3)-1:hi(3)+1,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,hi(3)+1,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  do i=lo(1)-1,hi(1)+1
                     unew(i,j,k,m) = (unew(i,j,k,m) - uold(i,j,k,m)) * dtinv
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL

      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         !$OMP PARALLEL PRIVATE(i,j,k,m)
         do m=1,3
            !$OMP DO
            do k=lo(3)-1,hi(3)+1
               do j=lo(2)-1,hi(2)+1
                  do i=lo(1)-1,hi(1)+1
                     unew(i,j,k,m) = unew(i,j,k,m) + dt*gp(i,j,k,m)/rhohalf(i,j,k)
                  end do
               end do
            end do
            !$OMP END DO NOWAIT
         end do
         !$OMP END PARALLEL

      else

          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(lo(1)-1,:,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(hi(1)+1,:,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,lo(2)-1,:,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,hi(2)+1,:,:) = ZERO
      if (phys_bc(3,1)==SLIP_WALL .or. phys_bc(3,1)==NO_SLIP_WALL) unew(:,:,lo(3)-1,:) = ZERO
      if (phys_bc(3,2)==SLIP_WALL .or. phys_bc(3,2)==NO_SLIP_WALL) unew(:,:,hi(3)+1,:) = ZERO

    end subroutine create_uvec_3d

    !   ********************************************************************************* !

    subroutine mkgphi_2d(lo,hi,gp,ng_g,phi,ng_p,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_g,ng_p
      real(kind=dp_t), intent(inout) ::  gp(lo(1)-ng_g:,lo(2)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j
      real(kind=dp_t) :: dxinv(2)

      dxinv(1) = 1.d0/dx(1)
      dxinv(2) = 1.d0/dx(2)

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            gp(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                              phi(i  ,j) - phi(i  ,j+1) ) * dxinv(1)
            gp(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                              phi(i,j  ) - phi(i+1,j  ) ) * dxinv(2)
         end do
      end do

    end subroutine mkgphi_2d

    !   ******************************************************************************** !

    subroutine mkgphi_3d(lo,hi,gp,ng_g,phi,ng_p,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_g,ng_p
      real(kind=dp_t), intent(inout) ::  gp(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k
      real(kind=dp_t) :: dxinv(3)

      dxinv(1) = 1.d0/dx(1)
      dxinv(2) = 1.d0/dx(2)
      dxinv(3) = 1.d0/dx(3)


      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               gp(i,j,k,1) = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                    -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) * dxinv(1)
               gp(i,j,k,2) = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) * dxinv(2)
               gp(i,j,k,3) = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) * dxinv(3)
            end do
         end do
      end do

    end subroutine mkgphi_3d

    !   ****************************************************************************** !

    subroutine hg_update_2d(proj_type,unew,uold,ng_u,gp,ng_g,gphi,ng_h, &
                            rhohalf,ng_r,p,ng_p,phi,ng_i,dt,lo,hi)

      use proj_parameters

      integer        , intent(in   ) :: proj_type,lo(:),hi(:)
      integer        , intent(in   ) :: ng_u,ng_g,ng_h,ng_r,ng_p,ng_i
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_u:,lo(2)-ng_u:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_u:,lo(2)-ng_u:,:)
      real(kind=dp_t), intent(inout) ::      gp(lo(1)-ng_g:,lo(2)-ng_g:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(lo(1)-ng_h:,lo(2)-ng_h:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_r:,lo(2)-ng_r:)
      real(kind=dp_t), intent(inout) ::       p(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) ::     phi(lo(1)-ng_i:,lo(2)-ng_i:)
      real(kind=dp_t), intent(in   ) :: dt

      real(kind=dp_t) :: dtinv

      dtinv = 1.d0 / dt

      !     Subtract off the density-weighted gradient.
      unew(lo(1):hi(1),lo(2):hi(2),1) = &
           unew(lo(1):hi(1),lo(2):hi(2),1) - gphi(lo(1):hi(1),lo(2):hi(2),1)/rhohalf(lo(1):hi(1),lo(2):hi(2)) 
      unew(lo(1):hi(1),lo(2):hi(2),2) = &
           unew(lo(1):hi(1),lo(2):hi(2),2) - gphi(lo(1):hi(1),lo(2):hi(2),2)/rhohalf(lo(1):hi(1),lo(2):hi(2)) 

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(lo(1):hi(1),lo(2):hi(2),:) = &
           uold(lo(1):hi(1),lo(2):hi(2),:) + dt * unew(lo(1):hi(1),lo(2):hi(2),:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(lo(1):hi(1)  ,lo(2):hi(2),:) = gp(lo(1):hi(1)  ,lo(2):hi(2),:) + &
                                         gphi(lo(1):hi(1),lo(2):hi(2),:)
          p(lo(1):hi(1)+1,lo(2):hi(2)+1) =  p(lo(1):hi(1)+1,lo(2):hi(2)+1) + &
                                          phi(lo(1):hi(1)+1,lo(2):hi(2)+1)

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(lo(1):hi(1)  ,lo(2):hi(2),:) = dtinv * gphi(lo(1):hi(1)  ,lo(2):hi(2),:)
          p(lo(1):hi(1)+1,lo(2):hi(2)+1) = dtinv *  phi(lo(1):hi(1)+1,lo(2):hi(2)+1)

      end if

    end subroutine hg_update_2d

    !   ******************************************************************************* !

    subroutine hg_update_3d(proj_type,unew,uold,ng_u,gp,ng_g,gphi,ng_h, &
                            rhohalf,ng_r,p,ng_p,phi,ng_i,dt,lo,hi)

      use proj_parameters

      integer        , intent(in   ) :: proj_type,lo(:),hi(:)
      integer        , intent(in   ) :: ng_u,ng_g,ng_h,ng_r,ng_p,ng_i
      real(kind=dp_t), intent(inout) ::    unew(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
      real(kind=dp_t), intent(in   ) ::    uold(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
      real(kind=dp_t), intent(inout) ::      gp(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
      real(kind=dp_t), intent(inout) ::       p(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) ::     phi(lo(1)-ng_i:,lo(2)-ng_i:,lo(3)-ng_i:)
      real(kind=dp_t), intent(in   ) :: dt

      real(kind=dp_t) :: dtinv

      dtinv = 1.d0 / dt

      !     Subtract off the density-weighted gradient.
      unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = &
           unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) &
           - gphi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)/rhohalf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) 
      unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) = &
           unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2) &
           - gphi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2)/rhohalf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) 
      unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) = &
           unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3) &
           - gphi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3)/rhohalf(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) 

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = &
           uold(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) + dt * unew(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = gp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) &
              + gphi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)
          p(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1) =  p(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1) &
                                                      + phi(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3),:) = dtinv * &
             gphi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)
          p(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1) = dtinv * &
             phi(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)

      end if

    end subroutine hg_update_3d

  end subroutine hgproject

end module hgproject_module
