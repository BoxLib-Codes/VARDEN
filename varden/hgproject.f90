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
                       press_comp,divu_rhs,div_coeff_1d,div_coeff_3d)

    use bl_constants_module
    use bc_module
    use proj_parameters
    use nodal_stencil_module
    use multifab_fill_ghost_module , only : multifab_fill_ghost_cells
    use ml_restriction_module      , only : ml_cc_restriction
    use hg_multigrid_module        , only : hg_multigrid
    use hg_hypre_module            , only : hg_hypre
    use probin_module              , only : verbose, use_hypre

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
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_3d(:)

    ! Local  
    type(multifab) :: phi(mla%nlevel),gphi(mla%nlevel),rh(mla%nlevel)
    logical        :: nodal(mla%dim),use_div_coeff_1d, use_div_coeff_3d
    integer        :: n,nlevs,dm,ng,stencil_type
    real(dp_t)     :: umin,umax,vmin,vmax,wmin,wmax

    real(dp_t)      :: rel_solver_eps
    real(dp_t)      :: abs_solver_eps

    ! stencil_type = ST_DENSE
    stencil_type = ST_CROSS

    nlevs = mla%nlevel
    dm = mla%dim
    ng = nghost(unew(nlevs))

    nodal = .true.

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       print *,'PROJ_TYPE IN HGPROJECT:',proj_type
    endif

    use_div_coeff_1d = .false.
    if (present(div_coeff_1d)) use_div_coeff_1d = .true.

    use_div_coeff_3d = .false.
    if (present(div_coeff_3d)) use_div_coeff_3d = .true.

    if (use_div_coeff_1d .and. use_div_coeff_3d) &
       call bl_error('CANT HAVE 1D and 3D DIV_COEFF IN HGPROJECT ')

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

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(nlevs,unew,div_coeff_1d,.true.)
       call mult_by_1d_coeff(nlevs,rhohalf,div_coeff_1d,.false.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(nlevs,unew,div_coeff_3d,.true.)
       call mult_by_3d_coeff(nlevs,rhohalf,div_coeff_3d,.false.)
    end if

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

    if (use_div_coeff_1d) then
       call mult_by_1d_coeff(nlevs,unew,div_coeff_1d,.false.)
       call mult_by_1d_coeff(nlevs,rhohalf,div_coeff_1d,.true.)
    else if (use_div_coeff_3d) then
       call mult_by_3d_coeff(nlevs,unew,div_coeff_3d,.false.)
       call mult_by_3d_coeff(nlevs,rhohalf,div_coeff_3d,.true.)
    end if

    call mkgphi(nlevs,gphi,phi,dx)

    call hg_update(mla,proj_type,unew,uold,gp,gphi,rhohalf,  &
                   p,phi,ng,dt,the_bc_tower%bc_tower_array)

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
  
      integer :: i,n,dm,ng
  
      dm = get_dim(unew(nlevs))
      ng = nghost(unew(nlevs))
  
      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, nboxes(unew(n))
            if ( multifab_remote(unew(n), i) ) cycle
            unp => dataptr(unew(n)     , i) 
            uop => dataptr(uold(n)     , i) 
            gpp => dataptr(gp(n)       , i)
             rp => dataptr(  rhohalf(n), i)
            select case (dm)
               case (2)
                 call create_uvec_2d(unp(:,:,1,:), uop(:,:,1,:), rp(:,:,1,1), gpp(:,:,1,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng, proj_type)
               case (3)
                 call create_uvec_3d(unp(:,:,:,:), uop(:,:,:,:), rp(:,:,:,1), gpp(:,:,:,:), dt, &
                                     bc%phys_bc_level_array(i,:,:), ng, proj_type)
            end select
         end do
         call multifab_fill_boundary(unew(n))
      end do 

    end subroutine create_uvec_for_projection


    !   ******************************************************************************* !

    subroutine mkgphi(nlevs,gphi,phi,dx)

      integer       , intent(in   ) :: nlevs
      type(multifab), intent(inout) :: gphi(:)
      type(multifab), intent(in   ) :: phi(:)
      real(dp_t) :: dx(:,:)

      integer :: i,dm,n,lo(get_dim(phi(1))),hi(get_dim(phi(1)))

      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: pp(:,:,:,:) 

      dm = get_dim(phi(1))

      do n = 1, nlevs

         do i = 1, nboxes(phi(n))
            if ( multifab_remote(phi(n),i) ) cycle
            gph => dataptr(gphi(n),i)
            pp  => dataptr(phi(n),i)
            lo =  lwb(get_box(gphi(n),i))
            hi =  upb(get_box(gphi(n),i))
            select case (dm)
            case (2)
               call mkgphi_2d(lo, hi, gph(:,:,1,:), pp(:,:,1,1), dx(n,:))
            case (3)
               call mkgphi_3d(lo, hi, gph(:,:,:,:), pp(:,:,:,1), dx(n,:))
            end select
         end do

         call multifab_fill_boundary(gphi(n))

      end do

    end subroutine mkgphi

    !   ********************************************************************************** !

    subroutine hg_update(mla,proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt,the_bc_level)

      type(ml_layout), intent(in   ) :: mla
      integer        , intent(in   ) :: proj_type
      type(multifab) , intent(inout) :: unew(:)
      type(multifab) , intent(in   ) :: uold(:)
      type(multifab) , intent(inout) :: gp(:)
      type(multifab) , intent(in   ) :: gphi(:)
      type(multifab) , intent(in   ) :: rhohalf(:)
      type(multifab) , intent(inout) :: p(:)
      type(multifab) , intent(in   ) :: phi(:)
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(in   ) :: dt
      type(bc_level) , intent(in   ) :: the_bc_level(:)

      ! local
      integer :: i,dm,n,nlevs

      real(kind=dp_t), pointer :: upn(:,:,:,:) 
      real(kind=dp_t), pointer :: uon(:,:,:,:) 
      real(kind=dp_t), pointer :: gph(:,:,:,:) 
      real(kind=dp_t), pointer :: gpp(:,:,:,:) 
      real(kind=dp_t), pointer ::  rp(:,:,:,:) 
      real(kind=dp_t), pointer ::  ph(:,:,:,:) 
      real(kind=dp_t), pointer ::  pp(:,:,:,:) 

      nlevs = mla%nlevel
      dm    = mla%dim

      do n = 1, nlevs

         do i = 1, nboxes(unew(n))
            if ( multifab_remote(unew(n),i) ) cycle
            upn => dataptr(unew(n),i)
            uon => dataptr(uold(n),i)
            gpp => dataptr(gp(n),i)
            gph => dataptr(gphi(n),i)
            rp  => dataptr(rhohalf(n),i)
            pp  => dataptr(p(n),i)
            ph  => dataptr(phi(n),i)
            select case (dm)
            case (2)
               call hg_update_2d(proj_type, upn(:,:,1,:), uon(:,:,1,:), gpp(:,:,1,:), &
                                 gph(:,:,1,:),rp(:,:,1,1),pp(:,:,1,1), ph(:,:,1,1), ng, dt)
            case (3)
               call hg_update_3d(proj_type, upn(:,:,:,:), uon(:,:,:,:), gpp(:,:,:,:), &
                                 gph(:,:,:,:),rp(:,:,:,1),pp(:,:,:,1), ph(:,:,:,1), ng, dt)
            end select
         end do

         call multifab_fill_boundary(unew(n))
         call multifab_fill_boundary(gp(n))
         call multifab_fill_boundary(p(n))

      end do

      do n = nlevs, 2, -1
         call ml_cc_restriction(unew(n-1),unew(n),mla%mba%rr(n-1,:)) 
         call ml_cc_restriction(  gp(n-1),  gp(n),mla%mba%rr(n-1,:))
      end do

      do n = 2, nlevs
         call multifab_fill_ghost_cells(unew(n),unew(n-1),ng,mla%mba%rr(n-1,:), &
                                        the_bc_level(n-1),the_bc_level(n),1,1,dm)
      end do

    end subroutine hg_update

    !   ******************************************************************************** !

    subroutine create_uvec_2d(unew,uold,rhohalf,gp,dt,phys_bc,ng,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny
      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         unew(-1:nx,-1:ny,1) = ( unew(-1:nx,-1:ny,1) - uold(-1:nx,-1:ny,1) ) / dt
         unew(-1:nx,-1:ny,2) = ( unew(-1:nx,-1:ny,2) - uold(-1:nx,-1:ny,2) ) / dt
     
      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         unew(-1:nx,-1:ny,1) = unew(-1:nx,-1:ny,1) + dt*gp(-1:nx,-1:ny,1)/rhohalf(-1:nx,-1:ny)
         unew(-1:nx,-1:ny,2) = unew(-1:nx,-1:ny,2) + dt*gp(-1:nx,-1:ny,2)/rhohalf(-1:nx,-1:ny)

       else
     
          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(-1,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(nx,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,-1,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,ny,:) = ZERO

    end subroutine create_uvec_2d

!   *************************************************************************************** !

    subroutine create_uvec_3d(unew,uold,rhohalf,gp,dt,phys_bc,ng,proj_type)

      use proj_parameters

      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:, -1:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt
      integer        , intent(in   ) :: phys_bc(:,:)
      integer        , intent(in   ) :: proj_type

      integer :: nx,ny,nz

      nx = size(gp,dim=1) - 2
      ny = size(gp,dim=2) - 2
      nz = size(gp,dim=3) - 2

      if (phys_bc(1,1) .eq. INLET) gp(-1,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(1,2) .eq. INLET) gp(nx,-1:ny,-1:nz,:) = ZERO
      if (phys_bc(2,1) .eq. INLET) gp(-1:nx,-1,-1:nz,:) = ZERO
      if (phys_bc(2,2) .eq. INLET) gp(-1:nx,ny,-1:nz,:) = ZERO
      if (phys_bc(3,1) .eq. INLET) gp(-1:nx,-1:ny,-1,:) = ZERO
      if (phys_bc(3,2) .eq. INLET) gp(-1:nx,-1:ny,nz,:) = ZERO

      ! quantity projected is U
      if (proj_type .eq. initial_projection) then

      ! quantity projected is U
      else if (proj_type .eq. divu_iters) then

      ! quantity projected is (Ustar - Un)
      else if (proj_type .eq. pressure_iters) then

         unew(-1:nx,-1:ny,-1:nz,:) = ( unew(-1:nx,-1:ny,-1:nz,:) - uold(-1:nx,-1:ny,-1:nz,:) ) / dt

      ! quantity projected is Ustar + dt * (1/rho) Gp
      else if (proj_type .eq. regular_timestep) then

         unew(-1:nx,-1:ny,-1:nz,1) = unew(-1:nx,-1:ny,-1:nz,1) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,1)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,2) = unew(-1:nx,-1:ny,-1:nz,2) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,2)/rhohalf(-1:nx,-1:ny,-1:nz)
         unew(-1:nx,-1:ny,-1:nz,3) = unew(-1:nx,-1:ny,-1:nz,3) + &
                                    dt*gp(-1:nx,-1:ny,-1:nz,3)/rhohalf(-1:nx,-1:ny,-1:nz)

      else

          call bl_error('No proj_type by this number ')

      end if

      if (phys_bc(1,1)==SLIP_WALL .or. phys_bc(1,1)==NO_SLIP_WALL) unew(-1,:,:,:) = ZERO
      if (phys_bc(1,2)==SLIP_WALL .or. phys_bc(1,2)==NO_SLIP_WALL) unew(nx,:,:,:) = ZERO
      if (phys_bc(2,1)==SLIP_WALL .or. phys_bc(2,1)==NO_SLIP_WALL) unew(:,-1,:,:) = ZERO
      if (phys_bc(2,2)==SLIP_WALL .or. phys_bc(2,2)==NO_SLIP_WALL) unew(:,ny,:,:) = ZERO
      if (phys_bc(3,1)==SLIP_WALL .or. phys_bc(3,1)==NO_SLIP_WALL) unew(:,:,-1,:) = ZERO
      if (phys_bc(3,2)==SLIP_WALL .or. phys_bc(3,2)==NO_SLIP_WALL) unew(:,:,nz,:) = ZERO

    end subroutine create_uvec_3d

    !   ********************************************************************************* !

    subroutine mkgphi_2d(lo,hi,gp,phi,dx)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  gp(lo(1)  :,lo(2)  :,:)
      real(kind=dp_t), intent(inout) :: phi(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            gp(i,j,1) = HALF*(phi(i+1,j) + phi(i+1,j+1) - &
                              phi(i  ,j) - phi(i  ,j+1) ) /dx(1)
            gp(i,j,2) = HALF*(phi(i,j+1) + phi(i+1,j+1) - &
                              phi(i,j  ) - phi(i+1,j  ) ) /dx(2)
         end do
      end do

    end subroutine mkgphi_2d

    !   ******************************************************************************** !

    subroutine mkgphi_3d(lo,hi,gp,phi,dx)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  gp(lo(1):  ,lo(2):  ,lo(3):  ,1:)
      real(kind=dp_t), intent(inout) :: phi(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               gp(i,j,k,1) = FOURTH*(phi(i+1,j,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i+1,j,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i  ,j,k  ) - phi(i  ,j+1,k  ) &
                    -phi(i  ,j,k+1) - phi(i  ,j+1,k+1) ) /dx(1)
               gp(i,j,k,2) = FOURTH*(phi(i,j+1,k  ) + phi(i+1,j+1,k  ) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j  ,k+1) - phi(i+1,j  ,k+1) ) /dx(2)
               gp(i,j,k,3) = FOURTH*(phi(i,j  ,k+1) + phi(i+1,j  ,k+1) &
                    +phi(i,j+1,k+1) + phi(i+1,j+1,k+1) & 
                    -phi(i,j  ,k  ) - phi(i+1,j  ,k  ) &
                    -phi(i,j+1,k  ) - phi(i+1,j+1,k  ) ) /dx(3)
            end do
         end do
      end do

    end subroutine mkgphi_3d

    !   ****************************************************************************** !

    subroutine hg_update_2d(proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt)

      use proj_parameters

      integer        , intent(in   ) :: proj_type
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp( -1:, -1:,:)
      real(kind=dp_t), intent(in   ) ::    gphi(  0:,  0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf( -1:, -1:)
      real(kind=dp_t), intent(inout) ::       p( -1:, -1:)
      real(kind=dp_t), intent(in   ) ::     phi( -1:, -1:)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: nx,ny

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1

      !     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,1) = unew(0:nx,0:ny,1) - gphi(0:nx,0:ny,1)/rhohalf(0:nx,0:ny) 
      unew(0:nx,0:ny,2) = unew(0:nx,0:ny,2) - gphi(0:nx,0:ny,2)/rhohalf(0:nx,0:ny) 

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,:) = uold(0:nx,0:ny,:) + dt * unew(0:nx,0:ny,:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(0:nx,0:ny,:) = gp(0:nx,0:ny,:) +            gphi(0:nx,0:ny,:)
         p(0:nx,0:ny  ) =  p(0:nx,0:ny  ) +             phi(0:nx,0:ny  )

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(0:nx,0:ny,:) = (ONE/dt) * gphi(0:nx,0:ny,:)
         p(0:nx,0:ny  ) = (ONE/dt) *  phi(0:nx,0:ny  )

      end if

    end subroutine hg_update_2d

    !   ******************************************************************************* !

    subroutine hg_update_3d(proj_type,unew,uold,gp,gphi,rhohalf,p,phi,ng,dt)

      use proj_parameters

      integer        , intent(in   ) :: proj_type
      integer        , intent(in   ) :: ng
      real(kind=dp_t), intent(inout) ::    unew(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(in   ) ::    uold(-ng:,-ng:,-ng:,:)
      real(kind=dp_t), intent(inout) ::      gp(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) ::    gphi( 0:, 0:, 0:,:)
      real(kind=dp_t), intent(in   ) :: rhohalf(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) ::       p( -1:, -1:,-1:)
      real(kind=dp_t), intent(in   ) ::     phi( -1:, -1:,-1:)
      real(kind=dp_t), intent(in   ) :: dt

      integer         :: nx,ny,nz

      nx = size(gphi,dim=1)-1
      ny = size(gphi,dim=2)-1
      nz = size(gphi,dim=3)-1

      !     Subtract off the density-weighted gradient.
      unew(0:nx,0:ny,0:nz,1) = &
           unew(0:nx,0:ny,0:nz,1) - gphi(0:nx,0:ny,0:nz,1)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,2) = &
           unew(0:nx,0:ny,0:nz,2) - gphi(0:nx,0:ny,0:nz,2)/rhohalf(0:nx,0:ny,0:nz) 
      unew(0:nx,0:ny,0:nz,3) = &
           unew(0:nx,0:ny,0:nz,3) - gphi(0:nx,0:ny,0:nz,3)/rhohalf(0:nx,0:ny,0:nz) 

      if (proj_type .eq. pressure_iters) &    ! unew held the projection of (ustar-uold)
           unew(0:nx,0:ny,0:nz,:) = uold(0:nx,0:ny,0:nz,:) + dt * unew(0:nx,0:ny,0:nz,:)

      if ( (proj_type .eq. initial_projection) .or. (proj_type .eq. divu_iters) ) then

         gp = ZERO
         p = ZERO

      else if (proj_type .eq. pressure_iters) then

         !  phi held                 (change in pressure)
         ! gphi held the gradient of (change in pressure)
         gp(0:nx,0:ny,0:nz,:) = gp(0:nx,0:ny,0:nz,:) +            gphi(0:nx,0:ny,0:nz,:)
         p(0:nx,0:ny,0:nz  ) =  p(0:nx,0:ny,0:nz  ) +             phi(0:nx,0:ny,0:nz  )

      else if (proj_type .eq. regular_timestep) then

         !  phi held                 dt * (pressure)
         ! gphi held the gradient of dt * (pressure)
         gp(0:nx,0:ny,0:nz,:) = (ONE/dt) * gphi(0:nx,0:ny,0:nz,:)
         p(0:nx,0:ny,0:nz  ) = (ONE/dt) *  phi(0:nx,0:ny,0:nz  )

      end if

    end subroutine hg_update_3d

  end subroutine hgproject

  !   ********************************************************************************* !

  subroutine mult_by_1d_coeff(nlevs,u,div_coeff,do_mult)

    integer       , intent(in   )           :: nlevs
    type(multifab), intent(inout)           :: u(:)
    real(dp_t)    , intent(in   )           :: div_coeff(:,:)
    logical       , intent(in   ), optional :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    integer :: i,ng,n,dm
    integer :: lo(get_dim(u(1))),hi(get_dim(u(1)))
    logical :: local_do_mult

    local_do_mult = .true.
    if (present(do_mult)) local_do_mult = do_mult

    ng = nghost(u(1))
    dm = get_dim(u(1))

    do n = 1, nlevs

       ! Multiply u by div coeff
       do i = 1, nboxes(u(n))
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call mult_by_1d_coeff_2d(ump(:,:,1,:),div_coeff(n,:),lo,hi,ng,local_do_mult)
          case (3)
             call mult_by_1d_coeff_3d(ump(:,:,:,:),div_coeff(n,:),lo,hi,ng,local_do_mult)
          end select
       end do

    end do

  end subroutine mult_by_1d_coeff

  subroutine mult_by_1d_coeff_2d(u,div_coeff,lo,hi,ng,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng:,lo(2)-ng:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(0:)
    logical        , intent(in   ) :: do_mult

    integer :: j

    if (do_mult) then
       do j = lo(2),hi(2)
          u(:,j,:) = u(:,j,:) * div_coeff(j)
       end do
    else
       do j = lo(2),hi(2)
          u(:,j,:) = u(:,j,:) / div_coeff(j)
       end do
    end if

  end subroutine mult_by_1d_coeff_2d

  subroutine mult_by_1d_coeff_3d(u,div_coeff,lo,hi,ng,do_mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(0:)
    logical        , intent(in   ) :: do_mult

    integer :: k

    if (do_mult) then
       do k = lo(3),hi(3)
          u(:,:,k,:) = u(:,:,k,:) * div_coeff(k)
       end do
    else
       do k = lo(3),hi(3)
          u(:,:,k,:) = u(:,:,k,:) / div_coeff(k)
       end do
    end if

  end subroutine mult_by_1d_coeff_3d

  !   *********************************************************************************** !

  subroutine mult_by_3d_coeff(nlevs,u,div_coeff,do_mult)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    type(multifab) , intent(in   ) :: div_coeff(:)
    logical        , intent(in   ) :: do_mult

    ! local
    real(kind=dp_t), pointer :: ump(:,:,:,:) 
    real(kind=dp_t), pointer ::  dp(:,:,:,:) 
    integer :: i,ngu,ngd,n,dm

    ngu = nghost(u(1))
    ngd = nghost(div_coeff(1))
    dm = get_dim(u(1))

    do n = 1, nlevs
       ! Multiply u by div coeff
       do i = 1, nboxes(u(n))
          if ( multifab_remote(u(n),i) ) cycle
          ump => dataptr(u(n),i)
          dp => dataptr(div_coeff(n),i)
          select case (dm)
          case (3)
             call mult_by_3d_coeff_3d(ump(:,:,:,:), ngu, dp(:,:,:,1), ngd, do_mult)
          end select
       end do

    end do

  end subroutine mult_by_3d_coeff

  subroutine mult_by_3d_coeff_3d(u,ngu,div_coeff,ngd,do_mult)

    integer        , intent(in   ) :: ngu,ngd
    real(kind=dp_t), intent(inout) ::         u(-ngu:,-ngu:,-ngu:,:)
    real(dp_t)     , intent(in   ) :: div_coeff(-ngd:,-ngd:,-ngd:)
    logical        , intent(in   ) :: do_mult

    integer :: i,j,k,nx,ny,nz
    nx = size(u,dim=1) - 2*ngu
    ny = size(u,dim=2) - 2*ngu
    nz = size(u,dim=3) - 2*ngu

    if (do_mult) then
       do k = 0,nz-1 
          do j = 0,ny-1 
             do i = 0,nx-1 
                u(i,j,k,:) = u(i,j,k,:) * div_coeff(i,j,k)
             end do
          end do
       end do
    else
       do k = 0,nz-1 
          do j = 0,ny-1 
             do i = 0,nx-1 
                u(i,j,k,:) = u(i,j,k,:) / div_coeff(i,j,k)
             end do
          end do
       end do
    end if

  end subroutine mult_by_3d_coeff_3d

end module hgproject_module
