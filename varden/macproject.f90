module macproject_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use bl_error_module
  use sparse_solve_module
  use create_umac_grown_module

  implicit none

  private

  public :: macproject

contains 

  subroutine macproject(mla,umac,rho,dx,the_bc_tower,bc_comp,&
                        divu_rhs,div_coeff_1d,div_coeff_half_1d,div_coeff_3d)

    use probin_module, only: stencil_order

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: umac(:,:)
    type(multifab ), intent(inout) :: rho(:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp

    type(multifab ), intent(inout), optional :: divu_rhs(:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_1d(:,:)
    real(dp_t)     , intent(in   ), optional :: div_coeff_half_1d(:,:)
    type(multifab ), intent(in   ), optional :: div_coeff_3d(:)

    ! Local  
    type(multifab)  :: rh(mla%nlevel),phi(mla%nlevel)
    type(multifab)  :: alpha(mla%nlevel),beta(mla%nlevel,mla%dim)
    type(bndry_reg) :: fine_flx(2:mla%nlevel)
    real(dp_t)      :: umac_norm(mla%nlevel)
    integer         :: d,dm,i,n, nlevs
    logical         :: use_rhs, use_div_coeff_1d, use_div_coeff_3d

    nlevs = mla%nlevel
    dm = umac(nlevs,1)%dim

    use_rhs = .false.
    if (present(divu_rhs)) use_rhs = .true.

    use_div_coeff_1d = .false.
    if (present(div_coeff_1d)) use_div_coeff_1d = .true.

    use_div_coeff_3d = .false.
    if (present(div_coeff_3d)) use_div_coeff_3d = .true.

    if (use_div_coeff_1d .and. use_div_coeff_3d) &
       call bl_error('CANT HAVE 1D and 3D DIV_COEFF IN MACPROJECT ')

    do n = 1, nlevs
       call multifab_build(   rh(n), mla%la(n),  1, 0)
       call multifab_build(  phi(n), mla%la(n),  1, 1)
       call multifab_build(alpha(n), mla%la(n),  1, 1)
       do d = 1,dm
          call multifab_build( beta(n,d), mla%la(n), 1, 1, nodal = edge_nodal_flag(d,:))
       end do

       call setval(alpha(n),ZERO,all=.true.)
       call setval(  phi(n),ZERO,all=.true.)

    end do

    if (use_div_coeff_1d) then
       call mult_umac_by_1d_coeff(mla,nlevs,umac,div_coeff_1d,div_coeff_half_1d,.true.)
    else if (use_div_coeff_3d) then
       call mult_umac_by_3d_coeff(mla,nlevs,umac,div_coeff_3d,.true.)
    end if

    ! Compute umac_norm to be used inside the MG solver as part of a stopping criterion
    umac_norm = -1.0_dp_t
    do n = 1,nlevs
       do i = 1,dm
          umac_norm(n) = max(umac_norm(n),norm_inf(umac(n,i)))
       end do
    end do

    if (use_rhs) then
       call divumac(nlevs,umac,rh,dx,mla%mba%rr,.true.,divu_rhs)
    else
       call divumac(nlevs,umac,rh,dx,mla%mba%rr,.true.)
    end if

    call mk_mac_coeffs(nlevs,mla,rho,beta,the_bc_tower)

    if (use_div_coeff_1d) then
       call mult_beta_by_1d_coeff(nlevs,beta,div_coeff_1d,div_coeff_half_1d)
    else if (use_div_coeff_3d) then
       call mult_beta_by_3d_coeff(mla,nlevs,beta,div_coeff_3d)
    end if

    do n = 2,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                       stencil_order,mla%mba%rr,umac_norm)

    call mkumac(rh,umac,phi,beta,fine_flx,dx,the_bc_tower,bc_comp,mla%mba%rr)

    if (use_rhs) then
       call divumac(nlevs,umac,rh,dx,mla%mba%rr,.false.,divu_rhs)
    else
       call divumac(nlevs,umac,rh,dx,mla%mba%rr,.false.)
    end if

    if (use_div_coeff_1d) then
       call mult_umac_by_1d_coeff(mla,nlevs,umac,div_coeff_1d,div_coeff_half_1d,.false.)
    else if (use_div_coeff_3d) then
       call mult_umac_by_3d_coeff(mla,nlevs,umac,div_coeff_3d,.false.)
    end if

    if (nlevs .gt. 1) then
       do n=2,nlevs
          call create_umac_grown(n,umac(n,:),umac(n-1,:))
       end do
    else
       do n=1,nlevs
          do i=1,dm
             call multifab_fill_boundary(umac(n,i))
          end do
       end do
    end if

    do n = 1, nlevs
       call multifab_destroy(rh(n))
       call multifab_destroy(phi(n))
       call multifab_destroy(alpha(n))
       call multifab_destroy(beta(n))
    end do

    do n = 2,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

  contains

    subroutine divumac(nlevs,umac,rh,dx,ref_ratio,before,divu_rhs)

      use ml_restriction_module, only: ml_cc_restriction
      use probin_module, only: verbose

      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(in   ) :: umac(:,:)
      type(multifab) , intent(inout) :: rh(:)
      real(kind=dp_t), intent(in   ) :: dx(:,:)
      integer        , intent(in   ) :: ref_ratio(:,:)
      logical        , intent(in   ) :: before
      type(multifab ), intent(inout), optional :: divu_rhs(:)

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer :: rhp(:,:,:,:) 
      real(kind=dp_t)          :: rhmax
      integer :: i,dm,lo(rh(nlevs)%dim),hi(rh(nlevs)%dim)

      dm = rh(nlevs)%dim

      do n = 1,nlevs
         do i = 1, rh(n)%nboxes
            if ( multifab_remote(rh(n), i) ) cycle
            ump => dataptr(umac(n,1), i)
            vmp => dataptr(umac(n,2), i)
            rhp => dataptr(rh(n)  , i)
            lo =  lwb(get_box(rh(n), i))
            hi =  upb(get_box(rh(n), i))
            select case (dm)
            case (2)
               call divumac_2d(ump(:,:,1,1), vmp(:,:,1,1), rhp(:,:,1,1), dx(n,:),lo,hi)
            case (3)
               wmp => dataptr(umac(n,3), i)
               call divumac_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               rhp(:,:,:,1), dx(n,:),lo,hi)
            end select
         end do
      end do

      !     NOTE: the sign convention is because the elliptic solver solves
      !            (alpha MINUS del dot beta grad) phi = RHS
      !            Here alpha is zero.

      !     Do rh = divu_rhs - rh
      if (present(divu_rhs)) then
         do n = 1, nlevs
            call multifab_sub_sub(rh(n),divu_rhs(n))
         end do
      end if
      !     ... or rh = -rh
      do n = 1, nlevs
         call multifab_mult_mult_s(rh(n),-ONE)
      end do

      rhmax = norm_inf(rh(nlevs))
      do n = nlevs,2,-1
         call ml_cc_restriction(rh(n-1),rh(n),ref_ratio(n-1,:))
         rhmax = max(rhmax,norm_inf(rh(n-1)))
      end do

      if (parallel_IOProcessor() .and. verbose .ge. 1) then
         if (before) then 
            write(6,1000) 
            write(6,1001) rhmax
         else
            write(6,1002) rhmax
            write(6,1000) 
         end if
      end if

1000  format(' ')
1001  format('... before mac_projection: max of [div (coeff * UMAC) - RHS)]',e15.8)
1002  format('...  after mac_projection: max of [div (coeff * UMAC) - RHS)]',e15.8)

    end subroutine divumac

    subroutine divumac_2d(umac,vmac,rh,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: umac(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(in   ) :: vmac(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) ::   rh(lo(1)  :,lo(2)  :)
      real(kind=dp_t), intent(in   ) ::   dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rh(i,j) = (umac(i+1,j) - umac(i,j)) / dx(1) + &
                      (vmac(i,j+1) - vmac(i,j)) / dx(2)
         end do
      end do

    end subroutine divumac_2d

    subroutine divumac_3d(umac,vmac,wmac,rh,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(in   ) :: vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(inout) ::   rh(lo(1)  :,lo(2)  :,lo(3)  :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               rh(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) + &
                    (vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) + &
                    (wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
            end do
         end do
      end do

    end subroutine divumac_3d

    subroutine mult_umac_by_1d_coeff(mla,nlevs,umac,div_coeff,div_coeff_half,do_mult)

      use ml_restriction_module, only: ml_edge_restriction

      type(ml_layout), intent(in   ) :: mla
      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: umac(:,:)
      real(dp_t)     , intent(in   ) :: div_coeff(:,0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(:,0:)
      logical        , intent(in   ) :: do_mult

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      integer                  :: lo(umac(1,1)%dim)
      integer                  :: i,dm,n

      dm = umac(1,1)%dim

      do n=1,nlevs
         ! Multiply edge velocities by div coeff
         do i = 1, umac(n,1)%nboxes
            if ( multifab_remote(umac(n,1), i) ) cycle
            ump => dataptr(umac(n,1), i)
            vmp => dataptr(umac(n,2), i)
            lo =  lwb(get_box(umac(n,1), i))
            select case (dm)
            case (2)
               call mult_by_1d_coeff_2d(ump(:,:,1,1), vmp(:,:,1,1), &
                                        div_coeff(n,lo(dm):), div_coeff_half(n,lo(dm):), &
                                        do_mult)
            case (3)
               wmp => dataptr(umac(n,3), i)
               call mult_by_1d_coeff_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                        div_coeff(n,lo(dm):), div_coeff_half(n,lo(dm):), &
                                        do_mult)
            end select
         end do

         do i=1,dm
            call multifab_fill_boundary(umac(n,i))
         enddo
      end do

      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),mla%mba%rr(n-1,:),i)
         end do
      end do
      
    end subroutine mult_umac_by_1d_coeff

    subroutine mult_beta_by_1d_coeff(nlevs,beta,div_coeff,div_coeff_half)

      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: beta(:)
      real(dp_t)     , intent(in   ) :: div_coeff(:,0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(:,0:)

      real(kind=dp_t), pointer :: bp(:,:,:,:) 
      integer                  :: lo(beta(1)%dim)
      integer                  :: i,dm,n

      dm = beta(1)%dim

      do n=1,nlevs
         ! Multiply edge coefficients by div coeff
         do i = 1, beta(1)%nboxes
            if ( multifab_remote(beta(n), i) ) cycle
            bp => dataptr(beta(n),i)
            lo =  lwb(get_box(beta(n), i))
            select case (dm)
            case (2)
               call mult_by_1d_coeff_2d(bp(:,:,1,1), bp(:,:,1,2), &
                                        div_coeff(n,lo(dm):), div_coeff_half(n,lo(dm):), &
                                        .true.)
            case (3)
               call mult_by_1d_coeff_3d(bp(:,:,:,1), bp(:,:,:,2), bp(:,:,:,3), &
                                        div_coeff(n,lo(dm):), div_coeff_half(n,lo(dm):), &
                                        .true.)
            end select
         end do
      end do

    end subroutine mult_beta_by_1d_coeff

    subroutine mult_by_1d_coeff_2d(umac,vmac,div_coeff,div_coeff_half,do_mult)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      integer :: j,ny

      ny = size(umac,dim=2)-2

      if (do_mult) then
         do j = 0,ny-1
            umac(:,j) = umac(:,j) * div_coeff(j)
         end do
         do j = 0,ny
            vmac(:,j) = vmac(:,j) * div_coeff_half(j)
         end do
      else
         do j = 0,ny-1 
            umac(:,j) = umac(:,j) / div_coeff(j)
         end do
         do j = 0,ny
            vmac(:,j) = vmac(:,j) / div_coeff_half(j)
         end do
      end if

    end subroutine mult_by_1d_coeff_2d

    subroutine mult_by_1d_coeff_3d(umac,vmac,wmac,div_coeff,div_coeff_half,do_mult)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: wmac(-1:,-1:,-1:)
      real(dp_t)     , intent(in   ) :: div_coeff(0:)
      real(dp_t)     , intent(in   ) :: div_coeff_half(0:)
      logical        , intent(in   ) :: do_mult

      integer :: k,nz

      nz = size(umac,dim=3)-2

      if (do_mult) then
         do k = 0,nz-1 
            umac(:,:,k) = umac(:,:,k) * div_coeff(k)
         end do
         do k = 0,nz-1 
            vmac(:,:,k) = vmac(:,:,k) * div_coeff(k)
         end do
         do k = 0,nz
            wmac(:,:,k) = wmac(:,:,k) * div_coeff_half(k)
         end do
      else
         do k = 0,nz-1 
            umac(:,:,k) = umac(:,:,k) / div_coeff(k)
         end do
         do k = 0,nz-1
            vmac(:,:,k) = vmac(:,:,k) / div_coeff(k)
         end do
         do k = 0,nz
            wmac(:,:,k) = wmac(:,:,k) / div_coeff_half(k)
         end do
      end if

    end subroutine mult_by_1d_coeff_3d

    subroutine mult_umac_by_3d_coeff(mla,nlevs,umac,div_coeff,do_mult)

      use ml_restriction_module, only: ml_edge_restriction

      type(ml_layout), intent(in   ) :: mla
      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: umac(:,:)
      type(multifab) , intent(in   ) :: div_coeff(:)
      logical        , intent(in   ) :: do_mult

      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer ::  dp(:,:,:,:) 
      integer :: i,n,lo(umac(1,1)%dim),hi(umac(1,1)%dim)
      integer :: domlo(umac(1,1)%dim),domhi(umac(1,1)%dim)

      do n=1,nlevs
         domlo =  lwb(ml_layout_get_pd(mla,n))
         domhi =  upb(ml_layout_get_pd(mla,n))

         ! Multiply edge velocities by div coeff
         do i = 1, umac(n,1)%nboxes
            if ( multifab_remote(umac(n,1), i) ) cycle
            ump => dataptr(umac(n,1), i)
            vmp => dataptr(umac(n,2), i)
            wmp => dataptr(umac(n,3), i)
            dp => dataptr(div_coeff(n), i)
            lo =  lwb(get_box(umac(n,1), i))
            hi =  upb(get_box(umac(n,1), i))
            select case (umac(n,1)%dim)
            case (3)
               call mult_by_3d_coeff_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                        dp(:,:,:,1), lo, hi, domlo, domhi, do_mult)
            end select
         end do

         do i=1,dm
            call multifab_fill_boundary(umac(n,i))
         enddo
      enddo

      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),mla%mba%rr(n-1,:),i)
         end do
      end do

    end subroutine mult_umac_by_3d_coeff

    subroutine mult_beta_by_3d_coeff(mla,nlevs,beta,div_coeff)

      type(ml_layout), intent(in   ) :: mla
      integer        , intent(in   ) :: nlevs
      type(multifab) , intent(inout) :: beta(:)
      type(multifab) , intent(in   ) :: div_coeff(:)

      real(kind=dp_t), pointer :: bp(:,:,:,:) 
      real(kind=dp_t), pointer :: dp(:,:,:,:) 
      integer :: i,n,lo(beta(1)%dim),hi(beta(1)%dim)
      integer :: domlo(beta(1)%dim),domhi(beta(1)%dim)

      do n=1,nlevs
         domlo =  lwb(ml_layout_get_pd(mla,n))
         domhi =  upb(ml_layout_get_pd(mla,n))
         
         ! Multiply edge coefficients by div coeff
         do i = 1, beta(n)%nboxes
            if ( multifab_remote(beta(n), i) ) cycle
            bp => dataptr(     beta(n),i)
            dp => dataptr(div_coeff(n),i)
            lo =  lwb(get_box(beta(n), i))
            hi =  upb(get_box(beta(n), i))
            select case (beta(1)%dim)
            case (3)
               call mult_by_3d_coeff_3d(bp(:,:,:,1), bp(:,:,:,2), bp(:,:,:,3), &
                                        dp(:,:,:,1), lo, hi, domlo, domhi, .true.)
            end select
         end do
      enddo

    end subroutine mult_beta_by_3d_coeff

    subroutine mult_by_3d_coeff_3d(umac,vmac,wmac,div_coeff,lo,hi,domlo,domhi,do_mult)

      integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:)
      real(kind=dp_t), intent(inout) ::      umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(inout) ::      vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(inout) ::      wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(dp_t)     , intent(in   ) :: div_coeff(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      logical        , intent(in   ) :: do_mult

      integer :: i,j,k

      if (do_mult) then

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1)+1,hi(1)
                  umac(i,j,k) = umac(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i-1,j,k))
               end do
               if (lo(1).eq.domlo(1)) then
                  umac(lo(1),j,k) = umac(lo(1),j,k) * div_coeff(lo(1),j,k)
               else
                  umac(lo(1),j,k) = umac(lo(1),j,k) * HALF * (div_coeff(lo(1),j,k)+div_coeff(lo(1)-1,j,k))
               end if
               if (hi(1).eq.domhi(1)) then
                  umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) * div_coeff(hi(1),j,k)
               else
                  umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) * HALF * (div_coeff(hi(1)+1,j,k)+div_coeff(hi(1),j,k))
               end if
            end do
         end do

         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               do j = lo(2)+1,hi(2)
                  vmac(i,j,k) = vmac(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i,j-1,k))
               end do
               if (lo(2).eq.domlo(2)) then
                  vmac(i,lo(2),k) = vmac(i,lo(2),k) * div_coeff(i,lo(2),k)
               else
                  vmac(i,lo(2),k) = vmac(i,lo(2),k) * HALF * (div_coeff(i,lo(2),k)+div_coeff(i,lo(2)-1,k))
               end if
               if (hi(2).eq.domhi(2)) then
                  vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) * div_coeff(i,hi(2),k)
               else
                  vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) * HALF * (div_coeff(i,hi(2)+1,k)+div_coeff(i,hi(2),k))
               end if
            end do
         end do

         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               do k = lo(3)+1,hi(3)
                  wmac(i,j,k) = wmac(i,j,k) * HALF * (div_coeff(i,j,k)+div_coeff(i,j,k-1))
               end do
               if (lo(3).eq.domlo(3)) then
                  wmac(i,j,lo(3)) = wmac(i,j,lo(3)) * div_coeff(i,j,lo(3))
               else
                  wmac(i,j,lo(3)) = wmac(i,j,lo(3)) * HALF * (div_coeff(i,j,lo(3))+div_coeff(i,j,lo(3)-1))
               end if
               if (hi(3).eq.domhi(3)) then
                  wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) * div_coeff(i,j,hi(3))
               else
                  wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) * HALF * (div_coeff(i,j,hi(3)+1)+div_coeff(i,j,hi(3)))
               end if
            end do
         end do

      else

         do k = lo(3),hi(3)
            do j = lo(2),hi(2)
               do i = lo(1)+1,hi(1)
                  umac(i,j,k) = umac(i,j,k) / ( HALF * (div_coeff(i,j,k)+div_coeff(i-1,j,k)))
               end do
               if (lo(1).eq.domlo(1)) then
                  umac(lo(1),j,k) = umac(lo(1),j,k) / div_coeff(lo(1),j,k)
               else
                  umac(lo(1),j,k) = umac(lo(1),j,k) / ( HALF * (div_coeff(lo(1),j,k)+div_coeff(lo(1)-1,j,k)))
               end if
               if (hi(1).eq.domhi(1)) then
                  umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) / div_coeff(hi(1),j,k)
               else
                  umac(hi(1)+1,j,k) = umac(hi(1)+1,j,k) / ( HALF * (div_coeff(hi(1)+1,j,k)+div_coeff(hi(1),j,k)))
               end if
            end do
         end do

         do k = lo(3),hi(3)
            do i = lo(1),hi(1)
               do j = lo(2)+1,hi(2)
                  vmac(i,j,k) = vmac(i,j,k) / ( HALF * (div_coeff(i,j,k)+div_coeff(i,j-1,k)))
               end do
               if (lo(2).eq.domlo(2)) then
                  vmac(i,lo(2),k) = vmac(i,lo(2),k) / div_coeff(i,lo(2),k)
               else
                  vmac(i,lo(2),k) = vmac(i,lo(2),k) / ( HALF * (div_coeff(i,lo(2),k)+div_coeff(i,lo(2)-1,k)))
               end if
               if (hi(2).eq.domhi(2)) then
                  vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) / div_coeff(i,hi(2),k)
               else
                  vmac(i,hi(2)+1,k) = vmac(i,hi(2)+1,k) / ( HALF * (div_coeff(i,hi(2)+1,k)+div_coeff(i,hi(2),k)))
               end if
            end do
         end do

         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               do k = lo(3)+1,hi(3)
                  wmac(i,j,k) = wmac(i,j,k) / ( HALF * (div_coeff(i,j,k)+div_coeff(i,j,k-1)))
               end do
               if (lo(3).eq.domlo(3)) then
                  wmac(i,j,lo(3)) = wmac(i,j,lo(3)) / div_coeff(i,j,lo(3))
               else
                  wmac(i,j,lo(3)) = wmac(i,j,lo(3)) / ( HALF * (div_coeff(i,j,lo(3))+div_coeff(i,j,lo(3)-1)))
               end if
               if (hi(3).eq.domhi(3)) then
                  wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) / div_coeff(i,j,hi(3))
               else
                  wmac(i,j,hi(3)+1) = wmac(i,j,hi(3)+1) / ( HALF * (div_coeff(i,j,hi(3)+1)+div_coeff(i,j,hi(3))))
               end if
            end do
         end do

      end if

    end subroutine mult_by_3d_coeff_3d

    subroutine mk_mac_coeffs(nlevs,mla,rho,beta,the_bc_tower)

      use multifab_fill_ghost_module

      integer        , intent(in   ) :: nlevs
      type(ml_layout), intent(in   ) :: mla
      type(multifab ), intent(inout) :: rho(:)
      type(multifab ), intent(inout) :: beta(:)
      type(bc_tower ), intent(in   ) :: the_bc_tower

      real(kind=dp_t), pointer :: bp(:,:,:,:) 
      real(kind=dp_t), pointer :: rp(:,:,:,:) 
      integer :: i,dm,ng,ng_fill

      dm = rho(nlevs)%dim
      ng = rho(nlevs)%ng

      ng_fill = 1
      do n = 2, nlevs
         call multifab_fill_ghost_cells(rho(n),rho(n-1), &
                                        ng_fill,mla%mba%rr(n-1,:), &
                                        the_bc_tower%bc_tower_array(n-1), &
                                        the_bc_tower%bc_tower_array(n  ), &
                                        1,dm+1,1)
      end do

      do n = 1, nlevs
         call multifab_fill_boundary(rho(n))
         do i = 1, rho(n)%nboxes
            if ( multifab_remote(rho(n), i) ) cycle
            rp => dataptr(rho(n) , i)
            bp => dataptr(beta(n), i)
            select case (dm)
            case (2)
               call mk_mac_coeffs_2d(bp(:,:,1,:), rp(:,:,1,1), ng)
            case (3)
               call mk_mac_coeffs_3d(bp(:,:,:,:), rp(:,:,:,1), ng)
            end select
         end do
      end do

    end subroutine mk_mac_coeffs

    subroutine mk_mac_coeffs_2d(beta,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: beta( -1:, -1:,:)
      real(kind=dp_t), intent(inout) ::  rho(-ng:,-ng:)

      integer :: i,j
      integer :: nx,ny

      nx = size(beta,dim=1) - 2
      ny = size(beta,dim=2) - 2

      do j = 0,ny-1
         do i = 0,nx
            beta(i,j,1) = TWO / (rho(i,j) + rho(i-1,j))
         end do
      end do

      do j = 0,ny
         do i = 0,nx-1
            beta(i,j,2) = TWO / (rho(i,j) + rho(i,j-1))
         end do
      end do

    end subroutine mk_mac_coeffs_2d

    subroutine mk_mac_coeffs_3d(beta,rho,ng)

      integer :: ng
      real(kind=dp_t), intent(inout) :: beta( -1:, -1:, -1:,:)
      real(kind=dp_t), intent(inout) ::  rho(-ng:,-ng:,-ng:)

      integer :: i,j,k
      integer :: nx,ny,nz

      nx = size(beta,dim=1) - 2
      ny = size(beta,dim=2) - 2
      nz = size(beta,dim=3) - 2

      do k = 0,nz-1
         do j = 0,ny-1
            do i = 0,nx
               beta(i,j,k,1) = TWO / (rho(i,j,k) + rho(i-1,j,k))
            end do
         end do
      end do

      do k = 0,nz-1
         do j = 0,ny
            do i = 0,nx-1
               beta(i,j,k,2) = TWO / (rho(i,j,k) + rho(i,j-1,k))
            end do
         end do
      end do

      do k = 0,nz
         do j = 0,ny-1
            do i = 0,nx-1
               beta(i,j,k,3) = TWO / (rho(i,j,k) + rho(i,j,k-1))
            end do
         end do
      end do

    end subroutine mk_mac_coeffs_3d

    subroutine mkumac(rh,umac,phi,beta,fine_flx,dx,the_bc_tower,press_comp,ref_ratio)

      use ml_restriction_module, only: ml_edge_restriction

      type(multifab), intent(inout) :: umac(:,:)
      type(multifab), intent(inout) ::   rh(:)
      type(multifab), intent(in   ) ::  phi(:)
      type(multifab), intent(in   ) :: beta(:)
      type(bndry_reg),intent(in   ) :: fine_flx(2:)
      real(dp_t)    , intent(in   ) :: dx(:,:)
      type(bc_tower), intent(in   ) :: the_bc_tower
      integer       , intent(in   ) :: press_comp
      integer       , intent(in   ) :: ref_ratio(:,:)

      integer :: i,dm,nlevs

      type(bc_level)           :: bc
      real(kind=dp_t), pointer :: ump(:,:,:,:) 
      real(kind=dp_t), pointer :: vmp(:,:,:,:) 
      real(kind=dp_t), pointer :: wmp(:,:,:,:) 
      real(kind=dp_t), pointer :: php(:,:,:,:) 
      real(kind=dp_t), pointer ::  bp(:,:,:,:) 
      real(kind=dp_t), pointer :: lxp(:,:,:,:) 
      real(kind=dp_t), pointer :: hxp(:,:,:,:) 
      real(kind=dp_t), pointer :: lyp(:,:,:,:) 
      real(kind=dp_t), pointer :: hyp(:,:,:,:) 
      real(kind=dp_t), pointer :: lzp(:,:,:,:) 
      real(kind=dp_t), pointer :: hzp(:,:,:,:) 

      nlevs = size(rh,dim=1)
      dm = rh(nlevs)%dim

      do n = 1, nlevs
         bc = the_bc_tower%bc_tower_array(n)
         do i = 1, rh(n)%nboxes
            if ( multifab_remote(rh(n), i) ) cycle
            ump => dataptr(umac(n,1), i)
            vmp => dataptr(umac(n,2), i)
            php => dataptr( phi(n), i)
            bp => dataptr(beta(n), i)
            select case (dm)
            case (2)
               if (n > 1) then
                  lxp => dataptr(fine_flx(n)%bmf(1,0), i)
                  hxp => dataptr(fine_flx(n)%bmf(1,1), i)
                  lyp => dataptr(fine_flx(n)%bmf(2,0), i)
                  hyp => dataptr(fine_flx(n)%bmf(2,1), i)
                  call mkumac_2d(ump(:,:,1,1),vmp(:,:,1,1), &
                                 php(:,:,1,1), bp(:,:,1,:), &
                                 lxp(:,:,1,1),hxp(:,:,1,1),lyp(:,:,1,1),hyp(:,:,1,1), &
                                 dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
               else 
                  call mkumac_2d_base(ump(:,:,1,1),vmp(:,:,1,1), & 
                                      php(:,:,1,1), bp(:,:,1,:), &
                                      dx(n,:),bc%ell_bc_level_array(i,:,:,press_comp))
               end if
            case (3)
               wmp => dataptr(umac(n,3), i)
               if (n > 1) then
                  lxp => dataptr(fine_flx(n)%bmf(1,0), i)
                  hxp => dataptr(fine_flx(n)%bmf(1,1), i)
                  lyp => dataptr(fine_flx(n)%bmf(2,0), i)
                  hyp => dataptr(fine_flx(n)%bmf(2,1), i)
                  lzp => dataptr(fine_flx(n)%bmf(3,0), i)
                  hzp => dataptr(fine_flx(n)%bmf(3,1), i)
                  call mkumac_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                 php(:,:,:,1), bp(:,:,:,:), &
                                 lxp(:,:,:,1),hxp(:,:,:,1),lyp(:,:,:,1),hyp(:,:,:,1), &
                                 lzp(:,:,:,1),hzp(:,:,:,1),dx(n,:),&
                                 bc%ell_bc_level_array(i,:,:,press_comp))
               else
                  call mkumac_3d_base(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),& 
                                      php(:,:,:,1), bp(:,:,:,:), dx(n,:), &
                                      bc%ell_bc_level_array(i,:,:,press_comp))
               end if
            end select
         end do

         do d=1,dm
            call multifab_fill_boundary(umac(n,d))
         enddo
         
      end do
      
      do n = nlevs,2,-1
         do i = 1,dm
            call ml_edge_restriction(umac(n-1,i),umac(n,i),ref_ratio(n-1,:),i)
         end do
      end do

    end subroutine mkumac

    subroutine mkumac_2d_base(umac,vmac,phi,beta,dx,press_bc)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gpx,gpy
      integer :: i,j,nx,ny

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2

      if (press_bc(1,1) == BC_NEU) then
         do j = 0,ny-1
            phi(-1,j) = phi(0,j)
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do j = 0,ny-1
            phi(-1,j) = -TWO*phi(0,j) + THIRD * phi(1,j)
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do j = 0,ny-1
            phi(nx,j) = phi(nx-1,j)
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do j = 0,ny-1
            phi(nx,j) = -TWO*phi(nx-1,j) + THIRD * phi(nx-2,j)
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do i = 0,nx-1
            phi(i,-1) = phi(i,0)
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do i = 0,nx-1
            phi(i,-1) = -TWO*phi(i,0) + THIRD * phi(i,1)
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do i = 0,nx-1
            phi(i,ny) = phi(i,ny-1)
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do i = 0,nx-1
            phi(i,ny) = -TWO*phi(i,ny-1) + THIRD * phi(i,ny-2)
         end do
      end if

      do j = 0,ny-1
         do i = 0,nx
            gpx = (phi(i,j) - phi(i-1,j)) / dx(1)
            umac(i,j) = umac(i,j) - beta(i,j,1)*gpx
         end do
      end do

      do i = 0,nx-1
         do j = 0,ny
            gpy = (phi(i,j) - phi(i,j-1)) / dx(2)
            vmac(i,j) = vmac(i,j) - beta(i,j,2)*gpy
         end do
      end do

    end subroutine mkumac_2d_base

    subroutine mkumac_2d(umac,vmac,phi,beta, &
         lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx, &
         dx,press_bc)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,0:), lo_y_flx(0:,:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:,0:), hi_y_flx(0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gpx,gpy
      integer :: i,j,nx,ny

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2

      if (press_bc(1,1) == BC_NEU) then
         do j = 0,ny-1
            phi(-1,j) = phi(0,j)
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do j = 0,ny-1
            phi(-1,j) = -TWO*phi(0,j) + THIRD * phi(1,j)
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do j = 0,ny-1
            phi(nx,j) = phi(nx-1,j)
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do j = 0,ny-1
            phi(nx,j) = -TWO*phi(nx-1,j) + THIRD * phi(nx-2,j)
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do i = 0,nx-1
            phi(i,-1) = phi(i,0)
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do i = 0,nx-1
            phi(i,-1) = -TWO*phi(i,0) + THIRD * phi(i,1)
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do i = 0,nx-1
            phi(i,ny) = phi(i,ny-1)
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do i = 0,nx-1
            phi(i,ny) = -TWO*phi(i,ny-1) + THIRD * phi(i,ny-2)
         end do
      end if

      do j = 0,ny-1
         umac( 0,j) = umac( 0,j) - lo_x_flx(1,j) * dx(1)
         umac(nx,j) = umac(nx,j) + hi_x_flx(1,j) * dx(1)
         do i = 1,nx-1
            gpx = (phi(i,j) - phi(i-1,j)) / dx(1)
            umac(i,j) = umac(i,j) - beta(i,j,1)*gpx
         end do
      end do


      do i = 0,nx-1
         vmac(i, 0) = vmac(i, 0) - lo_y_flx(i,1) * dx(2)
         vmac(i,ny) = vmac(i,ny) + hi_y_flx(i,1) * dx(2)
         do j = 1,ny-1
            gpy = (phi(i,j) - phi(i,j-1)) / dx(2)
            vmac(i,j) = vmac(i,j) - beta(i,j,2)*gpy
         end do
      end do

    end subroutine mkumac_2d

    subroutine mkumac_3d_base(umac,vmac,wmac,phi,beta,dx,press_bc)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: wmac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gpx,gpy,gpz
      integer :: i,j,k,nx,ny,nz

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2
      nz = size(phi,dim=3) - 2

      if (press_bc(1,1) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = phi(0,j,k)
            end do
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = -TWO*phi(0,j,k) + THIRD * phi(1,j,k)
            end do
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = phi(nx-1,j,k)
            end do
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = -TWO*phi(nx-1,j,k) + THIRD * phi(nx-2,j,k)
            end do
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = phi(i,0,k)
            end do
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = -TWO*phi(i,0,k) + THIRD * phi(i,1,k)
            end do
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = phi(i,ny-1,k)
            end do
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = -TWO*phi(i,ny-1,k) + THIRD * phi(i,ny-2,k)
            end do
         end do
      end if
      if (press_bc(3,1) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = phi(i,j,0)
            end do
         end do
      else if (press_bc(3,1) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = -TWO*phi(i,j,0) + THIRD * phi(i,j,1)
            end do
         end do
      end if
      if (press_bc(3,2) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = phi(i,j,nz-1)
            end do
         end do
      else if (press_bc(3,2) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = -TWO*phi(i,j,nz-1) + THIRD * phi(i,j,nz-2)
            end do
         end do
      end if

      do k = 0,nz-1
         do j = 0,ny-1
            do i = 0,nx
               gpx = (phi(i,j,k) - phi(i-1,j,k)) / dx(1)
               umac(i,j,k) = umac(i,j,k) - beta(i,j,k,1)*gpx
            end do
         end do
      end do

      do k = 0,nz-1
         do j = 0,ny
            do i = 0,nx-1
               gpy = (phi(i,j,k) - phi(i,j-1,k)) / dx(2)
               vmac(i,j,k) = vmac(i,j,k) - beta(i,j,k,2)*gpy
            end do
         end do
      end do

      do k = 0,nz
         do j = 0,ny-1
            do i = 0,nx-1
               gpz = (phi(i,j,k) - phi(i,j,k-1)) / dx(3)
               wmac(i,j,k) = wmac(i,j,k) - beta(i,j,k,3)*gpz
            end do
         end do
      end do

    end subroutine mkumac_3d_base

    subroutine mkumac_3d(umac,vmac,wmac,phi,beta,lo_x_flx,hi_x_flx,lo_y_flx,hi_y_flx, &
         lo_z_flx,hi_z_flx,dx,press_bc)

      real(kind=dp_t), intent(inout) :: umac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: vmac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) :: wmac(-1:,-1:,-1:)
      real(kind=dp_t), intent(inout) ::  phi(-1:,-1:,-1:)
      real(kind=dp_t), intent(in   ) :: beta(-1:,-1:,-1:,:)
      real(kind=dp_t), intent(in   ) :: lo_x_flx(:,0:,0:), lo_y_flx(0:,:,0:), lo_z_flx(0:,0:,:)
      real(kind=dp_t), intent(in   ) :: hi_x_flx(:,0:,0:), hi_y_flx(0:,:,0:), hi_z_flx(0:,0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: press_bc(:,:)

      real(kind=dp_t) :: gpx,gpy,gpz
      integer :: i,j,k,nx,ny,nz

      nx = size(phi,dim=1) - 2
      ny = size(phi,dim=2) - 2
      nz = size(phi,dim=3) - 2

      if (press_bc(1,1) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = phi(0,j,k)
            end do
         end do
      else if (press_bc(1,1) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(-1,j,k) = -TWO*phi(0,j,k) + THIRD * phi(1,j,k)
            end do
         end do
      end if
      if (press_bc(1,2) == BC_NEU) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = phi(nx-1,j,k)
            end do
         end do
      else if (press_bc(1,2) == BC_DIR) then
         do k = 0,nz-1
            do j = 0,ny-1
               phi(nx,j,k) = -TWO*phi(nx-1,j,k) + THIRD * phi(nx-2,j,k)
            end do
         end do
      end if
      if (press_bc(2,1) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = phi(i,0,k)
            end do
         end do
      else if (press_bc(2,1) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,-1,k) = -TWO*phi(i,0,k) + THIRD * phi(i,1,k)
            end do
         end do
      end if
      if (press_bc(2,2) == BC_NEU) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = phi(i,ny-1,k)
            end do
         end do
      else if (press_bc(2,2) == BC_DIR) then
         do k = 0,nz-1
            do i = 0,nx-1
               phi(i,ny,k) = -TWO*phi(i,ny-1,k) + THIRD * phi(i,ny-2,k)
            end do
         end do
      end if
      if (press_bc(3,1) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = phi(i,j,0)
            end do
         end do
      else if (press_bc(3,1) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,-1) = -TWO*phi(i,j,0) + THIRD * phi(i,j,1)
            end do
         end do
      end if
      if (press_bc(3,2) == BC_NEU) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = phi(i,j,nz-1)
            end do
         end do
      else if (press_bc(3,2) == BC_DIR) then
         do j = 0,ny-1
            do i = 0,nx-1
               phi(i,j,nz) = -TWO*phi(i,j,nz-1) + THIRD * phi(i,j,nz-2)
            end do
         end do
      end if

      do k = 0,nz-1
         do j = 0,ny-1
            umac( 0,j,k) = umac( 0,j,k) - lo_x_flx(1,j,k) * dx(1)
            umac(nx,j,k) = umac(nx,j,k) + hi_x_flx(1,j,k) * dx(1)
            do i = 1,nx-1
               gpx = (phi(i,j,k) - phi(i-1,j,k)) / dx(1)
               umac(i,j,k) = umac(i,j,k) - beta(i,j,k,1)*gpx
            end do
         end do
      end do

      do k = 0,nz-1
         do i = 0,nx-1
            vmac(i, 0,k) = vmac(i, 0,k) - lo_y_flx(i,1,k) * dx(2)
            vmac(i,ny,k) = vmac(i,ny,k) + hi_y_flx(i,1,k) * dx(2)
            do j = 1,ny-1
               gpy = (phi(i,j,k) - phi(i,j-1,k)) / dx(2)
               vmac(i,j,k) = vmac(i,j,k) - beta(i,j,k,2)*gpy
            end do
         end do
      end do

      do j = 0,ny-1
         do i = 0,nx-1
            wmac(i,j, 0) = wmac(i,j, 0) - lo_z_flx(i,j,1) * dx(3)
            wmac(i,j,nz) = wmac(i,j,nz) + hi_z_flx(i,j,1) * dx(3)
            do k = 1,nz-1
               gpz = (phi(i,j,k) - phi(i,j,k-1)) / dx(3)
               wmac(i,j,k) = wmac(i,j,k) - beta(i,j,k,3)*gpz
            end do
         end do
      end do

    end subroutine mkumac_3d

  end subroutine macproject

end module macproject_module
