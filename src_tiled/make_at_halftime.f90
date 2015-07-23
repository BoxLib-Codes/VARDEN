module rhohalf_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: make_at_halftime

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime(mla,rhohalf,sold,snew,in_comp,out_comp,the_bc_level)

    use multifab_physbc_module
    use ml_restrict_fill_module

    type(multifab) , intent(inout) :: rhohalf(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(in   ) :: snew(:)
    integer        , intent(in   ) :: in_comp,out_comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    real(kind=dp_t), pointer:: rhp(:,:,:,:)
    real(kind=dp_t), pointer:: rop(:,:,:,:)
    real(kind=dp_t), pointer:: rnp(:,:,:,:)
    integer   :: lo(get_dim(rhohalf(1))),hi(get_dim(rhohalf(1)))
    integer   :: nlevs,ng_h,ng_o,dm,i,n

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"make_at_halftime")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_h = nghost(rhohalf(1))
    ng_o = nghost(sold(1))

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(rhp,rop,rnp,lo,hi)

    do n = 1, nlevs
       call mfiter_build(mfi, rhohalf(n), tiling=.true.)
       do while (more_tile(mfi))
          i = get_fab_index(mfi)
          
          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i = 1, nfabs(rhohalf(n))
          rhp => dataptr(rhohalf(n), i)
          rop => dataptr(sold(n), i)
          rnp => dataptr(snew(n), i)
          lo =  lwb(get_box(rhohalf(n), i))
          hi =  upb(get_box(rhohalf(n), i))
          select case (dm)
          case (2)
             call make_at_halftime_2d(rhp(:,:,1,out_comp),rop(:,:,1,in_comp), &
                                      rnp(:,:,1,in_comp),lo,hi,ng_h,ng_o,tlo,thi)
          case (3)
             call make_at_halftime_3d(rhp(:,:,:,out_comp),rop(:,:,:,in_comp), &
                                      rnp(:,:,:,in_comp),lo,hi,ng_h,ng_o,tlo,thi)
          end select
       end do
    end do
    !$omp end parallel

    call ml_restrict_and_fill(nlevs, rhohalf, mla%mba%rr, the_bc_level, &
         icomp=out_comp, bcomp=dm+in_comp, nc=1, ng=ng_h)

    call destroy(bpt)

  end subroutine make_at_halftime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_2d(rhohalf,rhoold,rhonew,glo,ghi,ng_half,ng_old,tlo,thi)

    use bl_constants_module

    integer         , intent(in   ) :: glo(:),ghi(:),ng_half,ng_old,tlo(:),thi(:)
    real (kind=dp_t), intent(  out) :: rhohalf(glo(1)-ng_half:,glo(2)-ng_half:)
    real (kind=dp_t), intent(in   ) :: rhoold(glo(1)-ng_old:,glo(2)-ng_old:)
    real (kind=dp_t), intent(in   ) :: rhonew(glo(1)-ng_old:,glo(2)-ng_old:)

    !  Local variables
    integer :: i, j

    do j = tlo(2),thi(2)
       do i = tlo(1),thi(1)
          rhohalf(i,j) = HALF * (rhoold(i,j) + rhonew(i,j))
       end do
    end do

  end subroutine make_at_halftime_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_at_halftime_3d(rhohalf,rhoold,rhonew,glo,ghi,ng_half,ng_old,tlo,thi)

    use bl_constants_module

    integer         , intent(in   ) :: glo(:),ghi(:),ng_half,ng_old,tlo(:),thi(:)
    real (kind=dp_t), intent(  out) :: rhohalf(glo(1)-ng_half:,glo(2)-ng_half:,glo(3)-ng_half:)
    real (kind=dp_t), intent(in   ) :: rhoold(glo(1)-ng_old:,glo(2)-ng_old:,glo(3)-ng_old:)
    real (kind=dp_t), intent(in   ) :: rhonew(glo(1)-ng_old:,glo(2)-ng_old:,glo(3)-ng_old:)

    ! Local variables
    integer :: i, j, k

    do k = tlo(3),thi(3)
       do j = tlo(2),thi(2)
          do i = tlo(1),thi(1)
             rhohalf(i,j,k) = HALF * (rhoold(i,j,k) + rhonew(i,j,k))
          end do
       end do
    end do

  end subroutine make_at_halftime_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module rhohalf_module
