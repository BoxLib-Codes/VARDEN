module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: fill_restart_data

contains

  subroutine fill_restart_data(restart_int,mba,chkdata,chk_p,time,dt)

    use checkpoint_module
    use probin_module, only : MAX_ALLOWED_LEVS, check_base_name

    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chk_p(:)
    character(len=5)                  :: check_index
    character(len=256)                :: check_file_name

    integer                           :: n,nlevs,dm
    integer                           :: rrs(MAX_ALLOWED_LEVS)

    write(unit=check_index,fmt='(i5.5)') restart_int
    check_file_name = trim(check_base_name) // check_index
    if ( parallel_IOProcessor() ) &
       print *,'Reading ',trim(check_file_name),' to get state data for restart'
    call checkpoint_read(chkdata, chk_p, check_file_name, rrs, time, dt, nlevs)

    dm = get_dim(chkdata(1))

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = rrs(n-1)
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

  end subroutine fill_restart_data

end module restart_module
