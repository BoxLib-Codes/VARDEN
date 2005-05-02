module restart_module

  use bl_error_module
  use bl_string_module
  use bl_IO_module
  use bl_types
  use box_util_module
  use define_bc_module
  use setbc_module
  use fab_module
  use fabio_module
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use parallel

  implicit none

contains

  subroutine fill_restart_data(nscal,ng_cell,restart_int,uold,sold,gp,p,dx,time,dt)

    integer          , intent(in    ) :: nscal,ng_cell,restart_int
    real(dp_t)       , intent(   out) :: time,dt
    real(dp_t)       , pointer        :: dx(:,:)
    type(multifab)   , pointer        :: uold(:),sold(:),gp(:),p(:)

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chk_p(:)
    character(len=7)                  :: sd_name
    integer          , pointer        :: rrs(:)
    integer                           :: n,nlevs,dm


    write(unit=sd_name,fmt='("chk",i4.4)') restart_int
    print *,'Reading ',sd_name,' to get state data for restart'
    call checkpoint_read(chkdata, chk_p, sd_name, time, dt, nlevs)

    dm = chkdata(1)%dim

    do n = 1,nlevs
       call multifab_copy_c(uold(n),1,chkdata(n),1     ,dm)
       call multifab_copy_c(sold(n),1,chkdata(n),1+  dm,nscal)
       call multifab_copy_c(  gp(n),1,chkdata(n),1+2*dm,dm)
       call multifab_copy_c(   p(n),1,  chk_p(n),1     ,1)
       call multifab_destroy(chkdata(n))
       call multifab_destroy(chk_p(n))
    end do

    call destroy(chkdata(1)%la)
    call destroy(chk_p(1)%la)

    deallocate(chkdata,chk_p)

    ! Synchronize incoming data
    do n = 1,nlevs
       call multifab_fill_boundary(uold(n))
       call multifab_fill_boundary(sold(n))
       call multifab_fill_boundary(  gp(n))
       call multifab_fill_boundary(   p(n))
    end do

  end subroutine fill_restart_data

end module restart_module
