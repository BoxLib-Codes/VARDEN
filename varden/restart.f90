module restart_module

  use bl_error_module
  use bl_string_module
  use bl_IO_module
  use bl_types
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

  subroutine do_restart(nscal,ng_cell,restart_int,mla,uold,sold,gp,p,dx,time,dt)

    integer          , intent(in    ) :: nscal,ng_cell,restart_int
    real(dp_t)       , intent(   out) :: time,dt
    real(dp_t)       , pointer        :: dx(:,:)
    type(multifab)   , pointer        :: uold(:),sold(:),gp(:),p(:)
    type(ml_layout)                   :: mla

    type(ml_boxarray)                 :: mba
    type(multifab)   , pointer        :: chkdata(:)
    character(len=7)                  :: sd_name
    integer          , pointer        :: rrs(:)
    integer                           :: n,nlevs,dm

    write(unit=sd_name,fmt='("chk",i4.4)') restart_int
    print *,'Reading ',sd_name,' to get state data for restart'

    call checkpoint_read(mba, chkdata, p, sd_name, time, dt, nlevs, rrs, dx)
    call ml_layout_build(mla,mba)

    dm = mba%dim

    allocate(uold(nlevs),sold(nlevs),gp(nlevs))

    do n = 1,nlevs
     call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
     call multifab_build(   sold(n), mla%la(n),    dm, ng_cell)
     call multifab_build(     gp(n), mla%la(n),    dm,       1)
    end do

    do n = 1,nlevs
     call setval(uold(n),0.0_dp_t,all=.true.)
     call setval(sold(n),0.0_dp_t,all=.true.)
     call setval(  gp(n),0.0_dp_t,all=.true.)
    end do

    do n = 1,nlevs
       call multifab_copy_c(uold(n),1,chkdata(n),1     ,dm)
       call multifab_copy_c(sold(n),1,chkdata(n),1+  dm,nscal)
       call multifab_copy_c(  gp(n),1,chkdata(n),1+2*dm,dm)
    end do

    do n = 1,nlevs
!      Synchronize incoming data
       call multifab_fill_boundary(  gp(n))
    end do
  end subroutine do_restart

end module restart_module
