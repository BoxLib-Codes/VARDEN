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
  use checkpoint_module
  use parallel

  implicit none

contains

  subroutine do_restart(dm,nscal,ng_cell,restart_int,la_tower,uold,sold,gp,p,rrs,dx,time,dt)

    integer   , intent(in    ) :: dm,nscal,ng_cell,restart_int
    real(dp_t), intent(   out) :: time,dt
    integer   , pointer  :: rrs(:)
    real(dp_t), pointer  :: dx(:,:)
    type(multifab), pointer  :: uold(:),sold(:),gp(:),p(:)
    type(layout)  , pointer  :: la_tower(:)
    type(multifab), pointer        :: chkdata(:)
    character(len=7)               :: sd_name
    integer                        :: n_chk_comps
    integer                        :: n,nlevs
    integer                        :: ref_ratio(dm)
    logical                        :: nodal(2)

    write(unit=sd_name,fmt='("chk",i4.4)') restart_int
    print *,'Reading ',sd_name,' to get state data for restart'
    n_chk_comps = 2*dm + nscal

    call checkpoint_read(chkdata, p, sd_name, time, dt, nlevs, rrs, dx)

    allocate(la_tower(nlevs))
    call build(la_tower(1),get_boxarray(chkdata(1)))
    ref_ratio = 2
    do n = 2,nlevs
      call layout_build_pn(la_tower(n),la_tower(n-1),get_boxarray(chkdata(n)),ref_ratio)
    end do

    allocate(uold(nlevs),sold(nlevs),gp(nlevs))

    do n = 1,nlevs
     call multifab_build(   uold(n), la_tower(n),    dm, ng_cell)
     call multifab_build(   sold(n), la_tower(n),    dm, ng_cell)
     call multifab_build(     gp(n), la_tower(n),    dm,       1)
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
