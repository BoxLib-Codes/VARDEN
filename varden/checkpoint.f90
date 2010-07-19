module checkpoint_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(nlevs_in, dirname, mfs, mfs_nodal, rrs, time_in, dt_in, verbose)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
    use parallel

    integer       , intent(in) :: nlevs_in
    type(multifab), intent(in) :: mfs(:), mfs_nodal(:)
    integer       , intent(in) :: rrs(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t), intent(in) :: time_in, dt_in
    integer        , intent(in) :: verbose
    integer :: n
    character(len=128) :: header, sd_name, sd_name_nodal
    integer :: nc, un, dm
    integer, allocatable ::  lo(:),  hi(:)
    type(box) :: lbbox

    integer         :: nlevs
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ time
    namelist /chkpoint/ dt
    namelist /chkpoint/ nlevs

    if ( parallel_IOProcessor() ) call fabio_mkdir(dirname)

    call parallel_barrier() ! All CPUs have to wait till the directory is built.

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

    write(unit=sd_name_nodal, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal, rrs(:,1), sd_name_nodal)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      print *,'Writing    state to checkpoint file ',trim(sd_name)
      print *,'Writing pressure to checkpoint file ',trim(sd_name_nodal)
      print *,' '
    end if
    
    nc = ncomp(mfs(1))
    dm = get_dim(mfs(1))
    allocate(lo(dm),hi(dm))
    lbbox = bbox(get_boxarray(mfs(1)))

    lo = lwb(lbbox); hi = upb(lbbox)

    time = time_in
      dt =   dt_in

    if (parallel_IOProcessor()) then
       header = "Header"
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       nlevs = nlevs_in
       write(unit=un, nml = chkpoint)
       do n = 1,nlevs-1
          write(unit=un,fmt=*) rrs(n,1)
       end do
       close(un)
    end if
    
    deallocate(lo,hi)

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodal, dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab  ),                pointer :: mfs(:), mfs_nodal(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    integer          ,               pointer :: rrs(:)

    integer :: n
    character(len=128) :: header, sd_name
    integer :: nc, un, nl, dm
    type(box) :: lbbox

    integer         :: nlevs
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)
    allocate(rrs(nlevs-1))
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

!   Read the pressure data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_nodal, sd_name)
    
    nl = nlevs
    nc = ncomp(mfs(1))
    dm = get_dim(mfs(1))
    lbbox = bbox(get_boxarray(mfs(1)))

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

  end subroutine checkpoint_read

end module checkpoint_module
