module pre_advance_module

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use velpred_module
  use mkforce_module
  use setbc_module

  implicit none

contains

  subroutine advance_premac(nlevs,uold,sold,lapu,umac,gp,ext_vel_force,dx,time,dt, &
                            the_bc_level,visc_coef,use_godunov_debug,use_minion)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(in   ) :: lapu(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: gp(:)
    type(multifab) , intent(in   ) :: ext_vel_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time,dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: visc_coef
    logical        , intent(in)    :: use_godunov_debug
    logical        , intent(in)    :: use_minion

    type(multifab), allocatable :: vel_force(:)
    integer                     :: dm,n
    real(kind=dp_t)             :: visc_fac

    allocate(vel_force(nlevs))

    dm = uold(1)%dim
    visc_fac = 1.0d0

    do n = 1, nlevs
       call multifab_build(vel_force(n),ext_vel_force(n)%la,dm,1)
       call setval(vel_force(n),0.0_dp_t,all=.true.)
    enddo

    call mkvelforce(nlevs,vel_force,ext_vel_force,sold,gp,uold,lapu,dx,visc_coef,visc_fac)

    if(use_godunov_debug) then
       call velpred(nlevs,uold,umac,vel_force,dx,dt,the_bc_level,use_minion)
    else
       call velpred_debug(nlevs,uold,umac,vel_force,dx,dt,the_bc_level,use_minion)
    endif

    do n = 1, nlevs
       call multifab_destroy(vel_force(n))
    enddo

  end subroutine advance_premac

end module pre_advance_module
