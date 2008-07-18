module pre_advance_module

  use bl_types

  use define_bc_module
  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: advance_premac

contains

  subroutine advance_premac(mla,uold,sold,lapu,umac,gp,ext_vel_force,dx,dt, &
                            the_bc_level)

    use velpred_module
    use mkforce_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(in   ) :: lapu(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: gp(:)
    type(multifab) , intent(in   ) :: ext_vel_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab)  :: vel_force(mla%nlevel)
    integer         :: dm,n,nlevs
    real(kind=dp_t) :: visc_fac

    dm    = uold(1)%dim
    nlevs = mla%nlevel

    do n = 1, nlevs
       call multifab_build(vel_force(n),ext_vel_force(n)%la,dm,1)
       call setval(vel_force(n),0.0_dp_t,all=.true.)
    enddo

    visc_fac = 1.0d0
    call mkvelforce(nlevs,vel_force,ext_vel_force,sold,gp,lapu,visc_fac)

    call velpred(nlevs,uold,umac,vel_force,dx,dt,the_bc_level,mla)

    do n = 1, nlevs
       call multifab_destroy(vel_force(n))
    enddo

  end subroutine advance_premac

end module pre_advance_module
