module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module

  implicit none

  type bc_level

     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: max_level_built = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()
     integer       , pointer :: domain_bc(:,:) => Null()

  end type bc_tower

  private

  public :: bc_level, bc_tower, bc_tower_init, bc_tower_level_build, bc_tower_destroy

  contains

  subroutine bc_tower_init(bct,num_levs,dm,phys_bc_in)

    type(bc_tower ), intent(  out) :: bct
    integer        , intent(in   ) :: num_levs
    integer        , intent(in   ) :: dm
    integer        , intent(in   ) :: phys_bc_in(:,:)

    allocate(bct%bc_tower_array(num_levs))
    allocate(bct%domain_bc(dm,2))

    bct%domain_bc(:,:) = phys_bc_in(:,:)

  end subroutine bc_tower_init

  subroutine bc_tower_level_build(bct,n,la)

    use probin_module, only : nscal

    type(bc_tower ), intent(inout) :: bct
    integer        , intent(in   ) :: n
    type(layout)   , intent(in   ) :: la

    integer :: ngrids,default_value,dm

    if (associated(bct%bc_tower_array(n)%phys_bc_level_array)) then
      deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
      deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
      deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
      bct%bc_tower_array(n)%phys_bc_level_array => NULL()
      bct%bc_tower_array(n)%adv_bc_level_array => NULL()
      bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end if

    ngrids = layout_nlocal(la)
    dm = layout_dim(la)

    allocate(bct%bc_tower_array(n)%phys_bc_level_array(0:ngrids,dm,2))
    default_value = INTERIOR
    call phys_bc_level_build(bct%bc_tower_array(n)%phys_bc_level_array,la, &
                             bct%domain_bc,default_value)

    ! Here we allocate dm components for velocity
    !                    nscal components for scalars
    !                        1 component  for generic extrap
    !                        1 component  for pressure
    allocate(bct%bc_tower_array(n)%adv_bc_level_array(0:ngrids,dm,2,dm+nscal+2))
    default_value = INTERIOR
    call adv_bc_level_build(bct%bc_tower_array(n)%adv_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

    allocate(bct%bc_tower_array(n)%ell_bc_level_array(0:ngrids,dm,2,dm+nscal+1))
    default_value = BC_INT
    call ell_bc_level_build(bct%bc_tower_array(n)%ell_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

     bct%max_level_built = n

  end subroutine bc_tower_level_build

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: n

    do n = 1,bct%max_level_built
       deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
       bct%bc_tower_array(n)%phys_bc_level_array => NULL()
       bct%bc_tower_array(n)%adv_bc_level_array => NULL()
       bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end do
    deallocate(bct%bc_tower_array)

    deallocate(bct%domain_bc)

  end subroutine bc_tower_destroy

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: d,i

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do d = 1,layout_dim(la_level)
       phys_bc_level(i,d,1) = domain_bc(d,1)
       phys_bc_level(i,d,2) = domain_bc(d,2)
    end do

    do i = 1,layout_nlocal(la_level)
       bx = layout_get_box(la_level,global_index(la_level,i))
       do d = 1,layout_dim(la_level)
          if (lwb(bx,d) == lwb(pd,d)) phys_bc_level(i,d,1) = domain_bc(d,1)
          if (upb(bx,d) == upb(pd,d)) phys_bc_level(i,d,2) = domain_bc(d,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value)

    use probin_module, only : nscal,extrap_comp

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm, ns
    integer :: igrid,d,lohi
    integer :: press_comp

    adv_bc_level = default_value

!    *** 2-D ***
! last index represents
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : density
!   COMP = 4 - (4+nscal-1) : tracer_i
!   next-to-last COMP   : pressure
!   last COMP   : generic extrap

!    *** 3-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : z-velocity
!   COMP = 4  : density
!   COMP = 5 - (5+nscal-1) : tracer_i
!   next-to-last COMP   : pressure
!   last COMP   : generic extrap
 
    dm = size(adv_bc_level,dim=2)

    press_comp  = dm+nscal+1
    extrap_comp = press_comp+1

    do igrid = 0, size(adv_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2

       if (phys_bc_level(igrid,d,lohi) == SLIP_WALL) then

          adv_bc_level(igrid,d,lohi,1:dm) = HOEXTRAP      ! tangential vel.
          adv_bc_level(igrid,d,lohi,   d) = EXT_DIR       ! normal vel.
          do ns = 1, nscal
             adv_bc_level(igrid,d,lohi,dm+ns) = HOEXTRAP      ! density and tracers
          enddo
          adv_bc_level(igrid,d,lohi, press_comp) = FOEXTRAP      ! pressure
          adv_bc_level(igrid,d,lohi,extrap_comp) = FOEXTRAP      ! generic extrap bc

       else if (phys_bc_level(igrid,d,lohi) == NO_SLIP_WALL) then

          adv_bc_level(igrid,d,lohi,1:dm) = EXT_DIR       ! vel.
          do ns = 1, nscal
             adv_bc_level(igrid,d,lohi,dm+ns) = HOEXTRAP      ! density and tracers
          enddo
          adv_bc_level(igrid,d,lohi, press_comp) = FOEXTRAP      ! pressure
          adv_bc_level(igrid,d,lohi,extrap_comp) = FOEXTRAP      ! generic extrap bc

       else if (phys_bc_level(igrid,d,lohi) == INLET) then

          adv_bc_level(igrid,d,lohi,1:dm) = EXT_DIR       ! vel.
          do ns = 1, nscal
             adv_bc_level(igrid,d,lohi,dm+ns) = EXT_DIR       ! density and tracers
          enddo
          adv_bc_level(igrid,d,lohi,press_comp) = FOEXTRAP      ! pressure
          adv_bc_level(igrid,d,lohi,extrap_comp) = FOEXTRAP      ! generic extrap bc

       else if (phys_bc_level(igrid,d,lohi) == OUTLET) then

          adv_bc_level(igrid,d,lohi,1:dm) = FOEXTRAP      ! vel.
          do ns = 1, nscal
             adv_bc_level(igrid,d,lohi,dm+ns) = FOEXTRAP      ! density and tracers
          enddo
          adv_bc_level(igrid,d,lohi,press_comp) = EXT_DIR       ! pressure
          adv_bc_level(igrid,d,lohi,extrap_comp) = FOEXTRAP      ! generic extrap bc

       else if (phys_bc_level(igrid,d,lohi) == SYMMETRY) then

          adv_bc_level(igrid,d,lohi,1:dm) = REFLECT_EVEN  ! tangential vel.
          adv_bc_level(igrid,d,lohi,   d) = REFLECT_ODD   ! normal vel.
          do ns = 1, nscal
             adv_bc_level(igrid,d,lohi,dm+ns) = REFLECT_EVEN  ! density and tracers
          enddo
          adv_bc_level(igrid,d,lohi,press_comp) = EXT_DIR       ! pressure
          adv_bc_level(igrid,d,lohi,extrap_comp) = REFLECT_EVEN ! generic extrap bc

       end if

    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value)

    use probin_module, only : nscal

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm,ns
    integer :: comp,d,lohi
    integer :: press_comp

    ell_bc_level = default_value
 
    dm = size(ell_bc_level,dim=2)

    press_comp = dm+nscal+1

!    *** 2-D ***
! last index represents
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : density
!   COMP = 4 - (4+nscal-1) : tracer_i
!   last COMP   : pressure

!    *** 3-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : z-velocity
!   COMP = 4  : density
!   COMP = 5 - (5+nscal-1) : tracer_i
!   last COMP   : pressure

    do comp = 0, size(ell_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2
       if (phys_bc_level(comp,d,lohi) == SLIP_WALL) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_NEU   ! tangential vel.
          ell_bc_level(comp,d,lohi,   d) = BC_DIR   ! normal vel.
          do ns = 1, nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_NEU   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_NEU   ! pressure

       else if (phys_bc_level(comp,d,lohi) == NO_SLIP_WALL) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_DIR   ! vel.
          do ns = 1, nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_NEU   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_NEU   ! pressure

       else if (phys_bc_level(comp,d,lohi) == INLET) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_DIR   ! vel.
          do ns = 1,nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_DIR   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_NEU   ! pressure

       else if (phys_bc_level(comp,d,lohi) == OUTLET) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_NEU   ! vel.
          do ns = 1, nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_NEU   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_DIR   ! pressure

       else if (phys_bc_level(comp,d,lohi) == SYMMETRY) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_NEU   ! tangential vel.
          ell_bc_level(comp,d,lohi,   d) = BC_DIR   ! normal vel.
          do ns = 1, nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_NEU   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_NEU   ! pressure

       else if (phys_bc_level(comp,d,lohi) == PERIODIC) then
          ell_bc_level(comp,d,lohi,1:dm) = BC_PER   ! tangential vel.
          ell_bc_level(comp,d,lohi,   d) = BC_PER   ! normal vel.
          do ns = 1, nscal
             ell_bc_level(comp,d,lohi,dm+ns) = BC_PER   ! density and tracers
          enddo
          ell_bc_level(comp,d,lohi,press_comp) = BC_PER   ! pressure
       end if
    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
