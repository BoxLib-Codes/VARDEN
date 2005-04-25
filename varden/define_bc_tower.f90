module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module

  type bc_level

     integer   :: dim    = 0
     integer   :: ngrids = 0
     type(box) :: domain 
     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: dim     = 0
     integer :: nlevels = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()

  end type bc_tower

  contains

  subroutine bc_tower_build(bct,mla,domain_bc,domain_box,nscal)

    type(bc_tower ), intent(  out) :: bct
    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: domain_bc(:,:)
    type(box)      , intent(in   ) :: domain_box(:)
    integer        , intent(in   ) :: nscal

    integer :: ngrids
    integer :: default_value

    bct%nlevels = mla%nlevel
    bct%dim     = mla%dim

    allocate(bct%bc_tower_array(bct%nlevels))
    do i = 1,bct%nlevels
       ngrids = layout_nboxes(mla%la(i))
       bct%bc_tower_array(i)%dim    = bct%dim
       bct%bc_tower_array(i)%ngrids = ngrids
       bct%bc_tower_array(i)%domain = domain_box(i)

       allocate(bct%bc_tower_array(i)%phys_bc_level_array(0:ngrids,bct%dim,2))
       default_value = INTERIOR
       call phys_bc_level_build(bct%bc_tower_array(i)%phys_bc_level_array,mla%la(i), &
                                domain_bc,default_value)

       allocate(bct%bc_tower_array(i)%adv_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+1))
       default_value = INTERIOR
       call adv_bc_level_build(bct%bc_tower_array(i)%adv_bc_level_array, &
                               bct%bc_tower_array(i)%phys_bc_level_array,default_value)

       allocate(bct%bc_tower_array(i)%ell_bc_level_array(0:ngrids,bct%dim,2,bct%dim+nscal+1))
       default_value = BC_INT
       call ell_bc_level_build(bct%bc_tower_array(i)%ell_bc_level_array, &
                               bct%bc_tower_array(i)%phys_bc_level_array,default_value)
    end do

  end subroutine bc_tower_build

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: i,dm

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do dm = 1,layout_dim(la_level)
       phys_bc_level(i,dm,1) = domain_bc(dm,1)
       phys_bc_level(i,dm,2) = domain_bc(dm,2)
    end do

    do i = 1,layout_nboxes(la_level)
       bx = layout_get_box(la_level,i)
       do dm = 1,layout_dim(la_level)
          if (bx%lo(dm) == pd%lo(dm)) phys_bc_level(i,dm,1) = domain_bc(dm,1)
          if (bx%hi(dm) == pd%hi(dm)) phys_bc_level(i,dm,2) = domain_bc(dm,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm
    integer :: n,d,i

    adv_bc_level = default_value

!    *** 2-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : density
!   COMP = 4  : tracer

!    *** 3-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : z-velocity
!   COMP = 4  : density
!   COMP = 5  : tracer

    dm = size(adv_bc_level,dim=2)

    do n  = 0,size(adv_bc_level,dim=1)-1
    do d  = 1,dm
    do i  = 1,2
       if (phys_bc_level(n,d,i) == SLIP_WALL) then
          adv_bc_level(n,d,i,1:dm) = HOEXTRAP      ! tangential vel.
          adv_bc_level(n,d,i,   d) = EXT_DIR       ! normal vel.
          adv_bc_level(n,d,i,dm+1) = HOEXTRAP      ! density
          adv_bc_level(n,d,i,dm+2) = HOEXTRAP      ! tracer
          adv_bc_level(n,d,i,dm+3) = FOEXTRAP      ! pressure

       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          adv_bc_level(n,d,i,1:dm) = EXT_DIR       ! vel.
          adv_bc_level(n,d,i,dm+1) = HOEXTRAP      ! density
          adv_bc_level(n,d,i,dm+2) = HOEXTRAP      ! tracer
          adv_bc_level(n,d,i,dm+3) = FOEXTRAP      ! pressure
   
       else if (phys_bc_level(n,d,i) == INLET) then
          adv_bc_level(n,d,i,1:dm) = EXT_DIR       ! vel.
          adv_bc_level(n,d,i,dm+1) = EXT_DIR       ! density
          adv_bc_level(n,d,i,dm+2) = EXT_DIR       ! tracer
          adv_bc_level(n,d,i,dm+3) = FOEXTRAP      ! pressure

       else if (phys_bc_level(n,d,i) == OUTLET) then
          adv_bc_level(n,d,i,1:dm) = FOEXTRAP      ! vel.
          adv_bc_level(n,d,i,dm+1) = FOEXTRAP      ! density
          adv_bc_level(n,d,i,dm+2) = FOEXTRAP      ! tracer
          adv_bc_level(n,d,i,dm+3) = EXT_DIR       ! pressure

       else if (phys_bc_level(n,d,i) == SYMMETRY) then
          adv_bc_level(n,d,i,1:dm) = REFLECT_EVEN  ! tangential vel.
          adv_bc_level(n,d,i,   d) = REFLECT_ODD   ! normal vel.
          adv_bc_level(n,d,i,dm+1) = REFLECT_EVEN  ! density
          adv_bc_level(n,d,i,dm+2) = REFLECT_EVEN  ! tracer
          adv_bc_level(n,d,i,dm+3) = REFLECT_EVEN  ! pressure
       end if
    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm
    integer :: n,d,i

    ell_bc_level = default_value

!    *** 2-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : density
!   COMP = 4  : tracer
!   COMP = 5  : pressure

!    *** 3-D ***
!   COMP = 1  : x-velocity
!   COMP = 2  : y-velocity
!   COMP = 3  : z-velocity
!   COMP = 4  : density
!   COMP = 5  : tracer
!   COMP = 6  : pressure

    dm = size(ell_bc_level,dim=2)

    do n = 0,size(ell_bc_level,dim=1)-1
    do d = 1,dm
    do i = 1,2
       if (phys_bc_level(n,d,i) == SLIP_WALL) then
          ell_bc_level(n,d,i,1:dm) = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,   d) = BC_DIR   ! normal vel.
          ell_bc_level(n,d,i,dm+1) = BC_NEU   ! density
          ell_bc_level(n,d,i,dm+2) = BC_NEU   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_NEU   ! pressure
       else if (phys_bc_level(n,d,i) == NO_SLIP_WALL) then
          ell_bc_level(n,d,i,1:dm) = BC_DIR   ! vel.
          ell_bc_level(n,d,i,dm+1) = BC_NEU   ! density
          ell_bc_level(n,d,i,dm+2) = BC_NEU   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_NEU   ! pressure
       else if (phys_bc_level(n,d,i) == INLET) then
          ell_bc_level(n,d,i,1:dm) = BC_DIR   ! vel.
          ell_bc_level(n,d,i,dm+1) = BC_DIR   ! density
          ell_bc_level(n,d,i,dm+2) = BC_DIR   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_NEU   ! pressure
       else if (phys_bc_level(n,d,i) == OUTLET) then
          ell_bc_level(n,d,i,1:dm) = BC_NEU   ! vel.
          ell_bc_level(n,d,i,dm+1) = BC_NEU   ! density
          ell_bc_level(n,d,i,dm+2) = BC_NEU   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_DIR   ! pressure
       else if (phys_bc_level(n,d,i) == SYMMETRY) then
          ell_bc_level(n,d,i,1:dm) = BC_NEU   ! tangential vel.
          ell_bc_level(n,d,i,   d) = BC_DIR   ! normal vel.
          ell_bc_level(n,d,i,dm+1) = BC_NEU   ! density
          ell_bc_level(n,d,i,dm+2) = BC_NEU   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_NEU   ! pressure
       else if (phys_bc_level(n,d,i) == PERIODIC) then
          ell_bc_level(n,d,i,1:dm) = BC_PER   ! tangential vel.
          ell_bc_level(n,d,i,   d) = BC_PER   ! normal vel.
          ell_bc_level(n,d,i,dm+1) = BC_PER   ! density
          ell_bc_level(n,d,i,dm+2) = BC_PER   ! tracer
          ell_bc_level(n,d,i,dm+3) = BC_PER   ! pressure
       end if
    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
