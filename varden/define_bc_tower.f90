module define_bc_module

  use bl_types
  use layout_module
  use bc_module

  type bc_level

     integer :: dim    = 0
     integer :: ngrids = 0
     integer, pointer :: bc_level_array(:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: dim     = 0
     integer :: nlevels = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()

  end type bc_tower

  contains

  subroutine bc_tower_build(bct,la_tower,domain_bc,default_value)

    type(bc_tower), intent(  out) :: bct
    type(layout)  , intent(in   ) :: la_tower(:)
    integer       , intent(in   ) :: domain_bc(:,:)
    integer       , intent(in   ) :: default_value

    integer :: ngrids

    bct%nlevels = size(la_tower,dim=1)
    bct%dim     = layout_dim(la_tower(1))

    allocate(bct%bc_tower_array(bct%nlevels))
    do i = 1,bct%nlevels
       ngrids = layout_nboxes(la_tower(i))
       bct%bc_tower_array(i)%dim    = bct%dim
       bct%bc_tower_array(i)%ngrids = ngrids
       allocate(bct%bc_tower_array(i)%bc_level_array(ngrids,bct%dim,2))
       call bc_level_build(bct%bc_tower_array(i)%bc_level_array,la_tower(i), &
                           domain_bc,default_value)
    end do

  end subroutine bc_tower_build

  subroutine bc_level_build(bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: bc_level(:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: i,dm

    pd = layout_get_pd(la_level) 

    bc_level = default_value

    do i = 1,layout_nboxes(la_level)
       bx = layout_get_box(la_level,i)
       do dm = 1,layout_dim(la_level)
          if (bx%lo(dm) == pd%lo(dm)) bc_level(i,dm,1) = domain_bc(dm,1)
          if (bx%hi(dm) == pd%hi(dm)) bc_level(i,dm,2) = domain_bc(dm,2)
       end do
    end do

  end subroutine bc_level_build

end module define_bc_module
