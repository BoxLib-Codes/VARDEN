module regrid_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use box_util_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,uold,sold,gp,p,dx,the_bc_tower)

    use multifab_physbc_module
    use multifab_fill_ghost_module
    use ml_cc_restriction_module
    use make_new_grids_module

    use probin_module, only : verbose, nodal, pmask, &
         amr_buf_width, ref_ratio, max_levs, nscal, nlevs, max_grid_size, &
         ng_cell, ng_grow

    type(ml_layout), intent(inout) :: mla
    type(multifab),  pointer       :: uold(:),sold(:),gp(:),p(:)
    real(dp_t)    ,  pointer       :: dx(:,:)
    type(bc_tower),  intent(inout) :: the_bc_tower

    ! local
    logical           :: new_grid
    integer           :: n, nl, dm, ng_buffer
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba

    ! These are copies to hold the old data.
    type(multifab) :: uold_temp(max_levs), sold_temp(max_levs), gp_temp(max_levs)
    type(multifab) :: p_temp(max_levs)

    dm    = mla%dim

    if (verbose .ge. 1) then
       if (parallel_IOProcessor()) print*,'Calling regrid'
    end if

    if (max_levs < 2) then
       call bl_error('Dont call regrid with max_levs < 2')
    end if

    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n = 1,nlevs

       ! Create copies of the old data.
       call multifab_build(  uold_temp(n),mla_old%la(n),   dm, ng_cell)
       call multifab_build(  sold_temp(n),mla_old%la(n),nscal, ng_cell)
       call multifab_build(    gp_temp(n),mla_old%la(n),   dm, ng_grow)
       call multifab_build(     p_temp(n),mla_old%la(n),    1, ng_grow, nodal)

       ! This is to take care of the ghost cells.
       call setval(p_temp(n),0.d0,all=.true.)

       call multifab_copy_c(  uold_temp(n),1,  uold(n),1,   dm)
       call multifab_copy_c(  sold_temp(n),1,  sold(n),1,nscal)
       call multifab_copy_c(    gp_temp(n),1,    gp(n),1,   dm)
       call multifab_copy_c(     p_temp(n),1,     p(n),1,    1)

       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy(  uold(n))
       call multifab_destroy(  sold(n))
       call multifab_destroy(    gp(n))
       call multifab_destroy(     p(n))

    end do

    call destroy(mla)

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might 
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    do n = 1, max_levs-1
       mba%rr(n,:) = ref_ratio
    enddo

    if (associated(uold)) then
       deallocate(uold,sold,p,gp)
    end if

    allocate(uold(max_levs),sold(max_levs),p(max_levs),gp(max_levs))

    ! Copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! Copy the pd(:)
    mba%pd(1) = mla_old%mba%pd(1)
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

    ! This makes sure the boundary conditions are properly defined everywhere
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    ! Build and fill the level 1 data only.
    call build_and_fill_data(1,la_array(1),mla_old, &
                             uold     ,sold     ,gp     ,p     , &
                             uold_temp,sold_temp,gp_temp,p_temp, &
                             the_bc_tower,dm,mba%rr(1,:))

    nl       = 1
    new_grid = .true.
 
    ng_buffer = 4

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       ! Need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(sold(nl))
       call multifab_physbc(sold(nl),1,dm+1,nscal, &
                            the_bc_tower%bc_tower_array(nl))

       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                           amr_buf_width,ref_ratio,nl,max_grid_size)

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Need to enforce proper nesting within the grid creation procedure 
          !  so that we can fillpatch the new levels.
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, 2, nl+1)) then

                do n = 2,nl
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(  sold(n))
                   call destroy(  uold(n))
                   call destroy(    gp(n))
                   call destroy(     p(n))
                end do

                call enforce_proper_nesting(mba,la_array,max_grid_size)

                ! Loop over all the lower levels which we might have changed when we enforced proper nesting.
                do n = 2,nl
   
                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n))
   
                   ! Rebuild the lower level data again if it changed.
                   call build_and_fill_data(n,la_array(n),mla_old, &
                                            uold     ,sold     ,gp     ,p     , &
                                            uold_temp,sold_temp,gp_temp,p_temp, &
                                            the_bc_tower,dm,mba%rr(n-1,:))
                end do
             end if
          end if

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call build_and_fill_data(nl+1,la_array(nl+1),mla_old, &
                                   uold     ,sold     ,gp     ,p     , &
                                   uold_temp,sold_temp,gp_temp,p_temp, &
                                   the_bc_tower,dm,mba%rr(nl,:))

          nlevs = nl+1
          nl = nl+1

       endif

    enddo

    nlevs = nl

    ! Note: This build actually sets mla%la(n) = la_array(n) so we mustn't delete la_array(n).
    !       Doing the build this way means we don't have to re-create all the multifabs because
    !       we have kept the same layouts.
    call build(mla,mba,la_array,pmask,nlevs)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n = 1, nlevs
       call bc_tower_level_build(the_bc_tower,n,la_array(n))
    end do

    do nl = 1, mla_old%nlevel
       call destroy(  uold_temp(nl))
       call destroy(  sold_temp(nl))
       call destroy(    gp_temp(nl))
       call destroy(     p_temp(nl))
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(uold(nlevs))
       call multifab_fill_boundary(sold(nlevs))
       call multifab_fill_boundary(  gp(nlevs))
       call multifab_fill_boundary(   p(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(uold(nlevs),1,1,dm,the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(sold(nlevs),1,dm+1,nscal, &
                            the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(p(nlevs),1,1,1,the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(gp(nlevs),1,1,dm,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(uold(n-1),uold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(  gp(n-1),  gp(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(uold(n),uold(n-1),ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,dm,fill_crse_input=.false.)
          call multifab_fill_ghost_cells(sold(n),sold(n-1),ng_cell,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,dm+1,nscal,fill_crse_input=.false.)
          call multifab_fill_ghost_cells(gp(n),gp(n-1),ng_grow,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,1,dm,fill_crse_input=.false.)

       enddo

    end if

    call destroy(mba)

    call destroy(mla_old)

  end subroutine regrid

  subroutine build_and_fill_data(lev,la,mla_old, &
                                 uold     ,sold     ,gp     ,p     , &
                                 uold_temp,sold_temp,gp_temp,p_temp, &
                                 the_bc_tower,dm,rr)

    use fillpatch_module
    use ml_prolongation_module
    use probin_module, only : nodal, nscal, ng_cell, ng_grow

    integer                    , intent(in   ) :: lev, dm, rr(:)
    type(layout)               , intent(in   ) :: la
    type(ml_layout)            , intent(in   ) :: mla_old
    type(bc_tower)             , intent(inout) :: the_bc_tower
    type(multifab)   ,  pointer, intent(inout) :: uold(:),sold(:),gp(:),p(:)
    type(multifab)             , intent(in   ) :: uold_temp(:),sold_temp(:),gp_temp(:),p_temp(:)
 
    ! Build the level lev data only.
    call multifab_build(  uold(lev), la,    dm, ng_cell)
    call multifab_build(  sold(lev), la, nscal, ng_cell)
    call multifab_build(    gp(lev), la,    dm, ng_grow)
    call multifab_build(     p(lev), la,     1, ng_grow, nodal)

    ! We need this to take care of the ghost cells.
    call setval(p(lev),0.d0,all=.true.)

    ! Fill the data in the new level lev state -- first from the coarser data if lev > 1.

    if (lev .gt. 1) then

       call fillpatch(uold(lev),uold(lev-1), &
                      ng_cell,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,1,dm)
       call fillpatch(sold(lev),sold(lev-1), &
                      ng_cell,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,dm+1,nscal)
       call fillpatch(gp(lev),gp(lev-1), &
                      ng_grow,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,1,dm)
       ! We interpolate p differently because it is nodal, not cell-centered
       call ml_nodal_prolongation(p(lev), p(lev-1), rr)

    end if

    ! Copy from old data at current level, if it exists
    if (mla_old%nlevel .ge. lev) then
       call multifab_copy_c(  uold(lev),1,  uold_temp(lev),1,   dm)
       call multifab_copy_c(  sold(lev),1,  sold_temp(lev),1,nscal)
       call multifab_copy_c(    gp(lev),1,    gp_temp(lev),1,   dm)
       call multifab_copy_c(     p(lev),1,     p_temp(lev),1,    1)
    end if

  end subroutine build_and_fill_data

end module regrid_module
