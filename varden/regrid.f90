module regrid_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use init_module
  use box_util_module
  use make_new_grids_module
  use initialize_module
  use fillpatch_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla_old,old_u,old_s,old_gp,old_p,mla_new,new_u,new_s,new_gp,new_p,dx,the_bc_tower)

     use probin_module, only : dim_in, nlevs, nscal, ng_cell, ng_grow, nodal, &
                               pmask, regrid_int, max_grid_size, ref_ratio, max_levs, &
                               verbose

     type(ml_layout),intent(in   )   :: mla_old
     type(ml_layout),intent(  out)   :: mla_new
     real(dp_t)    , pointer       :: dx(:,:)
     type(multifab), intent(inout) :: old_u(:),old_s(:),old_gp(:),old_p(:)
     type(multifab), pointer       :: new_u(:),new_s(:),new_gp(:),new_p(:)
     type(bc_tower), intent(inout) :: the_bc_tower

     integer                   :: buf_wid
     type(layout), allocatable :: la_array(:)
     type(boxarray)            :: ba_new,ba_new_comp,ba_old_comp
     type(boxarray)            :: ba_newest
     type(ml_boxarray)         :: mba
     type(list_box)            :: bl

     logical              :: new_grid
     integer              :: i, n, nl, dm

     print *,'DOING REGRID '
     call print(mla_old,'OLD MLA')
   
     if (max_levs < 2) &
       call bl_error('Dont call regrid with max_levs < 2')

     dm = old_u(1)%dim

     buf_wid = regrid_int

     ! mba is big enough to hold max_levs levels
     ! even though we know we had nlevs last time, we might 
     !   want more or fewer levels after regrid (if nlevs < max_levs)
     call ml_boxarray_build_n(mba,max_levs,dm)

     do n = 1, max_levs-1
        mba%rr(n,:) = mla_old%mba%rr(n,:)
     enddo

     allocate(la_array(max_levs))
     allocate(new_u(max_levs),new_s(max_levs),new_p(max_levs),new_gp(max_levs))

     ! Copy the level 1 boxarray
     call copy(mba%bas(1),mla_old%mba%bas(1))

     ! Copy the pd(:)
     do n = 1, max_levs
        mba%pd(n) = mla_old%mba%pd(n)
     enddo

     ! Build the level 1 layout.
     call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

     ! Build the level 1 data only.
     call make_new_state(la_array(1),new_u(1),new_s(1),new_gp(1),new_p(1)) 

     ! Copy the level 1 data.
     call multifab_copy_c(new_u(1) ,1,old_u(1) ,1,   dm)
     call multifab_copy_c(new_s(1) ,1,old_s(1) ,1,nscal)
     call multifab_copy_c(new_gp(1),1,old_gp(1),1,   dm)
     call multifab_copy_c(new_p(1) ,1,old_p(1) ,1,    1)

     new_grid = .true.
     nl = 1

     do while ( (nl .lt. max_levs) .and. (new_grid) )

        ! Do we need finer grids?
        call make_new_grids(la_array(nl),la_array(nl+1),new_s(nl),dx(nl,1),buf_wid,&
                            ref_ratio,nl,max_grid_size,new_grid)
        
        if (new_grid) then

           call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

           ! Build the level nl+1 data only.
           call make_new_state(la_array(nl+1),new_u(nl+1),new_s(nl+1),new_gp(nl+1),new_p(nl+1)) 

           ! Define bc_tower at level nl+1.
           call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
         
           ! Fill the data in the new level nl+1 state -- first from the coarser data.
           call fillpatch(new_u(nl+1),new_u(nl), &
                      ng_cell,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)
            call fillpatch(new_s(nl+1),new_s(nl), &
                      ng_cell,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,dm+1,nscal)
            call fillpatch(new_gp(nl+1),new_gp(nl), &
                      ng_grow,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)
!           call fillpatch(new_p(nl+1),new_p(nl), &
!                     ng_grow,mba%rr(nl,:), &
!                     the_bc_tower%bc_tower_array(nl  ), &
!                     the_bc_tower%bc_tower_array(nl+1), &
!                     1,1,1,1)
            
           ! Copy from old data at current level, if it exists
           if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(new_u(nl+1) ,1,old_u(nl) ,1,   dm)
             call multifab_copy_c(new_s(nl+1) ,1,old_s(nl) ,1,nscal)
             call multifab_copy_c(new_gp(nl+1),1,old_gp(nl),1,   dm)
             call multifab_copy_c(new_p(nl+1) ,1,old_p(nl) ,1,    1)
           end if

           nlevs = nl+1
           nl = nl + 1

        endif ! if (new_grid) 

     enddo ! end while

      do n = 1,nl
         call destroy(new_s(n))
         call destroy(new_u(n))
         call destroy(new_gp(n))
         call destroy(new_p(n))
      end do

      nlevs = nl

      call ml_layout_restricted_build(mla_new,mba,nlevs,pmask)

      ! check for proper nesting
      if (nlevs .ge. 3) then
        
         nl = nlevs - 1
         new_grid = .true.

         do while ( (nl .ge. 2) .and. (new_grid) )

            if (.not. ml_boxarray_properly_nested(mba, ng_cell, pmask, nl, nl+1)) then

                new_grid = .true.
               
                if ( parallel_IOProcessor() .and. verbose .ge. 1) then
                   print *,' '
                   print *,'LEVEL ',nl+1,' grids are not properly nested '
                   print *,' '
                end if

                ! Buffer returns a boxarray "ba_new" that contains everything at level nl 
                !  that the level nl+1 level will need for proper nesting
                call buffer(nl,la_array(nl+1),la_array(nl),la_array(nl-1),&
                            ba_new,ref_ratio,ng_cell)
                call boxarray_intersection(ba_new,mba%pd(nl))

                ! Make sure the new grids start on even and end on odd
                call boxarray_coarsen(ba_new, ref_ratio)
                call boxarray_refine (ba_new, ref_ratio)

                ! Merge the new boxarray "ba_new" with the existing box_array 
                ! mba%bas(nl) so that we get the union of points.
                call boxarray_complementIn(ba_old_comp,mba%pd(nl),mba%bas(nl))
                call boxarray_intersection(ba_old_comp,ba_new)
                do i = 1, mba%bas(nl)%nboxes
                   call push_back(bl,  mba%bas(nl)%bxs(i))
                end do
                do i = 1, ba_old_comp%nboxes
                   call push_back(bl, ba_old_comp%bxs(i))
                end do
                call build(ba_newest,bl)
                call destroy(bl)
                call boxarray_simplify(ba_newest)
                call boxarray_maxsize(ba_newest,max_grid_size)

                ! Replace mba%bas(nl) by ba_new
                call destroy(mba%bas(nl))
                call copy(mba%bas(nl),ba_newest)
                call destroy(ba_newest)

                ! Double check we got the proper nesting right
                if (.not. ml_boxarray_properly_nested(mba, ng_cell, pmask, nl, nl+1)) &
!                 call bl_error('Still not properly nested, darn it')
                  print *,'WARNING WARNING -- NOT PROPERLY NESTED'

                ! Destroy the old layout and build a new one.
                call destroy(la_array(nl))
                call layout_build_ba(la_array(nl),mba%bas(nl),mba%pd(nl),pmask)

            else

                new_grid = .false.

            endif  !if not properly nested

            nl = nl - 1

         enddo ! do while
      end if ! if (nlevs .ge. 3)

   call ml_layout_restricted_build(mla_new,mba,nlevs,pmask)

   nlevs = mla_new%nlevel

   ! Now make the data for the final time.
   do nl = 1,nlevs-1

      call make_new_state(mla_new%la(nl+1),new_u(nl+1),new_s(nl+1),new_gp(nl+1),new_p(nl+1)) 

      ! Fill the data in the new level nl+1 state -- first from the coarser data.
       call fillpatch(new_u(nl+1),old_u(nl), &
                      ng_cell,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)
        call fillpatch(new_s(nl+1),old_s(nl), &
                      ng_cell,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,dm+1,nscal)
        call fillpatch(new_gp(nl+1),old_gp(nl), &
                      ng_grow,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)
!       call fillpatch(new_p(nl+1),old_p(nl), &
!                     ng_grow,mba%rr(nl,:), &
!                     the_bc_tower%bc_tower_array(nl  ), &
!                     the_bc_tower%bc_tower_array(nl+1), &
!                     1,1,1,1)
 
        ! Copy from old data at current level, if it exists
        if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(new_u(nl+1) ,1,old_u(nl) ,1,   dm)
             call multifab_copy_c(new_s(nl+1) ,1,old_s(nl) ,1,nscal)
             call multifab_copy_c(new_gp(nl+1),1,old_gp(nl),1,   dm)
             call multifab_copy_c(new_p(nl+1) ,1,old_p(nl) ,1,    1)
        end if

   end do

   call destroy(mba)

   do n = 1,nl
      call destroy(old_s(n))
      call destroy(old_u(n))
      call destroy(old_gp(n))
      call destroy(old_p(n))
   end do

   call print(mla_new,'NEW MLA')

  end subroutine regrid

end module regrid_module
