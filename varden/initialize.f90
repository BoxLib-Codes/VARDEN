module initialize_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use init_module
  use box_util_module
  use make_new_grids_module

  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, initialize_with_adaptive_grids, &
            make_new_state

contains

  subroutine initialize_from_restart(mla,restart,time,dt,dx,pmask,uold,sold,gp,p,the_bc_tower)
 
     use probin_module, only : dim_in, nlevs, nscal, ng_cell, ng_grow, nodal

     type(ml_layout),intent(out)   :: mla
     integer       , intent(in   ) :: restart
     real(dp_t)    , intent(  out) :: time,dt
     real(dp_t)    , pointer       :: dx(:,:)
     logical       , intent(in   ) :: pmask(:)
     type(multifab), pointer       :: sold(:),uold(:),gp(:),p(:)
     type(bc_tower), intent(  out) :: the_bc_tower

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(multifab), pointer   :: chk_p(:)

     integer :: n,dm

     dm = dim_in

     call fill_restart_data(restart,mba,chkdata,chk_p,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     allocate(uold(nlevs),sold(nlevs),gp(nlevs),p(nlevs))

     do n = 1,nlevs
        call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
        call multifab_build(   sold(n), mla%la(n), nscal, ng_cell)
        call multifab_build(     gp(n), mla%la(n),    dm, ng_grow)
        call multifab_build(      p(n), mla%la(n),     1, ng_grow, nodal)
     end do
     do n = 1,nlevs
        call multifab_copy_c(uold(n),1,chkdata(n),1         ,dm)
        call multifab_copy_c(sold(n),1,chkdata(n),1+dm      ,nscal)
        call multifab_copy_c(  gp(n),1,chkdata(n),1+dm+nscal,dm)
        call multifab_copy_c(   p(n),1,  chk_p(n),1         ,1)
        !
        ! The layouts for chkdata and chk_p are built standalone, level
        ! by level, and need to be destroy()d as such as well.
        !
        call destroy(chkdata(n)%la)
        call destroy(chk_p(n)%la)
        call multifab_destroy(chkdata(n))
        call multifab_destroy(chk_p(n))
     end do
     deallocate(chkdata,chk_p)

     call initialize_dx(dx,mba,nlevs)

     call initialize_bc(the_bc_tower,nlevs,pmask)
     do n = 1,nlevs
        call bc_tower_level_build( the_bc_tower,n,mla%la(n))
     end do

  end subroutine initialize_from_restart

  subroutine initialize_with_fixed_grids(mla,pmask,dx,uold,sold,gp,p,the_bc_tower)

     use probin_module, only : dim_in, nlevs, nscal, ng_cell, ng_grow, nodal, fixed_grids

     type(ml_layout),intent(out)   :: mla
     logical       , intent(in   ) :: pmask(:)
     real(dp_t)    , pointer       :: dx(:,:)
     type(multifab), pointer       :: uold(:),sold(:),gp(:),p(:)
     type(bc_tower), intent(  out) :: the_bc_tower

     type(ml_boxarray)         :: mba

     integer :: n,dm

     dm = dim_in

     call read_a_hgproj_grid(mba, fixed_grids)
     call ml_layout_build(mla,mba,pmask)

     ! check for proper nesting
     if (.not. ml_boxarray_properly_nested(mla%mba, ng_cell, pmask)) &
         call bl_error('fixed_grids not properly nested')

     nlevs = mla%nlevel
     allocate(uold(nlevs),sold(nlevs),p(nlevs),gp(nlevs))

     do n = 1,nlevs
        call multifab_build(   uold(n), mla%la(n),    dm, ng_cell)
        call multifab_build(   sold(n), mla%la(n), nscal, ng_cell)
        call multifab_build(     gp(n), mla%la(n),    dm, ng_grow)
        call multifab_build(      p(n), mla%la(n),     1, ng_grow, nodal)
     end do

     call initialize_dx(dx,mba,nlevs)

     call initialize_bc(the_bc_tower,nlevs,pmask)
     do n = 1,nlevs
        call bc_tower_level_build( the_bc_tower,n,mla%la(n))
     end do

     call initdata(nlevs,uold,sold,dx,the_bc_tower%bc_tower_array,mla)

  end subroutine initialize_with_fixed_grids

  subroutine initialize_with_adaptive_grids(mla,pmask,dx,uold,sold,gp,p,the_bc_tower)

     use probin_module, only : dim_in, nlevs, nscal, ng_cell, ng_grow, nodal, &
                               n_cellx, n_celly, n_cellz, &
                               regrid_int, max_grid_size, ref_ratio, max_levs, &
                               verbose

     type(ml_layout),intent(out)   :: mla
     logical       , intent(in   ) :: pmask(:)
     real(dp_t)    , pointer       :: dx(:,:)
     type(multifab), pointer       :: uold(:),sold(:),gp(:),p(:)
     type(bc_tower), intent(  out) :: the_bc_tower

     integer                   :: buf_wid
     type(layout), allocatable :: la_array(:)
     type(box)   , allocatable :: bxs(:)
     type(boxarray)            :: ba_new,ba_new_comp,ba_old_comp
     type(boxarray)            :: ba_newest
     type(ml_boxarray)         :: mba
     type(list_box)            :: bl

     logical              :: new_grid
     integer, allocatable :: lo(:), hi(:)
     integer              :: i, n, nl, dm

     dm = dim_in

!    buf_wid = regrid_int
     buf_wid = 0

     ! set up hi & lo to carry indexing info
     allocate(lo(dm),hi(dm))
     lo(:) = 0
     hi(1) = n_cellx-1
     if (dm > 1) then   
        hi(2) = n_celly - 1        
        if (dm > 2)  then
           hi(3) = n_cellz -1
        endif
     endif

     ! mba is big enough to hold max_levs levels
     call ml_boxarray_build_n(mba,max_levs,dm)
     do n = 1, max_levs-1
        mba%rr(n,:) = ref_ratio
     enddo

     if (max_levs > 1) allocate(la_array(max_levs))
     allocate(bxs(max_levs))
     allocate(uold(max_levs),sold(max_levs),p(max_levs),gp(max_levs))

       ! Build the level 1 boxarray
     call box_build_2(bxs(1),lo,hi)
     call boxarray_build_bx(mba%bas(1),bxs(1))
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! build pd(:)
     mba%pd(1) = bxs(1)
     do n = 2, max_levs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     enddo

     ! Need to build pd before making dx
     call initialize_dx(dx,mba,max_levs)

     ! Initialize bc's.
     call initialize_bc(the_bc_tower,max_levs,pmask)

     if (max_levs > 1) then

        ! Build the level 1 layout.
        call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

        ! Build the level 1 data only.
        call make_new_state(la_array(1),uold(1),sold(1),gp(1),p(1)) 

        ! Define bc_tower at level 1.
        call bc_tower_level_build(the_bc_tower,1,la_array(1))

        ! Initialize the level 1 data only.
        call initdata_on_level(uold(1),sold(1),dx(1,:),the_bc_tower%bc_tower_array(1),la_array(1))

        new_grid = .true.
        nl = 1

        do while ( (nl .lt. max_levs) .and. (new_grid) )

           ! Do we need finer grids?
           call make_new_grids(la_array(nl),la_array(nl+1),sold(nl),dx(nl,1),buf_wid,&
                               ref_ratio,nl,max_grid_size,new_grid)
        
           if (new_grid) then

              call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

              ! Build the level nl+1 data only.
              call make_new_state(la_array(nl+1),uold(nl+1),sold(nl+1),gp(nl+1),p(nl+1)) 

              ! Define bc_tower at level nl+1.
              call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
            
             ! fills the physical region of each level with problem data (blob now)
              call initdata_on_level(uold(nl+1),sold(nl+1),dx(nl+1,:),the_bc_tower%bc_tower_array(nl+1),la_array(nl+1))

              nlevs = nl+1
              nl = nl + 1

           endif ! if (new_grid) 

      enddo          

      do n = 1,nl
         call destroy(sold(n))
         call destroy(uold(n))
         call destroy(gp(n))
         call destroy(p(n))
      end do

      nlevs = nl

      ! check for proper nesting
      if (nlevs .ge. 3) then
        
         nl = nlevs - 1
         new_grid = .true.

         do while ( (nl .ge. 2) .and. (new_grid) )

            if (.not. ml_boxarray_properly_nested(mba, ng_cell, pmask, nl, nl+1)) then

                new_grid = .true.
               
                ! Buffer returns a boxarray "ba_new" that contains everything at level nl 
                !  that the level nl+1 level will need for proper nesting
                call buffer(nl,la_array(nl+1),ba_new,ref_ratio,ng_cell)
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
                  print *,'WARNING WARNING -- LEVEL ',nl+1,' GRIDS ARE NOT PROPERLY NESTED'

                ! Destroy the old layout and build a new one.
                call destroy(la_array(nl))
                call layout_build_ba(la_array(nl),mba%bas(nl),mba%pd(nl),pmask)

            else

                new_grid = .false.

            endif  !if not properly nested

            nl = nl - 1

         enddo ! do while
      end if ! if (nlevs .ge. 3)

   end if ! end if (maxlev > 1)

   call ml_layout_restricted_build(mla,mba,nlevs,pmask)

   nlevs = mla%nlevel

   do n = 1,nlevs
      call make_new_state(mla%la(n),uold(n),sold(n),gp(n),p(n)) 
   end do

   call initdata(nlevs,uold,sold,dx,the_bc_tower%bc_tower_array,mla)

   call destroy(mba)
   deallocate(bxs)

  end subroutine initialize_with_adaptive_grids

  subroutine make_new_state(la_loc,uold_loc,sold_loc,gp_loc,p_loc)

    use probin_module, only : dim_in, nscal, ng_cell, ng_grow, nodal
    use bl_constants_module

    type(layout),intent(in   ) :: la_loc
    type(multifab ),intent(inout) :: uold_loc,sold_loc,gp_loc,p_loc
 
    integer :: dm

    dm = dim_in

    call multifab_build(   uold_loc, la_loc,    dm, ng_cell)
    call multifab_build(   sold_loc, la_loc, nscal, ng_cell)
    call multifab_build(     gp_loc, la_loc,    dm, ng_grow)
    call multifab_build(      p_loc, la_loc,     1, ng_grow, nodal)

    call setval(  uold_loc,ZERO, all=.true.)
    call setval(  sold_loc,ZERO, all=.true.)
    call setval(    gp_loc,ZERO, all=.true.)
    call setval(     p_loc,ZERO, all=.true.)

  end subroutine make_new_state

  subroutine initialize_bc(the_bc_tower,num_levs,pmask)

     use bc_module

     use probin_module, only : dim_in, &
                               bcx_lo, bcx_hi, bcy_lo, bcy_hi, bcz_lo, bcz_hi

     type(bc_tower), intent(  out) :: the_bc_tower
     integer       , intent(in   ) :: num_levs
     logical       , intent(in   ) :: pmask(:)

     integer, allocatable :: domain_phys_bc(:,:)

     integer :: dm

     dm = dim_in

     ! Define the physical boundary conditions on the domain
     allocate(domain_phys_bc(dm,2))
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bcx_lo
     domain_phys_bc(1,2) = bcx_hi
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bcy_lo
        domain_phys_bc(2,2) = bcy_hi
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bcz_lo
        domain_phys_bc(3,2) = bcz_hi
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.

     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

     deallocate(domain_phys_bc)

  end subroutine initialize_bc

  subroutine initialize_dx(dx,mba,num_levs)

     use probin_module, only : dim_in, prob_hi
  
     real(dp_t)       , pointer     :: dx(:,:)
     type(ml_boxarray), intent(in ) :: mba
     integer          , intent(in ) :: num_levs

     integer :: i,n,dm

     dm = dim_in

     allocate(dx(num_levs,dm))

     do i = 1,dm
        dx(1,i) = prob_hi(i) / float(mba%pd(1)%hi(i)-mba%pd(1)%lo(i)+1)
     end do
     do n = 2,num_levs
        dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
     end do

  end subroutine initialize_dx

end module initialize_module
