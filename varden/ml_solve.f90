module ml_solve_module

   use bl_types
   use define_bc_module
   use multifab_module
   use boxarray_module
   use stencil_module
   use mg_module
   use list_box_module
   use mboxarray_module
   use itsol_module
   use sparse_solve_module
   use bl_mem_stat_module
   use bl_timer_module
   use box_util_module
   use bl_IO_module
   use fabio_module

   use ml_restriction_module
   use ml_prolongation_module
   use ml_interface_stencil_module
   use ml_util_module
   use bndry_reg_module
   use flux_reg_module
 
   implicit none

   real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
   real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t
   real (kind = dp_t), private, parameter :: ONE   = 1.0_dp_t

contains

   subroutine ml_cc_solve(la_tower,mgt,rh,soln,fine_flx,bc,stencil_order)

      type(layout  ), intent(inout) :: la_tower(:)
      type(mg_tower), intent(inout) :: mgt(:)
      type(multifab), intent(inout) :: rh(:)
      type(multifab), intent(inout) :: soln(:)
      type(flux_reg), intent(inout) :: fine_flx(2:)
      integer       , intent(in   ) :: bc(:,:)
      integer       , intent(in   ) :: stencil_order

      type(boxarray) :: pdv
      type(box     ) :: pd,pdc,pdf
      type(boxarray) :: bac
      
      type(multifab), pointer :: uu(:)      => Null()
      type(multifab), pointer :: uu_hold(:) => Null()
      type(multifab), pointer :: res(:)     => Null()
      type(multifab), pointer :: temp_res(:)=> Null()
      type(lmultifab), pointer :: fine_mask(:) => Null()

      type(bndry_reg), pointer ::  brs_flx(:)   => Null()
      type(bndry_reg), pointer ::  brs_bcs(:)   => Null()

      type(layout) :: la
      integer, allocatable :: lo(:), hi(:)

      integer  :: i,n,dm,ns
      integer  :: mglev, mglev_crse, iter, it
      integer  :: nlevs
      integer  :: do_diagnostics

      integer, allocatable :: ref_ratio(:),hi_fine(:),hi_crse(:)

      real(dp_t) :: Anorm, bnorm, eps
      real(dp_t) :: res_norm,soln_norm
      real(dp_t) :: snrm(2)

      dm = layout_dim(la_tower(1))

      allocate(lo(dm), hi(dm))
      nlevs = size(la_tower)

      do_diagnostics = 0

      eps = 1.d-12

      print *, 'NLEVS ', nlevs
      print *, 'EPS   ', eps

      allocate(uu(nlevs), uu_hold(nlevs-1), res(nlevs))
      allocate(temp_res(nlevs))
      allocate(fine_mask(nlevs))

      allocate(brs_flx(2:nlevs))
      allocate(brs_bcs(2:nlevs))

      allocate(ref_ratio(dm))
      allocate(hi_fine(dm), hi_crse(dm))

  ! NOTE THIS CHANGE: we now have stencil values which reach outside the
  !  grid in the case of Dirichlet bc's and skewed stencils
  ns = 1 + dm*3 

  do n = nlevs, 1, -1

     pd = layout_get_pd(la_tower(n))
     print *,'LEVEL ',n
     print *,'PD ', (extent(pd,dim=i), i=1,dm)
     la = la_tower(n)
     call multifab_build(uu(n), la, 1, 1)
     if ( n < nlevs ) then
        call multifab_build(uu_hold(n), la, 1, 1)
        call setval(uu_hold(n), ZERO, all=.true.)
     end if
     call multifab_build(res(n), la, 1, 0)
     call multifab_build(temp_res(n), la, 1, 0)
     call lmultifab_build(fine_mask(n), la, 1, 0)
     call setval(uu(n), ZERO, all=.true.)
     call multifab_copy(res(n),rh(n),all=.true.)

     call setval(fine_mask(n), val = .true., all = .true.)
     if ( n < nlevs ) then
        hi_fine = upb(layout_get_pd(la_tower(n+1))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n  ))) + 1
        ref_ratio = hi_fine / hi_crse
        call copy(bac, get_boxarray(la_tower(n+1)))
        call boxarray_coarsen(bac, ref_ratio)
        call setval(fine_mask(n), .false., bac)
        call destroy(bac)
     end if

     if ( n == 1 ) exit

     ! Build the (coarse resolution) flux registers to be used in computing
     !  the residual at a non-finest AMR level.

     hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
     hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
     ref_ratio = hi_fine / hi_crse
     pdc = layout_get_pd(la_tower(n-1))
     pdf = layout_get_pd(la_tower(n  ))
     call bndry_reg_build(  brs_flx(n), la, ref_ratio, pdc, width = 0)
     call bndry_reg_build(  brs_bcs(n), la, ref_ratio, pdc, width = 2)
  end do

  !!
  !! Solver Starts Here
  !!
  
  mglev = mgt(nlevs)%nlevels

  Anorm = stencil_norm(mgt(nlevs)%ss(mglev))
  bnorm = norm_inf(rh(nlevs))

  do n = 1, nlevs-1
     bnorm = max(norm_inf(rh(n), fine_mask(n)), bnorm)
     mglev = mgt(n)%nlevels
     Anorm = max(stencil_norm(mgt(n)%ss(mglev), fine_mask(n)), Anorm)
  end do

  print *, 'bNorm = ', bNorm
  print *, 'ANorm = ', ANorm

  do iter = 1, mgt(nlevs)%max_iter

     if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps) ) exit

     ! Initialize corrections at all levels to zero.
     do n = 1,nlevs
        call setval(uu(n)     , ZERO, all=.true.)
     end do

     do n = 1,nlevs-1
        call setval(uu_hold(n), ZERO, all=.true.)
     end do

     !   Down the V-cycle
     do n = nlevs,1,-1

        mglev = mgt(n)%nlevels

        if (do_diagnostics == 1) then
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
           print *,'DWN: BEFORE GSRB: RES AT LEVEL ',n, norm_inf(temp_res(n))
        end if

        call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
             uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
             mgt(n)%gamma)

        if (n > 1) then
           hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
           hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
           ref_ratio = hi_fine / hi_crse
           mglev_crse = mgt(n-1)%nlevels

           ! There are three steps in defining the RHS at the non-finest
           ! AMR level:
           ! (1) Use the existing residual away from the finer grids
           ! (already defined)
           call mg_defect(mgt(n-1)%ss(mglev_crse), res(n-1), &
                          rh(n-1), soln(n-1), mgt(n-1)%mm(mglev_crse))

           pdc = layout_get_pd(la_tower(n-1))
           call crse_fine_residual(n,mgt,uu,brs_flx(n),pdc,ref_ratio)

           ! (3) Restrict the residual from the finer level (impt to do
           !     this last so we overwrite anything extra which may have
           !     been defined by (2) near fine-fine interfaces)
           mglev = mgt(n)%nlevels
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))

           call multifab_copy(res(n), temp_res(n), all=.true.)
           mglev_crse = mgt(n-1)%nlevels
           call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev), &
                mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)

           if (n < nlevs) call multifab_copy(uu_hold(n), uu(n), all=.true.)

           call saxpy(soln(n), ONE, uu(n))
           call setval(uu(n), ZERO, all=.true.)

        else
           call saxpy(soln(n), ONE, uu(n))
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
           call multifab_copy(res(n), temp_res(n), all=.true.)
           if (nlevs == 1) &
                call mg_defect(mgt(n)%ss(mglev), res(n), rh(n), soln(n), mgt(n)%mm(mglev))
        end if

        if (do_diagnostics == 1) then
           print *,'DWN: AFTER  GSRB: RES AT LEVEL ', n, norm_inf(temp_res(n))
        end if

     end do

     !   Back up the V-cycle
     do n = 2, nlevs

        pd = layout_get_pd(la_tower(n))
        hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
        ref_ratio = hi_fine / hi_crse
        call ml_prolongation(uu(n), uu(n-1), pd, ref_ratio)

        call bndry_reg_copy(brs_bcs(n), uu(n-1))
        do i = 1, dm
           call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio, -i)
           call ml_interp_bcs(uu(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio, +i)
        end do

        mglev = mgt(n)%nlevels

        if (do_diagnostics == 1) then
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
           print *,'UP: BEFORE GSRB: RES AT LEVEL ',n, norm_inf(temp_res(n))
        end if

        call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
             uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
             mgt(n)%gamma)

        call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
        call multifab_copy(res(n), temp_res(n), all=.true.)

        if (do_diagnostics == 1) &
           print *,'UP: AFTER  GSRB: RES AT LEVEL ',n, norm_inf(temp_res(n))

        call saxpy(soln(n), ONE, uu(n), all = .true.)
        if (n < nlevs) call saxpy(uu(n), ONE, uu_hold(n), all = .true.)

        ! Only do this as long as tangential interp looks under fine grids
        mglev_crse = mgt(n-1)%nlevels
        call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
             mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)

     end do

     do n = nlevs,2,-1
        hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
        ref_ratio = hi_fine / hi_crse
!       call mg_defect(mgt(n)%ss(mglev), temp_res(n), res(n), uu(n), mgt(n)%mm(mglev))
!       call multifab_copy(res(n), temp_res(n), all=.true.)
        call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),& 
                            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)
        pdc = layout_get_pd(la_tower(n-1))
        call crse_fine_residual(n,mgt,uu,brs_flx(n),pdc,ref_ratio)
     end do

     res_norm = ZERO
     do n = 1,nlevs
        res_norm = max(res_norm,norm_inf(res(n)))
     end do

     soln_norm = ZERO
     do n = 1,nlevs
        soln_norm = max(soln_norm,norm_inf(soln(n)))
     end do

     if ( mgt(nlevs)%verbose > 0 .and. parallel_IOProcessor() ) then
        write(unit=*, fmt='(i3,": Ninf(defect) = ",g15.8)') iter, res_norm
     end if

  end do

  ! Make sure even the coarse cells under fine cells have the best answer.
  do n = nlevs,2,-1
    mglev      = mgt(n  )%nlevels
    mglev_crse = mgt(n-1)%nlevels
    hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
    hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
    ref_ratio = hi_fine / hi_crse
    call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
         mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)
  end do 

! Interpolate boundary conditions of soln in order to get correct grad(phi) at 
!   crse-fine boundaries
  do n = 2,nlevs
     hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
     hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
     ref_ratio = hi_fine / hi_crse
     call bndry_reg_copy(brs_bcs(n), soln(n-1))
     do i = 1, dm
        call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,0), pd, ref_ratio, -i)
        call ml_interp_bcs(soln(n), brs_bcs(n)%bmf(i,1), pd, ref_ratio, +i)
     end do 
     mglev = mgt(n)%nlevels
     do i = 1, dm
        call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,0), &
             soln(n), mgt(n)%mm(mglev), -1, i)
        call ml_fill_fine_fluxes(mgt(n)%ss(mglev), fine_flx(n)%bmf(i,1), &
             soln(n), mgt(n)%mm(mglev),  1, i)
     end do
  end do

  do n = 1,nlevs
    call multifab_fill_boundary(soln(n))
  end do 

  if ( mgt(nlevs)%verbose > 0 .AND. parallel_IOProcessor() ) then
     write(unit=*, fmt='("MG finished at ", i3, " iterations")') iter
  end  if

  snrm(1) = ml_norm_l2 (soln, fine_mask)
  snrm(2) = ml_norm_inf(soln, fine_mask)
  if ( parallel_IOProcessor() ) &
       write(unit=*, fmt='("solution norm = ",es24.16, "/",es24.16, " iter = ", i0)') snrm, iter

      do n = 1,nlevs
         call multifab_destroy(uu(n))
         if ( n < nlevs) then
            call multifab_destroy(uu_hold(n))
         end if
         call multifab_destroy(res(n))
         call multifab_destroy(temp_res(n))
         call lmultifab_destroy(fine_mask(n))
         if ( n > 1 ) then
            call bndry_reg_destroy(brs_flx(n))
            call bndry_reg_destroy(brs_bcs(n))
         end if
      end do

      deallocate(uu, uu_hold, res, temp_res)
      deallocate(fine_mask, brs_flx, brs_bcs)

   contains

     subroutine crse_fine_residual(n,mgt,uu,brs_flx,pdc,ref_ratio)

        integer        , intent(in   ) :: n
        type(mg_tower) , intent(inout) :: mgt(:)
        type(bndry_reg), intent(inout) :: brs_flx
        type(multifab) , intent(inout) :: uu(:)
        type(box)      , intent(in   ) :: pdc
        integer        , intent(in   ) :: ref_ratio(:)

        integer :: i,dm,mglev
 
        dm = brs_flx%dim
        mglev = mgt(n)%nlevels

        ! (2) Use a coarse-fine stencil at crse cells adjacent to fine grids
        do i = 1, dm
           call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,0), &
                uu(n), mgt(n)%mm(mglev), ref_ratio(i), -1, i)
           call ml_interface(res(n-1), brs_flx%bmf(i,0), uu(n-1), &
                mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, -1, i)

           call ml_fill_fluxes(mgt(n)%ss(mglev), brs_flx%bmf(i,1), &
                uu(n), mgt(n)%mm(mglev), ref_ratio(i), 1, i)
           call ml_interface(res(n-1), brs_flx%bmf(i,1), uu(n-1), &
                mgt(n-1)%ss(mgt(n-1)%nlevels), pdc, +1, i)
        end do

     end subroutine crse_fine_residual

   end subroutine ml_cc_solve

   subroutine ml_nd_solve(la_tower,mgt,rh,soln)

      type(layout)   , intent(inout) :: la_tower(:)
      type(mg_tower) , intent(inout) :: mgt(:)
      type(multifab) , intent(inout) :: rh(:)
      type(multifab) , intent(inout) :: soln(:)

       type(box)      :: pd,pdc
       type(boxarray) :: pdv

       type(boxarray) :: bac


       integer :: nlevs
       type(multifab), pointer ::       uu(:) => Null()
       type(multifab), pointer ::  uu_hold(:) => Null()
       type(multifab), pointer ::      res(:) => Null()
       type(multifab), pointer :: temp_res(:) => Null()
       type(lmultifab), pointer :: fine_mask(:) => Null()

       type(bndry_reg), pointer :: brs_flx(:) => Null()

       type(layout) :: la
            integer i, n, dm, ns
       integer :: ng_for_res
       integer mglev, mglev_crse, iter, it
       integer, allocatable :: lo(:), hi(:)

       real(dp_t) :: Anorm, bnorm
       real(dp_t) :: eps

       integer, allocatable :: ref_ratio(:)
       integer, allocatable :: hi_fine(:)
       integer, allocatable :: hi_crse(:)
       real(dp_t) :: snrm(2), ni
       real(dp_t) :: res_norm,soln_norm

       logical :: all_done
       logical, allocatable :: nodal(:)

       integer              :: do_diagnostics 

       do_diagnostics = 1

       eps = 1.d-12

       nlevs = size(la_tower)
       dm = rh(nlevs)%dim
       allocate(nodal(dm), lo(dm), hi(dm))
       nodal = .true.

       allocate(uu(nlevs), uu_hold(nlevs), res(nlevs))
       allocate(temp_res(nlevs))
       allocate(fine_mask(nlevs))

       allocate(brs_flx(2:nlevs))

       allocate(ref_ratio(dm))
       allocate(hi_fine(dm))
       allocate(hi_crse(dm))

!      We are only considering the dense stencils here (3 in 1d, 9 in 2d, 27 in 3d)
       ns = 3**dm

!      We carry ghost cells in the residual for the nodal case only.
       ng_for_res = 1

       do n = nlevs, 1, -1

          la = la_tower(n)
          pd = layout_get_pd(la)
          if ( parallel_ioprocessor() ) then
             print *,'LEVEL ',n
             print *,'PD ',extent(pd,dim=1),extent(pd,dim=2)
          end if

          call multifab_build(uu(n), la, 1, 1, nodal)
          call multifab_build(uu_hold(n), la, 1, 1, nodal)
          call multifab_build(res(n), la, 1, ng_for_res, nodal)
          call multifab_build(temp_res(n), la, 1, ng_for_res, nodal)
          call lmultifab_build(fine_mask(n), la, 1, 0, nodal)
          call setval(  uu(n), ZERO,all=.true.)
          call setval(  uu_hold(n), ZERO,all=.true.)
          call setval( res(n), ZERO,all=.true.)
          call setval(temp_res(n), ZERO,all=.true.)
 
          call multifab_copy(res(n),rh(n),all=.true.)

          if ( n < nlevs ) then
             call create_nodal_mask(n,fine_mask(n), &
                                    mgt(n  )%mm(mgt(n  )%nlevels), &
                                    mgt(n+1)%mm(mgt(n+1)%nlevels), &
                                    la_tower)
          else
             call setval(fine_mask(n), val = .true., all = .true.)
          endif
          if ( n == 1 ) exit

!        Build the (coarse resolution) flux registers to be used in computing
!         the residual at a non-finest AMR level.
 
          pdc = layout_get_pd(la_tower(n-1))
          hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
          hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
          ref_ratio = hi_fine / hi_crse
          call bndry_reg_build(brs_flx(n), la, ref_ratio, pdc, nodal = nodal)

       end do

  !!
  !! Solver Starts Here
  !!
  
  n = nlevs
  mglev = mgt(n)%nlevels
  Anorm = stencil_norm(mgt(n)%ss(mglev))

  bnorm = norm_inf(rh(nlevs))
  do n = 1, nlevs-1
     bnorm = max(norm_inf(rh(n), fine_mask(n)), bnorm)
     mglev = mgt(n)%nlevels
     Anorm = max(stencil_norm(mgt(n)%ss(mglev), fine_mask(n)), Anorm)
  end do

  do iter = 1, mgt(nlevs)%max_iter

     if ( ml_converged(res, soln, fine_mask, bnorm, Anorm, eps) ) exit

     ! Initialize corrections at all levels to zero.
     do n = 1,nlevs
      call setval(uu(n), ZERO, all=.true.)
      call setval(uu_hold(n), ZERO, all=.true.)
     end do

     !   Down the V-cycle
     do n = nlevs,1,-1

        mglev = mgt(n)%nlevels

        if (do_diagnostics == 1) then
           call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
                res(n),uu(n),mgt(n)%mm(mglev))
           ni = norm_inf(temp_res(n))
           if ( parallel_ioprocessor() ) then
              print *,'DWN: BEFORE GSRB: RES AT LEVEL ',n, ni
           end if
        end if

        if (n > 1) then
          call mini_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
               mgt(n)%gamma)
        else 
          call mg_tower_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
               uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
               mgt(n)%gamma)
        end if

        if (n > 1) then
           hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
           hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
           ref_ratio = hi_fine / hi_crse
           mglev_crse = mgt(n-1)%nlevels

           ! There are three steps in defining the RHS at the non-finest AMR level:
           ! (1) Use the existing residual away from the finer grids (already defined)
           call mg_defect(mgt(n-1)%ss(mglev_crse),res(n-1), &
                          rh(n-1),soln(n-1),mgt(n-1)%mm(mglev_crse))

           call saxpy(soln(n),ONE,uu(n))

           ! (2) Restrict the residual from the finer level.
           mglev = mgt(n)%nlevels
           call mg_defect(mgt(n)%ss(mglev), temp_res(n), &
                          res(n),uu(n),mgt(n)%mm(mglev))
           call multifab_copy(res(n),temp_res(n),all=.true.)
           call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),& 
                               mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)

           ! (3) Compute a coarse-fine residual at the coarse-fine interface.

           pdc = layout_get_pd(la_tower(n-1))
           call crse_fine_residual(n,mgt,brs_flx(n),temp_res(n),pdc)

           if (n < nlevs) call multifab_copy(uu_hold(n),uu(n),all=.true.)
           call setval(uu(n),ZERO,all=.true.)

        else

           call saxpy(soln(n), ONE, uu(n))

        end if

        call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
                       res(n),uu(n),mgt(n)%mm(mglev))
        call multifab_copy(res(n),temp_res(n),all=.true.)

        if (do_diagnostics == 1 .and.  parallel_ioprocessor() ) &
           print *,'DWN: AFTER  GSRB: RES AT LEVEL ',n, norm_inf(res(n))

     end do

     !   Back up the V-cycle
     do n = 2, nlevs

        pd = layout_get_pd(la_tower(n))
        hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
        ref_ratio = hi_fine / hi_crse
        mglev = mgt(n)%nlevels

        ! Interpolate the correction and add to "uu", "uu_hold" and "soln".
        call ml_prolongation(uu(n), uu(n-1), pd, ref_ratio)
        call saxpy(soln(n), ONE, uu(n), .true.)
        call saxpy(uu_hold (n), ONE, uu(n), .true.)

        ! Calculate the new residual and set "uu" to zero.
        call mg_defect(mgt(n)%ss(mglev),temp_res(n), &
             res(n),uu(n),mgt(n)%mm(mglev))
        call multifab_copy(res(n),temp_res(n),all=.true.)
        call setval(uu(n),ZERO,all=.true.)

        ! Relax ...
        if ( do_diagnostics == 1 .and. parallel_ioprocessor() ) then
           ni = norm_inf(res(n))
           print *,'UP: BEFORE GSRB: RES AT LEVEL ',n, ni
        end if
        call mini_cycle(mgt(n), mgt(n)%cycle, mglev, mgt(n)%ss(mglev), &
             uu(n), res(n), mgt(n)%mm(mglev), mgt(n)%nu1, mgt(n)%nu2, &
             mgt(n)%gamma)

        ! Add the new correction to "uu_hold" and "soln"
        call saxpy(soln(n), ONE,      uu(n), .true.)
        call saxpy(  uu(n), ONE, uu_hold(n), .true.)

        call mg_defect(mgt(n)%ss(mglev),temp_res(n),rh(n),soln(n), &
                       mgt(n)%mm(mglev))
        call multifab_copy(res(n),temp_res(n),all=.true.)

        if ( do_diagnostics == 1 .and. parallel_ioprocessor() ) &
           print *,'UP: AFTER  GSRB: RES AT LEVEL ',n, norm_inf(res(n))

     end do

     do n = nlevs,2,-1
        hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
        ref_ratio = hi_fine / hi_crse
        call ml_restriction(res(n-1), res(n), mgt(n)%mm(mglev),& 
                            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, ref_ratio)
        pdc = layout_get_pd(la_tower(n-1))
        call crse_fine_residual(n,mgt,brs_flx(n),temp_res(n),pdc)
     end do

     res_norm = ZERO
     do n = 1,nlevs
        res_norm = max(res_norm,norm_inf(res(n)))
     end do

     soln_norm = ZERO
     do n = 1,nlevs
        soln_norm = max(soln_norm,norm_inf(soln(n)))
     end do

     if ( mgt(nlevs)%verbose > 0 .and. parallel_IOProcessor() ) &
        write(unit=*, fmt='(i3,": Ninf(defect) = ",g15.8)') iter, res_norm

     all_done = res_norm < eps*(Anorm*soln_norm + bnorm)
     if (all_done) exit

      end do

!     Make sure even the coarse cells under fine cells have the best answer.
      do n = nlevs,2,-1
        mglev      = mgt(n)%nlevels
        mglev_crse = mgt(n-1)%nlevels
        hi_fine = upb(layout_get_pd(la_tower(n  ))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n-1))) + 1
        ref_ratio = hi_fine / hi_crse
        call ml_restriction(soln(n-1), soln(n), mgt(n)%mm(mglev), &
                            mgt(n-1)%mm(mglev_crse), mgt(n)%face_type, &
                            ref_ratio, inject = .true.)
      end do 

      do n = 1,nlevs
        call multifab_fill_boundary(soln(n))
      end do 

      if ( mgt(nlevs)%verbose > 0 .AND. parallel_IOProcessor() ) &
         write(unit=*, fmt='("MG finished at ", i3, " iterations")') iter

      if ( parallel_IOProcessor() ) then
         snrm(1) = ml_norm_l2(soln, fine_mask)
         snrm(2) = ml_norm_inf(soln, fine_mask)
         print *, 'SOLUTION MAX NORM ', snrm(2)
         print *, 'SOLUTION L2 NORM ', snrm(1)
      end if

      do n = 1,nlevs
         call multifab_destroy(uu(n))
         call multifab_destroy(res(n))
         call multifab_destroy(temp_res(n))
         call multifab_destroy(uu_hold(n))
         call lmultifab_destroy(fine_mask(n))
         if ( n > 1 ) then
            call bndry_reg_destroy(brs_flx(n))
         end if
      end do

      deallocate(uu, uu_hold, res, temp_res)
      deallocate(fine_mask, brs_flx)

   contains

        subroutine crse_fine_residual(n,mgt,brs_flx,temp_res,pdc)

           integer        , intent(in   ) :: n
           type(mg_tower) , intent(inout) :: mgt(:)
           type(bndry_reg), intent(inout) :: brs_flx
           type(multifab) , intent(inout) :: temp_res
           type(box)      , intent(in   ) :: pdc

           integer :: i,dm,mglev

           mglev = mgt(n)%nlevels
           dm = temp_res%dim

!          Zero out the flux registers which will hold the fine contributions
           call bndry_reg_setval(brs_flx, ZERO, all = .true.)

!          Compute the fine contributions at faces, edges and corners.

!          First compute a residual which only takes contributions from the
!             grid on which it is calculated.
           call grid_res(mgt(n),mglev,mgt(n)%ss(mglev),temp_res, &
                         rh(n),soln(n),mgt(n)%mm(mglev),mgt(n)%face_type)

           do i = 1,dm
              call ml_fine_contrib(brs_flx%bmf(i,0), &
                                   temp_res,mgt(n)%mm(mglev),ref_ratio,pdc,-i)
              call ml_fine_contrib(brs_flx%bmf(i,1), &
                                   temp_res,mgt(n)%mm(mglev),ref_ratio,pdc,+i)
           end do

!          Compute the crse contributions at edges and corners and add to res(n-1).
           do i = 1,dm
              call ml_crse_contrib(res(n-1), brs_flx%bmf(i,0), soln(n-1), &
                   mgt(n-1)%ss(mgt(n-1)%nlevels), &
                   mgt(n)%mm(mglev), &
                   pdc,ref_ratio, -i)
              call ml_crse_contrib(res(n-1), brs_flx%bmf(i,1), soln(n-1), &
                   mgt(n-1)%ss(mgt(n-1)%nlevels), &
                   mgt(n)%mm(mglev), &
                   pdc,ref_ratio, +i)
           end do
        end subroutine crse_fine_residual

        subroutine create_nodal_mask(n,mask,mm_crse,mm_fine,la_tower)

        integer        , intent(in   ) :: n
        type(lmultifab), intent(inout) :: mask
        type(imultifab), intent(inout) :: mm_crse
        type(imultifab), intent(inout) :: mm_fine
        type(layout   ), intent(inout) :: la_tower(:)
  
        integer        :: ref_ratio(mask%dim)
        integer        :: hi_fine(mask%dim)
        integer        :: hi_crse(mask%dim)
        type(box)      :: cbox,fbox

        logical, pointer :: mkp(:,:,:,:)
        integer, pointer :: cmp(:,:,:,:)
        integer, pointer :: fmp(:,:,:,:)

        integer :: loc(mask%dim),lof(mask%dim)
        integer :: i,j

        call setval(mask,.true.)

!       Note :          mm_fine is  in fine space
!       Note : mask and mm_crse are in crse space
  
        hi_fine = upb(layout_get_pd(la_tower(n+1))) + 1
        hi_crse = upb(layout_get_pd(la_tower(n  ))) + 1
        ref_ratio = hi_fine / hi_crse

        do j = 1,mask%nboxes

           cbox = get_ibox(mask,j)
           loc = lwb(cbox)

           do i = 1,mm_fine%nboxes

              fbox = get_ibox(mm_fine,i)
              lof = lwb(fbox)
              fbox = box_coarsen_v(fbox,ref_ratio)

              if (box_intersects(fbox,cbox)) then
                lo(:) = lwb(box_intersection(cbox,fbox))
                hi(:) = upb(box_intersection(cbox,fbox))

                mkp => dataptr(mask,j)
                cmp => dataptr(mm_crse,j)
                fmp => dataptr(mm_fine,i)
                select case (dm)
                case (2)
                   call create_nodal_mask_2d(mkp(:,:,1,1),cmp(:,:,1,1),loc,fmp(:,:,1,1),lof,lo,hi,ref_ratio)
                case (3)
                   call create_nodal_mask_3d(mkp(:,:,:,1),cmp(:,:,:,1),loc,fmp(:,:,:,1),lof,lo,hi,ref_ratio)
                end select
              end if
           end do
        end do

        end subroutine create_nodal_mask

        subroutine create_nodal_mask_2d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ref_ratio)

             integer, intent(in   ) :: loc(:),lof(:)
             logical, intent(inout) ::    mask(loc(1):,loc(2):)
             integer, intent(inout) :: mm_crse(loc(1):,loc(2):)
             integer, intent(inout) :: mm_fine(lof(1):,lof(2):)
             integer, intent(in   ) :: lo(:),hi(:)
             integer, intent(in   ) :: ref_ratio(:)

             integer :: i,j

             do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                if (.not.  bc_dirichlet(mm_fine(i*ref_ratio(1),j*ref_ratio(2)),1,0) .or. &
                           bc_dirichlet(mm_crse(i             ,j             ),1,0 ) ) &
                    mask(i,j) = .false.
             end do
             end do


        end subroutine create_nodal_mask_2d

        subroutine create_nodal_mask_3d(mask,mm_crse,loc,mm_fine,lof,lo,hi,ref_ratio)

             integer, intent(in   ) :: loc(:),lof(:)
             logical, intent(inout) ::    mask(loc(1):,loc(2):,loc(3):)
             integer, intent(inout) :: mm_crse(loc(1):,loc(2):,loc(3):)
             integer, intent(inout) :: mm_fine(lof(1):,lof(2):,lof(3):)
             integer, intent(in   ) :: lo(:),hi(:)
             integer, intent(in   ) :: ref_ratio(:)

             integer :: i,j,k,i_fine,j_fine,k_fine

             do k = lo(3),hi(3)
             do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                i_fine = i*ref_ratio(1)
                j_fine = j*ref_ratio(2)
                k_fine = k*ref_ratio(3)
                if (.not. bc_dirichlet(mm_fine(i_fine,j_fine,k_fine),1,0) .or. &
                          bc_dirichlet(mm_crse(i     ,j     ,k     ),1,0) ) &
                    mask(i,j,k) = .false.
             end do
             end do
             end do

        end subroutine create_nodal_mask_3d

   end subroutine ml_nd_solve

  function ml_converged(res, sol, mask, bnorm, Anorm, eps) result(r)
     logical :: r
     type(multifab), intent(in) :: res(:), sol(:)
     type(lmultifab), intent(in) :: mask(:)
     real(dp_t), intent(in) :: Anorm, eps, bnorm
     real(dp_t) :: ni_res, ni_sol
     ni_res = ml_norm_inf(res, mask)
     ni_sol = ml_norm_inf(sol, mask)
     r =  ni_res <= eps*(Anorm*ni_sol + bnorm) .or. &
          ni_res <= spacing(Anorm)
  end function ml_converged

  function ml_norm_inf(rr, mask) result(r)
     real(dp_t)  :: r
     type(multifab), intent(in) :: rr(:)
     type(lmultifab), intent(in) :: mask(:)
     integer :: n
     r = 0
     do n = 1, size(rr)
        r = max(norm_inf(rr(n),mask(n)), r)
        end do
  end function ml_norm_inf

  function ml_norm_l2(rr, mask) result(r)
     real(dp_t)  :: r
     type(multifab), intent(in) :: rr(:)
     type(lmultifab), intent(in) :: mask(:)
     integer :: n
     r = 0
     do n = 1, size(rr)
        r = max(norm_l2(rr(n),mask(n)), r)
     end do
  end function ml_norm_l2

end module ml_solve_module
