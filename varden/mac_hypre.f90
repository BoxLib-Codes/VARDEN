module mac_hypre_module
 
  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use fabio_module
 
  implicit none
 
  private
 
  public :: mac_hypre
 
contains

  subroutine mac_hypre(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                       stencil_order,ref_ratio,umac_norm)

    use stencil_fill_module, only: stencil_fill_cc_all_mglevels
    use mg_module, only: mg_tower, mg_tower_build, mg_tower_destroy
    use probin_module, only : verbose, pmask

    include 'HYPREf.h'

    type(ml_layout), intent(in   )        :: mla
    integer        , intent(in   )        :: stencil_order
    integer        , intent(in   )        :: ref_ratio(:,:)
    real(dp_t)     , intent(in)           :: dx(:,:)
    type(bc_tower) , intent(in)           :: the_bc_tower
    integer        , intent(in   )        :: bc_comp
    type(multifab) , intent(in   )        :: alpha(:), beta(:,:)
    type(multifab) , intent(inout)        ::    rh(:),  phi(:)
    type(bndry_reg), intent(inout)        :: fine_flx(2:)
    real(dp_t)     , intent(in), optional :: umac_norm(:)

    integer         :: lo(mla%dim),hi(mla%dim)
    integer         :: pdlo(mla%dim),pdhi(mla%dim)

    integer         :: d, i, n, ioff
    integer         :: dm, nlevs, ns_mg, ns_hy

    type(mg_tower)  :: mgt(mla%nlevel)

    type(layout)    :: la

    type(multifab), allocatable :: cell_coeffs(:)
    type(multifab), allocatable :: edge_coeffs(:,:)

    real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
    real(dp_t) :: pxa(mla%dim), pxb(mla%dim)
    real(dp_t) :: rel_eps, abs_eps

    integer(kind=8) :: grid
    integer(kind=8) :: A,b,x
    integer(kind=8) :: hypre_stencil
    integer(kind=8) :: solver,precond

    integer         :: ierr
    integer, allocatable :: stencil_indices(:)
    integer, allocatable :: offsets(:,:)
    integer, allocatable :: periodic_flag(:)
    type(box)       :: pd
    double precision :: tol
    double precision, allocatable :: values2d(:,:,:)
    double precision, allocatable :: values3d(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    real(kind=dp_t), pointer :: rhp(:,:,:,:), pp(:,:,:,:)
    real(kind=dp_t), pointer :: ssp(:,:,:,:)

    dm    = mla%dim
    nlevs = mla%nlevel

    if (nlevs > 1) &
       call bl_error('mac_hypre: not set up for nlevs > 1')

    if (nlevs .eq. 1) then
       rel_eps = 1.d-12
    else if (nlevs .eq. 2) then
       rel_eps = 1.d-11
    else
       rel_eps = 1.d-10
    endif
 
    abs_eps = -1.0_dp_t
    if (present(umac_norm)) then
       do n = 1,nlevs
          abs_eps = max(abs_eps, umac_norm(n) / dx(n,1))
       end do
       abs_eps = rel_eps * abs_eps
    end if

    ! We use standard cross stencils (with an extra place in each direction for boundary conditions in the mg stencil).
    ns_mg = 3*dm + 1
    ns_hy = 2*dm + 1

!   ******************************************************************************************************* 
!   Set up a mg_tower so that we can use that functionality to make the stencil
!   ******************************************************************************************************* 

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           dh = dx(n,:), &
                           ns = ns_mg, &
                           max_nlevel = 1, &
                           nodal = rh(nlevs)%nodal)

    end do

    do n = nlevs,1,-1

       allocate(cell_coeffs(mgt(n)%nlevels))
       allocate(edge_coeffs(mgt(n)%nlevels,dm))

       la = mla%la(n)

       call multifab_build(cell_coeffs(mgt(n)%nlevels), la, 1, 1)
       call multifab_copy_c(cell_coeffs(mgt(n)%nlevels),1,alpha(n),1, 1,ng=alpha(n)%ng)

       do d = 1, dm
          call multifab_build_edge(edge_coeffs(mgt(n)%nlevels,d),la,1,1,d)
          call multifab_copy_c(edge_coeffs(mgt(n)%nlevels,d),1,beta(n,d),1,1,ng=beta(n,d)%ng)
       end do

       if (n > 1) then
          xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
          xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
       else
          xa = ZERO
          xb = ZERO
       end if

       pxa = ZERO
       pxb = ZERO

       call stencil_fill_cc_all_mglevels(mgt(n), cell_coeffs, edge_coeffs, xa, xb, pxa, pxb, stencil_order, &
                                         the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))

       call destroy(cell_coeffs(mgt(n)%nlevels))
       deallocate(cell_coeffs)

       do d = 1, dm
          call destroy(edge_coeffs(mgt(n)%nlevels,d))
       end do
       deallocate(edge_coeffs)

    end do

!   ******************************************************************************************************* 
!   End of mg_tower stuff
!   ******************************************************************************************************* 

!   ******************************************************************************************************* 
!   Define a "grid" of dimension "dm" 
!   ******************************************************************************************************* 

    call HYPRE_StructGridCreate(MPI_COMM_WORLD,dm,grid,ierr) 

!   ******************************************************************************************************* 
!   Define the boxes with lo:hi
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, rh(n)%nboxes
         if ( multifab_remote(rh(n), i) ) cycle
         lo =  lwb(get_box(rh(n), i))
         hi =  upb(get_box(rh(n), i))
         call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
      end do
    end do

!   ******************************************************************************************************* 
!   Set the periodic flags correctly
!   ******************************************************************************************************* 

    pd = layout_get_pd(mla%la(1))
    pdlo =  lwb(pd)
    pdhi =  upb(pd)

    allocate(periodic_flag(3))
    periodic_flag(:) = 0
    do i = 1,dm
      if (pmask(i)) then
         periodic_flag(i) = pdhi(i)-pdlo(i)+1
!     else
!        periodic_flag(i) = 0
      end if
    end do

    call HYPRE_StructGridSetPeriodic(grid,periodic_flag)
    deallocate(periodic_flag)

!   ******************************************************************************************************* 
!   "Assemble" the grid.
!   ******************************************************************************************************* 

    call HYPRE_StructGridAssemble(grid,ierr)

!   ******************************************************************************************************* 
!   Define a stencil "hypre_stencil" with ns_hy entries in dm dimensions
!   ******************************************************************************************************* 

    call HYPRE_StructStencilCreate(dm,ns_hy,hypre_stencil,ierr)

!   ******************************************************************************************************* 
!   Define the offsets (locations) for the stencil entries
!   ******************************************************************************************************* 

    allocate(stencil_indices(ns_hy))
    allocate(offsets(dm,ns_hy))

    do i = 1, ns_hy
       stencil_indices(i) = i-1
    enddo
   
    ! Initialize to zero
    offsets(:,:) = 0
 
    ! Center
    ! All offsets are zero

    ! Left (lo i)
    offsets(1,2) = -1

    ! Right (hi i)
    offsets(1,3) =  1

    ! Down (lo j)
    offsets(2,4) = -1

    ! Up (hi j)
    offsets(2,5) =  1

    if (dm .eq. 3) then

       ! Back (lo k)
       offsets(3,6) = -1

       ! Front (hi k)
       offsets(3,7) = 1

    end if

    do i = 1, ns_hy
       call HYPRE_StructStencilSetElement(hypre_stencil,i-1,offsets(1:dm,i),ierr);
    end do
    deallocate(offsets)

!   ******************************************************************************************************* 
!   Link the matrix "A" to the stencil "hypre_stencil" and initialize
!   ******************************************************************************************************* 

    call HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,hypre_stencil,A,ierr)

    call HYPRE_StructMatrixInitialize(A,ierr)

!   ******************************************************************************************************* 
!   Define the elements of "A"
!   ******************************************************************************************************* 

    ! We have to use this offset to access into the mgt%ss array because it starts at 0 in stencil_fill codes...
    ioff = 1

    do n = 1,nlevs

      pd = layout_get_pd(mla%la(n))
      pdlo =  lwb(pd)
      pdhi =  upb(pd)

      do i = 1, mgt(n)%ss(1)%nboxes
         if ( multifab_remote(mgt(n)%ss(1), i) ) cycle
         lo =  lwb(get_box(mgt(n)%ss(1), i))
         hi =  upb(get_box(mgt(n)%ss(1), i))

         ssp => dataptr(mgt(n)%ss(1), i)

         if (dm.eq.2) then

            allocate(values2d(1:ns_hy,lo(1):hi(1),lo(2):hi(2)))
            values2d(:,:,:) = 0.d0

            ! Center
            values2d(1,lo(1):hi(1),lo(2):hi(2)) = ssp(lo(1):hi(1),lo(2):hi(2),1,ioff+0)

            ! Left (lo i)
            values2d(2,lo(1):hi(1),lo(2):hi(2)) = ssp(lo(1):hi(1),lo(2):hi(2),1,ioff+2)

            ! Right (hi i)
            values2d(3,lo(1):hi(1),lo(2):hi(2)) = ssp(lo(1):hi(1),lo(2):hi(2),1,ioff+1)

            ! Down (lo j)
            values2d(4,lo(1):hi(1),lo(2):hi(2)) = ssp(lo(1):hi(1),lo(2):hi(2),1,ioff+4)

            ! Up (hi j)
            values2d(5,lo(1):hi(1),lo(2):hi(2)) = ssp(lo(1):hi(1),lo(2):hi(2),1,ioff+3)

            call HYPRE_StructMatrixSetBoxValues(A,lo,hi,ns_hy,stencil_indices,values2d,ierr)
            deallocate(values2d)

         else if (dm.eq.3) then

            allocate(values3d(1:ns_hy,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
            values3d(:,:,:,:) = 0.d0

            ! Center
            values3d(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+0)

            ! Left (lo i)
            values3d(2,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+2)

            ! Right (hi i)
            values3d(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+1)

            ! Down (lo j)
            values3d(4,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+4)

            ! Up (hi j)
            values3d(5,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+3)

            ! Back (lo k)
            values3d(6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+6)

            ! Front (hi k)
            values3d(7,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ioff+5)

            call HYPRE_StructMatrixSetBoxValues(A,lo,hi,ns_hy,stencil_indices,values3d,ierr)
            deallocate(values3d)

         end if
       end do
    end do

    deallocate(stencil_indices)

!   ******************************************************************************************************* 
!   "Assemble" the matrix A
!   ******************************************************************************************************* 

    call HYPRE_StructMatrixAssemble(A,ierr)

!   ******************************************************************************************************* 
!   Create vectors "x" and "b".
!   ******************************************************************************************************* 

    ! Create an empty vector object
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, b, ierr)
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, x, ierr)
 
    ! Indicate that the vector coefficients are ready to be set
    call HYPRE_StructVectorInitialize(b,ierr)
    call HYPRE_StructVectorInitialize(x,ierr)

!   ******************************************************************************************************* 
!   Fill vector "b" with the "rh" multifab.
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, rh(n)%nboxes
         if ( multifab_remote(rh(n), i) ) cycle
         lo =  lwb(get_box(rh(n), i))
         hi =  upb(get_box(rh(n), i))

         rhp => dataptr(rh(n),i)

         ! Set RHS
         if (dm.eq.2) then
            allocate(values2d(1,lo(1):hi(1),lo(2):hi(2)))
            values2d(1,lo(1):hi(1),lo(2):hi(2)) = rhp(lo(1):hi(1),lo(2):hi(2),1,1)

            call HYPRE_StructVectorSetBoxValues(b, lo, hi, values2d, ierr)

            deallocate(values2d)

         else if (dm.eq.3) then

            allocate(values3d(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
            values3d(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = rhp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

            call HYPRE_StructVectorSetBoxValues(b, lo, hi, values3d, ierr)

            deallocate(values3d)
         end if

      end do
    end do

    ! This is a collective call finalizing the vector assembly.
    !  The vectors are now ``ready to be used''
    call HYPRE_StructVectorAssemble(b,ierr)
    call HYPRE_StructVectorAssemble(x,ierr)

!   Print the matrix "A"
!   call HYPRE_StructMatrixPrint(A,0,ierr)

!   Print the right-hand-side "b"
!   call HYPRE_StructVectorPrint(b,0,ierr)

!   Create an empty PCG Struct solver
    call HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ierr)

!   Set PCG parameters
    tol = rel_eps
    call HYPRE_StructPCGSetTol(solver, tol, ierr)
    call HYPRE_StructPCGSetPrintLevel(solver, 2, ierr)
    call HYPRE_StructPCGSetMaxIter(solver, 50, ierr)
 
!   Create the Struct PFMG solver for use as a preconditioner
    call HYPRE_StructPFMGCreate(MPI_COMM_WORLD, precond, ierr)
 
!   Set PFMG parameters
    call HYPRE_StructPFMGSetMaxIter(precond, 1, ierr)
    call HYPRE_StructPFMGSetTol(precond, 0.0, ierr)
    call HYPRE_StructPFMGSetZeroGuess(precond, ierr)
    call HYPRE_StructPFMGSetNumPreRelax(precond, 2, ierr)
    call HYPRE_StructPFMGSetNumPostRelax(precond, 2, ierr)
 
!   Non-Galerkin coarse grid (more efficient for this problem)
    call HYPRE_StructPFMGSetRAPType(precond, 1, ierr)
 
!   R/B Gauss-Seidel
    call HYPRE_StructPFMGSetRelaxType(precond, 2, ierr)
 
!   Skip relaxation on some levels (more efficient for this problem)
    call HYPRE_StructPFMGSetSkipRelax(precond, 1, ierr)

!   Set preconditioner (PFMG = 1) and solve
    call HYPRE_StructPCGSetPrecond(solver, 1, precond, ierr)
    call HYPRE_StructPCGSetup(solver, A, b, x, ierr)

    call HYPRE_StructPCGSolve(solver, A, b, x, ierr)

!   Print the solution "x"
!   call HYPRE_StructVectorPrint(x,0,ierr)

!   ******************************************************************************************************* 
!   Fill multifab "phi" from the vector solution "x".
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, phi(n)%nboxes
         if ( multifab_remote(phi(n), i) ) cycle
         lo =  lwb(get_box(phi(n), i))
         hi =  upb(get_box(phi(n), i))

         pp => dataptr(phi(n),i)

         ! Set RHS
         if (dm.eq.2) then
            allocate(values2d(1,lo(1):hi(1),lo(2):hi(2)))
            values2d(:,:,:) = 1.d20

            call HYPRE_StructVectorGetBoxValues(x, lo, hi, values2d, ierr)

            pp(lo(1):hi(1),lo(2):hi(2),1,1) = values2d(1,lo(1):hi(1),lo(2):hi(2)) 

            deallocate(values2d)

         else if (dm.eq.3) then

            allocate(values3d(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
            values3d(:,:,:,:) = 1.d20

            call HYPRE_StructVectorGetBoxValues(x, lo, hi, values3d, ierr)

            pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = values3d(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) 

            deallocate(values3d)
         end if

      end do
    end do

    call fabio_multifab_write_d(phi(1),'HYPRE_PHI','Phi')

!   Free memory
!   call HYPRE_StructPCGDestroy(solver);
    call HYPRE_StructPFMGDestroy(precond)

    call HYPRE_StructGridDestroy(grid);
    call HYPRE_StructStencilDestroy(hypre_stencil);

    call HYPRE_StructMatrixDestroy(A);
    call HYPRE_StructVectorDestroy(b);
    call HYPRE_StructVectorDestroy(x);

    call build(bpt, "mac_hypre")

    call destroy(bpt)
    stop

  end subroutine mac_hypre
end module mac_hypre_module