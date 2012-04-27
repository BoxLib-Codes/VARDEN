module hg_hypre_module
 
  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
 
  implicit none
 
  private
 
  public :: hg_hypre

contains

  subroutine hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                      press_comp,stencil_type, &
                      rel_solver_eps, abs_solver_eps, divu_rhs)

    use enforce_outflow_on_divu_module, only : enforce_outflow_on_divu_rhs

    use stencil_fill_module , only : stencil_fill_nodal_all_mglevels, stencil_fill_one_sided
    use nodal_divu_module   , only : divu, subtract_divu_from_rh
    use probin_module       , only : pmask
    use mg_module           

    include 'HYPREf.h'

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: press_comp
    integer        , intent(in   ) :: stencil_type
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    type(multifab ), intent(inout), optional :: divu_rhs(:)

    type(box     )  :: pd
    type(  layout)  :: la

    integer         :: lo(3),hi(3)
    integer         :: pdlo(mla%dim),pdhi(mla%dim)

    type(mg_tower) :: mgt(mla%nlevel)

    type(multifab) :: one_sided_ss(2:mla%nlevel)
    type(multifab), allocatable :: coeffs(:)

    logical :: nodal(mla%dim)
    integer :: i,n,dm,nlevs,ns_mg,ns_hy

    ! All the integers associated with Hypre are long.
    integer(kind=8) :: grid, graph
    integer(kind=8) :: A,b,x     ! These are SStruct versions
    integer(kind=8) :: sA,sb,sx  ! These are  Struct versions
    integer(kind=8) :: hypre_stencil
    integer(kind=8) :: solver,precond
 
    integer              :: vartypes(1)
    integer              :: var,nvars,part,nparts,object_type,ierr
    integer              :: periodic_flag(3)
    integer, allocatable :: stencil_indices(:)
    integer, allocatable :: offsets(:,:)
 
    double precision :: tol
    double precision, allocatable :: values(:,:,:,:)

!   This comes from 'sstruct_mv/HYPRE_sstruct_mv.h'
    integer, parameter ::  HYPRE_SSTRUCT_VARIABLE_NODE = 1 

    real(kind=dp_t), pointer :: rhp(:,:,:,:), pp(:,:,:,:)
    real(kind=dp_t), pointer :: ssp(:,:,:,:)

    !! Defaults:

    dm    = mla%dim
    nlevs = mla%nlevel

    call print(rhohalf(1),'RHOHALF')

    ! We dimension these for use in 2D
    lo(3) = 1
    hi(3) = 1

    var    = 0 ! We only have one variable-type (nodal)
    nvars  = 1 ! We only have one variable-type (nodal)

    part   = 0 ! We only have one region (though it's divided into parts)
    nparts = 1 ! We only have one region (though it's divided into parts)

    object_type = HYPRE_STRUCT

    nodal = .true.

    if (stencil_type .eq. ST_DENSE) then
       if (dm .eq. 3) then
          i = mgt(nlevs)%nlevels
          if ( (dx(nlevs,1) .eq. dx(nlevs,2)) .and. &
               (dx(nlevs,1) .eq. dx(nlevs,3)) ) then
             ns_mg = 21
             ns_hy = 21
          else
             ns_mg = 27
             ns_hy = 27
          end if
       else if (dm .eq. 2) then
          ns_mg = 9
          ns_hy = 9
       end if

    else if (stencil_type .eq. ST_CROSS) then

       ns_mg = 2*dm+1
       ns_hy = 2*dm+1
       do n = nlevs, 2, -1
          la = mla%la(n)
          call multifab_build(one_sided_ss(n), la, ns_mg, 0, nodal, stencil=.true.)
       end do

    else
       call bl_error("HG_HYPRE: dont know this stencil_type")
    end if

!   *******************************************************************************************************
!   Set up a mg_tower so that we can use that functionality to make the stencil
!   *******************************************************************************************************

    do n = nlevs, 1, -1

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,press_comp), &
                           dh = dx(n,:), &
                           ns = ns_mg, &
                           max_nlevel = 1, &
                           nodal = nodal)
       
    end do

    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1, 1)
       call setval(coeffs(mgt(n)%nlevels), 0.0_dp_t, 1, all=.true.)

       call mkcoeffs(rhohalf(n),coeffs(mgt(n)%nlevels))
       call multifab_fill_boundary(coeffs(mgt(n)%nlevels))

       call stencil_fill_nodal_all_mglevels(mgt(n), coeffs, stencil_type)

       if (stencil_type .eq. ST_CROSS .and. n .gt. 1) then
          i = mgt(n)%nlevels
          call stencil_fill_one_sided(one_sided_ss(n), coeffs(i), &
                                      mgt(n    )%dh(:,i), &
                                      mgt(n)%mm(i), mgt(n)%face_type)
       end if

       call destroy(coeffs(mgt(n)%nlevels))
       deallocate(coeffs)

    end do

!   *******************************************************************************************************
!   Create the rhs
!   *******************************************************************************************************

    call divu(nlevs,mgt,unew,rh,mla%mba%rr,nodal)
 
    ! Do rh = rh - divu_rhs (this routine preserves rh=0 on
    !  nodes which have bc_dirichlet = true.
    if (present(divu_rhs)) then
       call enforce_outflow_on_divu_rhs(divu_rhs,the_bc_tower)
       call subtract_divu_from_rh(nlevs,mgt,rh,divu_rhs)
    end if
 
!   *******************************************************************************************************
!   End of mg_tower stuff
!   *******************************************************************************************************

    do n = nlevs,1,-1
       call multifab_fill_boundary(phi(n))
    end do

!   ******************************************************************************************************* 
!   Define a "grid" of dimension "dm" based on nparts parts
!   ******************************************************************************************************* 

    call HYPRE_SStructGridCreate(MPI_COMM_WORLD,dm,nparts,grid,ierr) 

!   ******************************************************************************************************* 
!   Define the boxes with lo:hi
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, nboxes(rh(n))
         if ( multifab_remote(rh(n), i) ) cycle
         lo(1:dm) =  lwb(get_box(rh(n), i))
         hi(1:dm) =  upb(get_box(rh(n), i))
         print *,'CC LO ',lo(:)
         print *,'CC HI ',hi(:)
         call HYPRE_SStructGridSetExtents(grid,part,lo,hi,ierr)
      end do
    end do

!   Set the variable type to be nodal on this (only) part
    vartypes(1) = HYPRE_SSTRUCT_VARIABLE_NODE
    call HYPRE_SStructGridSetVariables(grid,part,nvars,vartypes,ierr)

!   ******************************************************************************************************* 
!   Set the periodic flags correctly
!   ******************************************************************************************************* 

    pd = layout_get_pd(mla%la(1))
    pdlo =  lwb(pd)
    pdhi =  upb(pd)

    periodic_flag(:) = 0
    do i = 1,dm
      if (pmask(i)) then
         periodic_flag(i) = pdhi(i)-pdlo(i)+1
      end if
    end do

    call HYPRE_SStructGridSetPeriodic(grid,part,periodic_flag,ierr)

!   ******************************************************************************************************* 
!   "Assemble" the grid.
!   ******************************************************************************************************* 

    call HYPRE_SStructGridAssemble(grid,ierr)

!   ******************************************************************************************************* 
!   Define a stencil "hypre_stencil" with ns_hy entries in dm dimensions
!   ******************************************************************************************************* 

    call HYPRE_SStructStencilCreate(dm,ns_hy,hypre_stencil,ierr)

!   ******************************************************************************************************* 
!   Define the offsets (locations) for the stencil entries
!   ******************************************************************************************************* 

    allocate(stencil_indices(ns_hy))
    allocate(offsets(dm,ns_hy))

    do i = 1, ns_hy
       stencil_indices(i) = i-1
    enddo
   
    ! Initialize to zero
    ! First index is i/j/k, second index is ordering of point in stencil
    offsets(:,:) = 0
 
    ! Center
    offsets(:,1) = 0

    if (dm.eq.2) then

       if (ns_hy .eq. 5) then

          ! Left (lo i, med j)
          offsets(1,2) = -1

          ! Right (hi i, med j)
          offsets(1,3) =  1

          ! Down (med i, lo j)
          offsets(2,4) = -1

          ! Up (med i, hi j)
          offsets(2,5) =  1

       else if (ns_hy .eq. 9) then

          ! (lo i, lo j)
          offsets(1,2) = -1
          offsets(2,2) = -1

          ! (med i, lo j)
          offsets(1,3) =  0
          offsets(2,3) = -1

          ! (hi i, lo j)
          offsets(1,4) =  1
          offsets(2,4) = -1

          ! (lo i, med j)
          offsets(1,5) = -1
          offsets(2,5) =  0

          ! (hi i, med j)
          offsets(1,6) =  1
          offsets(2,6) =  0

          ! (lo i, hi j)
          offsets(1,7) = -1
          offsets(2,7) =  1

          ! (med i, hi j)
          offsets(1,8) =  0
          offsets(2,8) =  1

          ! (hi i, hi j)
          offsets(1,9) =  1
          offsets(2,9) =  1

       else 

          call bl_error("HG_HYPRE: havent implemented stencil for this 2d ns_hy")

       end if
 
    else if (dm .eq. 3) then

       if (ns_hy .eq. 7) then

          ! Left (lo i, med j, med k)
          offsets(1,2) = -1

          ! Right (hi i, med j, med k)
          offsets(1,3) =  1

          ! Down (med i, lo j, med k)
          offsets(2,4) = -1

          ! Up (med i, hi j, med k)
          offsets(2,5) =  1

          ! Back (med i, med j, lo k)
          offsets(3,6) = -1

          ! Front (med i, med j, hi k)
          offsets(3,7) =  1

       else 
          call bl_error("HG_HYPRE: havent implemented stencil for this 3d ns_hy")
       end if

    end if

    do i = 1, ns_hy
       call HYPRE_SStructStencilSetEntry(hypre_stencil,i-1,offsets(1:dm,i),var,ierr);
    end do
    deallocate(offsets)

!   ******************************************************************************************************* 
!    Set up the Graph - this determines the non-zero structure of
!        the matrix and allows non-stencil relationships between the parts
!   ******************************************************************************************************* 
 
!     Create the graph object
      call HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, graph, ierr)
 
!     Now we need to tell the graph which stencil to use for each
!     variable on each part (we only have one variable and one part)
      call HYPRE_SStructGraphSetStencil(graph, part, var, hypre_stencil, ierr)
 
!     Assemble the graph
      call HYPRE_SStructGraphAssemble(graph, ierr)

!   ******************************************************************************************************* 
!   Create an empty matrix object and initialize
!   ******************************************************************************************************* 

    call HYPRE_SStructMatrixCreate(MPI_COMM_WORLD,graph,A,ierr)
    call HYPRE_SStructMatrixSetObjectTyp(A,object_type,ierr)
    call HYPRE_SStructMatrixInitialize(A,ierr)

!   ******************************************************************************************************* 
!   Define the elements of "A"
!   ******************************************************************************************************* 

    do n = 1,nlevs

      pd = layout_get_pd(mla%la(n))
      pdlo =  lwb(pd)
      pdhi =  upb(pd)

      do i = 1, nboxes(mgt(n)%ss(1))
         if ( multifab_remote(mgt(n)%ss(1), i) ) cycle
         lo(1:dm) =  lwb(get_ibox(mgt(n)%ss(1), i)) - 1  ! Subtract one because of the indexing convection
         hi(1:dm) =  upb(get_ibox(mgt(n)%ss(1), i)) - 1  ! Subtract one because of the indexing convection

         print *,'ND LO ',lo(:)
         print *,'ND HI ',hi(:)

         ssp => dataptr(mgt(n)%ss(1), i)

         allocate(values(1:ns_hy,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         values(:,:,:,:) = 0.d0

         if (dm.eq.2) then

           if (ns_hy .eq. 5) then

            ! Center
            values(1,lo(1):hi(1),lo(2):hi(2),1) = ssp(lo(1):hi(1),lo(2):hi(2),1,1)

            ! Left (lo i)
            values(2,lo(1):hi(1),lo(2):hi(2),1) = ssp(lo(1):hi(1),lo(2):hi(2),1,3)

            ! Right (hi i)
            values(3,lo(1):hi(1),lo(2):hi(2),1) = ssp(lo(1):hi(1),lo(2):hi(2),1,2)

            ! Down (lo j)
            values(4,lo(1):hi(1),lo(2):hi(2),1) = ssp(lo(1):hi(1),lo(2):hi(2),1,5)

            ! Up (hi j)
            values(5,lo(1):hi(1),lo(2):hi(2),1) = ssp(lo(1):hi(1),lo(2):hi(2),1,4)

           else
              call bl_error("HG_HYPRE: dont know this ns_hy when setting stencil")
           end if

         else if (dm.eq.3) then

           if (ns_hy .eq. 7) then

            ! Center
            values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

            ! Left (lo i)
            values(2,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3)

            ! Right (hi i)
            values(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),2)

            ! Down (lo j)
            values(4,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

            ! Up (hi j)
            values(5,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),4)

            ! Back (lo k)
            values(6,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),7)

            ! Front (hi k)
            values(7,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = ssp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),6)

           else
              call bl_error("HG_HYPRE: dont know this ns_hy when setting stencil")
           end if

         end if

         call HYPRE_SStructMatrixSetBoxValues(A,part,lo,hi,var,ns_hy,stencil_indices,values,ierr)
         deallocate(values)

       end do
    end do

    deallocate(stencil_indices)

!   ******************************************************************************************************* 
!   "Assemble" the matrix A
!   ******************************************************************************************************* 

    call HYPRE_SStructMatrixSetObjectTyp(A, object_type, ierr)

    call HYPRE_SStructMatrixAssemble(A,ierr)

!   ******************************************************************************************************* 
!   Create vectors "x" and "b".
!   ******************************************************************************************************* 

    ! Create an empty vector object
    call HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, b, ierr)
    call HYPRE_SStructVectorCreate(MPI_COMM_WORLD, grid, x, ierr)

!   As with the matrix, set the appropriate object type for the vectors
    call HYPRE_SStructVectorSetObjectTyp(b, object_type, ierr)
    call HYPRE_SStructVectorSetObjectTyp(x, object_type, ierr)

    ! Indicate that the vector coefficients are ready to be set
    call HYPRE_SStructVectorInitialize(b,ierr)
    call HYPRE_SStructVectorInitialize(x,ierr)

!   ******************************************************************************************************* 
!   Fill vector "b" with the "rh" multifab.
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, nboxes(rh(n))
         if ( multifab_remote(rh(n), i) ) cycle
         lo(1:dm) =  lwb(get_ibox(rh(n), i))
         hi(1:dm) =  upb(get_ibox(rh(n), i))

         rhp => dataptr(rh(n),i)

         allocate(values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

         ! Set RHS
         if (dm.eq.2) then
            values(1,lo(1):hi(1),lo(2):hi(2),1) = &
                 rhp(lo(1):hi(1),lo(2):hi(2),1,1)

         else if (dm.eq.3) then

            values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = &
                 rhp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)

         end if

         call HYPRE_SStructVectorSetBoxValues(b, part,lo, hi, var, values, ierr)
         deallocate(values)

      end do
    end do

    ! This is a collective call finalizing the vector assembly.
    !  The vectors are now ``ready to be used''
    call HYPRE_SStructVectorAssemble(b,ierr)
    call HYPRE_SStructVectorAssemble(x,ierr)

!   ******************************************************************************************************* 
!   Convert to Struct form so we can ues PFMG
!   ******************************************************************************************************* 

!   Because we are using a struct solver, we need to get the object 
!      of the matrix and vectors to pass in to the struct solvers  
    call HYPRE_SStructMatrixGetObject(A, sA, ierr) 
    call HYPRE_SStructVectorGetObject(b, sb, ierr) 
    call HYPRE_SStructVectorGetObject(x, sx, ierr) 

!   Print the matrix "A"
    call HYPRE_StructMatrixPrint(sA,0,ierr)

!   Print the right-hand-side "b"
!   call HYPRE_StructVectorPrint(sb,0,ierr)

!   Create an empty PCG Struct solver
    call HYPRE_StructPCGCreate(MPI_COMM_WORLD, solver, ierr)

!   Set PCG parameters
    tol = rel_solver_eps
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
    call HYPRE_StructPCGSetup(solver, sA, sb, sx, ierr)

    call HYPRE_StructPCGSolve(solver, sA, sb, sx, ierr)

!   ******************************************************************************************************* 
!   Fill multifab "phi" from the vector solution "x".
!   ******************************************************************************************************* 

    do n = 1,nlevs
      do i = 1, nboxes(phi(n))
         if ( multifab_remote(phi(n), i) ) cycle
         lo(1:dm) =  lwb(get_ibox(phi(n), i))
         hi(1:dm) =  upb(get_ibox(phi(n), i))

         pp => dataptr(phi(n),i)

         ! Set RHS
         if (dm.eq.2) then
            allocate(values(1,lo(1):hi(1),lo(2):hi(2),1))
            values(:,:,:,:) = 1.d20

            call HYPRE_StructVectorGetBoxValues(sx, lo, hi, values, ierr)

            pp(lo(1):hi(1),lo(2):hi(2),1,1) = values(1,lo(1):hi(1),lo(2):hi(2),1) 

         else if (dm.eq.3) then

            allocate(values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
            values(:,:,:,:) = 1.d20

            call HYPRE_StructVectorGetBoxValues(sx, lo, hi, values, ierr)

            pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1) = values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) 

         end if

         deallocate(values)

      end do
    end do

!   Free memory
    call HYPRE_StructPCGDestroy(solver,ierr);
    call HYPRE_StructPFMGDestroy(precond,ierr)

    call HYPRE_SStructGridDestroy(grid,ierr);
    call HYPRE_SStructStencilDestroy(hypre_stencil,ierr);

!   call HYPRE_SStructMatrixDestroy(A,ierr);
!   call HYPRE_SStructVectorDestroy(b,ierr);
!   call HYPRE_SStructVectorDestroy(x,ierr);

    call HYPRE_StructMatrixDestroy(sA,ierr);
    call HYPRE_StructVectorDestroy(sb,ierr);
    call HYPRE_StructVectorDestroy(sx,ierr);

    if (stencil_type .ne. ST_DENSE) then
       do n = nlevs, 2, -1
          call multifab_destroy(one_sided_ss(n))
       end do
    endif

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

  end subroutine hg_hypre

  !   ********************************************************************************** !

  subroutine mkcoeffs(rho,coeffs)

    type(multifab) , intent(in   ) :: rho
    type(multifab) , intent(inout) :: coeffs

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)
    integer :: i,dm,ng

    dm = get_dim(rho)
    ng = nghost(rho)

    do i = 1, nboxes(rho)
       if ( multifab_remote(rho, i) ) cycle
       rp => dataptr(rho   , i)
       cp => dataptr(coeffs, i)
       select case (dm)
       case (2)
          call mkcoeffs_2d(cp(:,:,1,1), rp(:,:,1,1), ng)
       case (3)
          call mkcoeffs_3d(cp(:,:,:,1), rp(:,:,:,1), ng)
       end select
    end do

  end subroutine mkcoeffs

  !   *********************************************************************************** !

  subroutine mkcoeffs_2d(coeffs,rho,ng)

      use bl_constants_module

    integer :: ng
    real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:)
    real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:)

    integer :: i,j
    integer :: nx,ny

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2

    do j = 1,ny
       do i = 1,nx
          coeffs(i,j) = ONE / rho(i,j)
       end do
    end do

  end subroutine mkcoeffs_2d

  !   ********************************************************************************** !

  subroutine mkcoeffs_3d(coeffs,rho,ng)

      use bl_constants_module

    integer :: ng
    real(kind=dp_t), intent(inout) :: coeffs(   0:,   0:, 0:)
    real(kind=dp_t), intent(in   ) ::  rho(1-ng:,1-ng:,1-ng:)

    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(coeffs,dim=1) - 2
    ny = size(coeffs,dim=2) - 2
    nz = size(coeffs,dim=3) - 2

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             coeffs(i,j,k) = ONE / rho(i,j,k)
          end do
       end do
    end do

  end subroutine mkcoeffs_3d

end module hg_hypre_module
