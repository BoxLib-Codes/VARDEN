module scalar_advance_sdc_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use explicit_diffusive_module
  use viscous_module
  use macproject_module
  use rxns_integrator
  use sdc_interpolation
  use probin_module, only : nscal, diff_coef, diffusion_type, stencil_order, &
                            verbose, mg_verbose, cg_verbose, sdc_iters, &
                            mass_fractions, nspec

  implicit none

  public :: scalar_advance_sdc


contains


  subroutine scalar_advance_sdc(mla,uold,sold,snew,umac, &
                                ext_scal_force,dx,dt,t,the_bc_tower)

    use mkflux_module
    use mkforce_module
    use update_module
    use bl_constants_module
    use reaction_fn

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,t
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    type(multifab), allocatable :: scal_force(:), divu(:)
    type(multifab), allocatable :: sflux(:,:), sedge(:,:)
    type(multifab), allocatable :: source(:),I_ADR(:)
    logical       , allocatable :: is_conservative(:)
    logical       , allocatable :: umac_nodal_flag(:)
    ! advection term:  adv_s(time,level) = - div(su) 
    ! these values are used to interpolate int(A+D+R)
    type(multifab) , allocatable :: adv_s(:,:)
    ! diffusion at various times:  D_s(time,level)
    type(multifab) , allocatable :: D_s(:,:)

    integer         :: i,j,k,n
    integer         :: comp,bc_comp,nlevs,dm
    logical         :: is_vel
    real(kind=dp_t) :: diff_fac,visc_mu
    real(kind=dp_t) :: smin,smax

    ! for convergence study--DELETE ME
    type(multifab), allocatable :: snew_old(:), difference(:)


    nlevs = mla%nlevel
    dm    = mla%dim

    allocate(umac_nodal_flag(dm))
    allocate(scal_force(nlevs),divu(nlevs),source(nlevs))
    allocate(sflux(nlevs,dm),sedge(nlevs,dm))
    allocate(is_conservative(nscal))
    allocate(D_s(0:n_diff-1,nlevs),adv_s(0:n_adv-1,nlevs))
    allocate(I_ADR(nlevs))    
    allocate(snew_old(nlevs),difference(nlevs))


    if (mass_fractions) then
       ! using mass fractions
       is_conservative(:) = .true.
    else
       is_conservative(1) = .true.
       ! using concentrations
       do comp = 2, nscal
          is_conservative(comp) = .false.
       end do
    end if

    is_vel  = .false.

    do n = 1, nlevs
       call multifab_build( scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build( divu(n),mla%la(n),    1,1)
! have no ghost cells--may need to change
       call multifab_build( source(n),mla%la(n),nspec,1)
       do j = 0, n_diff-1
          call multifab_build( D_s(j,n),mla%la(n),nspec,0)
          call setval(D_s(j,n),zero,all=.true.)
       enddo
       do j = 0, n_adv-1
          call multifab_build( adv_s(j,n),mla%la(n),nspec,0)
          call setval(adv_s(j,n),zero,all=.true.)
       enddo
       call multifab_build( I_ADR(n),mla%la(n),nspec,0)
       call setval(I_ADR(n),zero,all=.true.)

       do i = 1,dm
         umac_nodal_flag(:) = .false.
         umac_nodal_flag(i) = .true.
         call multifab_build(sflux(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call multifab_build(sedge(n,i),scal_force(n)%la,nscal,0,nodal = umac_nodal_flag)
         call setval(sflux(n,i),ZERO,all=.true.)
         call setval(sedge(n,i),ZERO,all=.true.)
       end do

       call setval(scal_force(n),0.0_dp_t,all=.true.)
       call setval(divu(n),0.0_dp_t,all=.true.)
       call setval(source(n),0.0_dp_t,all=.true.)

       call multifab_build(snew_old(n),mla%la(n),  nscal,0)
       call multifab_build(difference(n),mla%la(n),nscal,0)
       call copy(snew_old(n), sold(n), 0)
    enddo
   


   !--Compute provisional solution at t=n+1 using forward/backward euler------------

    !*********************************************
    ! Compute D(s) = laplacian(s) at time n
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(0,:),sold,comp,bc_comp,dx,the_bc_tower)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(0,n),diff_coef)
       enddo
    endif

    !************************************************************
    ! Create scalar force at time n.
    ! Provided as a source term to calculate advection term adv_s
    !************************************************************

    ! only source is ext_scal_force
    do n = 1, nlevs
       call multifab_copy_c(scal_force(n),1,ext_scal_force(n),1,&
                            nscal,ext_scal_force(n)%ng)
    enddo

    ! another option to do the same thing
    !diff_fac = zero     !in mkscalforce: diff_fac*source
    !call mkscalforce(nlevs,scal_force,ext_scal_force,source,diff_fac) 
    ! NOTE: use mksource when source mf only holds the reactive speices

    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************

    call mkflux(mla,sold,uold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_tower%bc_tower_array,is_vel,is_conservative)

    !*****************************************************************
    ! Update the scalars with conservative or convective differencing.
    !*****************************************************************

    call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s(0,:),dx,dt,t,& 
                is_vel,is_conservative,the_bc_tower%bc_tower_array)

    if (verbose .ge. 1) then
       do n = 1, nlevs
          do comp = 1, nscal
             smin = multifab_min_c(snew(n),comp) 
             smax = multifab_max_c(snew(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2000) n,smin,smax
             else if (comp .eq. 2) then
                if (parallel_IOProcessor()) write(6,2001) n,smin,smax
             end if
          end do
       end do
    end if

    !***************************************
    ! Solve u_t = A + D using implicit euler
    !***************************************

    do n = 1,nlevs
       ! just a dummy to pass in zero to diffusion solve
       call setval(scal_force(n),zero)
    end do

    if (diff_coef > ZERO) then
       do  comp = 2, nscal
           bc_comp = dm+comp
           visc_mu = dt*diff_coef
           call diff_scalar_solve(mla,snew,scal_force,dx,t,visc_mu,&
                                  the_bc_tower,comp,bc_comp)
        end do
    end if

    !*********************************************
    ! Compute a provisional D(s) at time n+1
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(1,:),snew,comp,bc_comp,&
                                           dx,the_bc_tower)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(1,n),diff_coef)
       enddo
    endif

    !***********************
    ! Integrate reactions
    !***********************
    
    do n = 1, nlevs
       call multifab_copy_c(snew(n),1,sold(n),1,nscal,sold(n)%ng)
    enddo

    call react(mla,the_bc_tower,snew,dx,dt,t,provisional,adv_s,D_s)

    !*********************************************
    ! Compute D(s) at time n+1
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(2,:),snew,comp,bc_comp,dx,the_bc_tower)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(2,n),diff_coef)
       enddo
    endif
    
    !*****************************
    ! Compute aprrox to int(A+D+R)
    !*****************************
    call provisional_intgrl(I_ADR,sold,snew,adv_s,D_s,dt,nlevs)

   !----SDC iters----------------------------------------------------
    do k = 1,sdc_iters
       
      ! Update advection:

       !******************************************************
       ! Create forcing term at time n.
       ! used as source term to calculate advection term adv_s

       diff_fac = ONE
       call mksource(nlevs,scal_force,ext_scal_force,I_ADR,diff_fac)

       !***********************************
       ! Create edge state scalars/fluxes
       
       call mkflux(mla,sold,uold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                   the_bc_tower%bc_tower_array,is_vel,is_conservative)
       
       !*************************************************************
       ! Create scalar force at time n+1/2
       ! Using implicit euler, is equal to the external forcing terms
  
       ! forcing term doesn't need to be changed
       !do n = 1, nlevs
       !   call multifab_copy_c(scal_force(n),1,ext_scal_force(n),1,&
       !                        nscal,ext_scal_force(n)%ng)
       !enddo

       !*****************************************************************
       ! Update the scalars with conservative or convective differencing.

       call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s(0,:),dx,dt,t,&
                   is_vel,is_conservative,the_bc_tower%bc_tower_array)

       if (verbose .ge. 1) then
          do n = 1, nlevs
             do comp = 1, nscal
                smin = multifab_min_c(snew(n),comp) 
                smax = multifab_max_c(snew(n),comp)
                if (comp .eq. 1) then
                   if (parallel_IOProcessor()) write(6,2000) n,smin,smax
                else if (comp .eq. 2) then
                   if (parallel_IOProcessor()) write(6,2001) n,smin,smax
                end if
             end do
          end do
       end if
       
      
      ! Update Diffusion:
       !*****************
       ! Solve snew = sold +dt*(delta A(s) + delta D(s))

       do n = 1, nlevs
          call multifab_copy_c(source(n),1,D_s(1,n),1,nscal-1)
          call multifab_mult_mult_s(source(n),-1.d0)
       enddo

       if (diff_coef > ZERO) then
          do  comp = 2, nscal
             bc_comp = dm+comp
             visc_mu = dt*diff_coef
             call diff_scalar_solve(mla,snew,source,dx,t,visc_mu,&
                                    the_bc_tower,comp,bc_comp)
          end do
       end if
       
       !*********************************************
       ! Compute a provisional D(s) at time n+1

       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(1,:),snew,comp,bc_comp,&
                                              dx,the_bc_tower)
          end do
          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(1,n),diff_coef)
          enddo
       endif
       
       !***********************
       ! Integrate reactions
       
       do n = 1, nlevs
          call multifab_copy_c(snew(n),1,sold(n),1,nscal,sold(n)%ng)
       enddo
       
       call react(mla,the_bc_tower,snew,dx,dt,t,provisional,adv_s,D_s)
       
       if (k > 1) then
          do n = 1, nlevs
             call multifab_copy_c(D_s(2,n),1,D_s(3,n),1,nscal-1)
          enddo
       end if

       !*********************************************
       ! Compute D(s) at time n+1

       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(3,:),snew,comp,bc_comp,&
                                              dx,the_bc_tower)
          end do
          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(3,n),diff_coef)
          enddo
       endif
    
       !*****************************
       ! Compute aprrox to int(A+D+R)

       call intgrl(I_ADR,sold,snew,adv_s,D_s,dt,nlevs)

 !Remove me:
 !*****************************
 ! Check convergence of SDC iters
 !      write(*,*)
 !      write(*,*)
 !      write(*,*) 'SDC corrections:  k=',k
 !      write(*,*) 's_k - s_k+1 min & max'
 !      do n = 1, nlevs
 !         call multifab_copy_c(difference(n),1,snew_old(n),1,nscal)
 !         call multifab_sub_sub_c(difference(n),1,snew(n),1,nscal)
 !         write(*,*)'LEVEL ',n
 !         do comp = 1,nscal
 !            smin = multifab_min_c(difference(n),comp) 
 !            smax = multifab_max_c(difference(n),comp)
 !            write(*,*)'component ',comp, ': ', smin, smax
 !         enddo
 !         write(*,*)
 !         call multifab_copy_c(snew_old(n),1,snew(n),1,nscal)
 !      enddo

       



    end do  ! sdc_iters loop

    ! snew is copied into sold in varden.f90

    deallocate(is_conservative)

    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(source(n))
       call multifab_destroy(I_ADR(n))
       call multifab_destroy(snew_old(n))
       call multifab_destroy(difference(n))
       do i = 1,dm
         call multifab_destroy(sflux(n,i))
         call multifab_destroy(sedge(n,i))
       end do
       do i = 1,n_adv
         call multifab_destroy(adv_s(i-1,n))
       end do
       do i = 1,n_diff
         call multifab_destroy(D_s(i-1,n))
       end do
    enddo

    deallocate(umac_nodal_flag)
    deallocate(scal_force,divu,sflux,sedge)
    deallocate(D_s, adv_s,I_ADR)
    deallocate(snew_old, difference)

2000 format('... level ', i2,' new min/max : density           ',e17.10,2x,e17.10)
2001 format('... level ', i2,' new min/max :  tracer           ',e17.10,2x,e17.10)

  end subroutine scalar_advance_sdc

 

end module scalar_advance_sdc_module
