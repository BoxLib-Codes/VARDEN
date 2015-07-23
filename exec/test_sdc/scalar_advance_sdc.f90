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
  use fabio_module
  use probin_module, only : nscal, diff_coef, diffusion_type, stencil_order, &
                            verbose, mg_verbose, cg_verbose, sdc_iters, &
                            mass_fractions, nspec, sdc_verbose

  implicit none

  public :: scalar_advance_sdc


contains

!************************************************************
! This does SDC with Gauss-Lobatto quadrature nodes 
! (i.e. both interval endpoints are used)
!************************************************************
  subroutine scalar_advance_sdc(mla,uold,sold,snew,umac, &
                                ext_scal_force,dx,dt,t,the_bc_tower)

    use mkflux_module
    use mkforce_module
    use update_module
    use bl_constants_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: ext_scal_force(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,t
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    type(multifab) :: scal_force(mla%nlevel), divu(mla%nlevel)
    type(multifab) :: sflux(mla%nlevel,mla%dim), sedge(mla%nlevel,mla%dim)
    type(multifab) :: source(mla%nlevel),I_ADR(mla%nlevel)
    logical        :: is_conservative(nscal)
    logical        :: umac_nodal_flag(mla%dim)
    ! these values are used to store the vlaue of the advective term and 
    ! the diffusive term at the quadrature nodes.
    ! D_s(level,quad point)
    type(multifab)  :: D_s(mla%nlevel,0:n_diff-1) 
    type(multifab)  :: adv_s(mla%nlevel,0:n_adv-1)
    type(multifab), allocatable  :: adv_rho(:)
 
    integer         :: i,j,k,n
    integer         :: comp,bc_comp,nlevs,dm
    logical         :: is_vel
    real(kind=dp_t) :: diff_fac,visc_mu
    real(kind=dp_t) :: smin,smax,norm

    type(multifab)  :: snew_old(mla%nlevel), difference(mla%nlevel)


    nlevs = mla%nlevel
    dm    = mla%dim

    if (mass_fractions) then
       is_conservative(:) = .true.
       allocate(adv_rho(mla%nlevel))
       do n = 1,nlevs
          call multifab_build(adv_rho(n),mla%la(n),1,0)
          call setval(adv_rho(n),zero,all=.true.)
       end do
    else
       is_conservative(1) = .true.
       do comp = 2, nscal
          is_conservative(comp) = .false.
       end do
    end if

    is_vel  = .false.

    do n = 1, nlevs
       call multifab_build( scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build( divu(n),mla%la(n),    1,1)
       call multifab_build( source(n),mla%la(n),nspec,0)
       do j = 0, n_diff-1
          call multifab_build( D_s(n,j),mla%la(n),nspec,0)
          call setval(D_s(n,j),zero,all=.true.)
       enddo
       do j = 0, n_adv-1
          call multifab_build( adv_s(n,j),mla%la(n),nspec,0)
          call setval(adv_s(n,j),zero,all=.true.)
       enddo
       call multifab_build( I_ADR(n),mla%la(n),nspec,0)
       call setval(I_ADR(n),zero,all=.true.)

       do i = 1,dm
         umac_nodal_flag(:) = .false.
         umac_nodal_flag(i) = .true.
         call multifab_build(sflux(n,i),scal_force(n)%la,nscal,0,&
                             nodal = umac_nodal_flag)
         call multifab_build(sedge(n,i),scal_force(n)%la,nscal,0,&
                             nodal = umac_nodal_flag)
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
  

!---------------------------------------------------------------
! Compute a provisional solution at t=n+1 using Godunov for the
! advection term and implicit euler for diffusion
!--------------------------------------------------------------- 
    !*********************************************
    ! Compute D(s) = laplacian(s) at time n
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(:,0),sold,comp,bc_comp,dx,&
                                           the_bc_tower,adj_index=.true.,&
                                           is_vel=.false.)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(n,0),diff_coef)
       enddo
    endif


    !************************************************************
    ! Create scalar force at time n.
    ! Provided as a source term to calculate advection term adv_s
    !************************************************************
    ! could use D_s(n,0) and/or time lagged I_R as a source
    diff_fac = zero     !in mksource: diff_fac*source
    call mksource(nlevs,scal_force,ext_scal_force,source,diff_fac) 

    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************
    call mkflux(mla,sold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_tower%bc_tower_array,is_vel,is_conservative)

    !*****************************************************************
    ! Update the scalars with conservative or convective differencing.
    !*****************************************************************
    call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s(:,0),dx,dt,t,& 
                is_vel,is_conservative,the_bc_tower%bc_tower_array,adv_rho)


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
    if (diff_coef > ZERO) then
       do n = 1,nlevs
          ! just a dummy to pass in zero to diffusion solve
          call setval(source(n),zero,all=.true.)
       end do

       do  comp = 2, nscal
           bc_comp = dm+comp
           visc_mu = dt*diff_coef
           call diff_scalar_solve(mla,snew,source,dx,t,visc_mu,&
                                  the_bc_tower,comp,bc_comp,adj_index=.true.)
        end do
    end if

    !***********************************************
    ! Compute a provisional Diffusion(s) at time n+1
    !***********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(:,1),snew,comp,bc_comp,&
                                           dx,the_bc_tower,adj_index=.true.,&
                                           is_vel=.false.)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(n,1),diff_coef)
       enddo
    endif

    !***********************
    ! Integrate reactions
    !***********************
     if(mass_fractions) then !rho gets integrated
       do n = 1, nlevs
          call multifab_copy_c(snew(n),1,sold(n),1,nscal,sold(n)%ng)
       enddo
    else !rho doesn't get integrated; leave updated rho
       do n = 1, nlevs
          call multifab_copy_c(snew(n),2,sold(n),2,nspec,sold(n)%ng)
       enddo
    end if

!    call multifab_print(snew(1))
    ! use sdc_flag=1 for provisional, =2 for SDC
    call react(mla,the_bc_tower,snew,dx,dt,t,adv_s,adv_rho,D_s,sdc_flag=1)

    !*********************************************
    ! Compute D(s) at time n+1
    !*********************************************
    if (diff_coef .gt. ZERO) then
       do comp = 2, nscal
          bc_comp = dm+comp
          call get_explicit_diffusive_term(mla,D_s(:,2),snew,comp,bc_comp,&
                                           dx,the_bc_tower,adj_index=.true.,&
                                           is_vel=.false.)
       end do
       do n = 1, nlevs
          call multifab_mult_mult_s(D_s(n,2),diff_coef)
       enddo
    endif
    
    !****************************************************
    ! Use quadrature to compute approx to integral(A+D+R)
    !****************************************************
    call provisional_intgrl(I_ADR,sold,snew,adv_s,D_s,dt,nlevs)



    !******************************
    ! Check changes over a timestep
    !******************************
    if (sdc_verbose > 0) then
       if (parallel_IOProcessor()) write(*,*)  'Initial/Provisional change in s'
       if (parallel_IOProcessor()) write(*,*) 's_new^0 - s_old min & max'
       do n = 1, nlevs
          call multifab_copy_c(difference(n),1,snew_old(n),1,nscal)
          call multifab_sub_sub_c(difference(n),1,snew(n),1,nscal)
          norm =  multifab_norm_l2(difference(n))*sqrt(dx(1,1))
          if (parallel_IOProcessor()) write(*,*) 'L2 norm',norm
          do comp = 1,nscal
             smin = multifab_min_c(difference(n),comp) 
             smax = multifab_max_c(difference(n),comp)
             if (comp .eq. 1) then
                if (parallel_IOProcessor()) write(6,2000) n,smin,smax
             else 
                if (parallel_IOProcessor()) write(6,2001) n,comp,smin,smax
             end if
          enddo
          call multifab_copy_c(snew_old(n),1,snew(n),1,nscal)
       enddo
    end if

!-----------------------------------------------------------------------
! SDC iterations
!-----------------------------------------------------------------------
    do k = 1,sdc_iters
      
       !**********************************
       ! Create new forcing term at time n
       !**********************************
       diff_fac = ONE
       call mksource(nlevs,scal_force,ext_scal_force,I_ADR,diff_fac)

       !*********************************
       ! Create edge state scalars/fluxes
       !*********************************
       call mkflux(mla,sold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                   the_bc_tower%bc_tower_array,is_vel,is_conservative)
       
       !****************************************************************
       ! Update the scalars with conservative or convective differencing
       !****************************************************************
       call update(mla,sold,umac,sedge,sflux,scal_force,snew,adv_s(:,0),&
                   dx,dt,t,is_vel,is_conservative,&
                   the_bc_tower%bc_tower_array,adv_rho)

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
       
       !********************************************************
       ! Update Diffusion
       ! Solve snew = sold +dt*(delta A(s) + delta D(s)) + I_ADR
       !********************************************************
       do n = 1, nlevs
          call multifab_copy_c(source(n),1,D_s(n,2),1,nscal-1)
          call multifab_mult_mult_s(source(n),-ONE/diff_coef)
       enddo

       if (diff_coef > ZERO) then
          do  comp = 2, nscal
             bc_comp = dm+comp
             visc_mu = dt*diff_coef
             call diff_scalar_solve(mla,snew,source,dx,t,visc_mu,&
                                    the_bc_tower,comp,bc_comp,adj_index=.true.)
          end do
       end if
       
       !********************************
       ! Compute approx D(s) at time n+1
       !********************************
       do n = 1,nlevs
          call setval(D_s(n,1),zero,all=.true.)
       end do

       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(:,1),snew,comp,bc_comp,&
                                              dx,the_bc_tower,adj_index=.true.,&
                                              is_vel=.false.)
          end do

          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(n,1),diff_coef)
          enddo
       endif
       
       !********************
       ! Integrate reactions
       !********************
       if(mass_fractions) then !rho gets integrated too
          do n = 1, nlevs
             call multifab_copy_c(snew(n),1,sold(n),1,nscal,sold(n)%ng)
          enddo
       else !rho doesn't get integrated; leave updated rho
          do n = 1, nlevs
             call multifab_copy_c(snew(n),2,sold(n),2,nspec,sold(n)%ng)
          enddo
       end if
       ! use sdc_flag=2 for SDC; =1 for provisional
       call react(mla,the_bc_tower,snew,dx,dt,t,adv_s,adv_rho,D_s,sdc_flag=2)

       !*********************************
       ! Compute updated D(s) at time n+1
       !*********************************
       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(:,3),snew,comp,bc_comp,&
                                              dx,the_bc_tower,adj_index=.true.,&
                                              is_vel=.false.)
          end do
          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(n,3),diff_coef)
          enddo
       endif
 
       !****************************************
       ! Compute quadrature aprrox to int(A+D+R)
       !****************************************
       call intgrl(I_ADR,sold,snew,adv_s,D_s,dt,nlevs)


       do n = 1, nlevs
          call multifab_copy_c(D_s(n,2),1,D_s(n,3),1,nscal-1)
       enddo
       

       !*******************************
       ! Check convergence of SDC iters
       !*******************************
       if (sdc_verbose >0) then
          if (parallel_IOProcessor()) write(*,*) 'SDC corrections:  k=',k
          if (parallel_IOProcessor()) write(*,*) 's_k - s_k+1 min & max'
          do n = 1, nlevs
             call multifab_copy_c(difference(n),1,snew_old(n),1,nscal)
             call multifab_sub_sub_c(difference(n),1,snew(n),1,nscal)
             norm = multifab_norm_l2(difference(n))*sqrt(dx(1,1))
             if (parallel_IOProcessor()) write(*,*) 'L2 norm', norm
             do comp = 1,nscal
                smin = multifab_min_c(difference(n),comp) 
                smax = multifab_max_c(difference(n),comp)
                if (comp .eq. 1) then
                   if (parallel_IOProcessor()) write(6,2000) n,smin,smax
                else 
                   if (parallel_IOProcessor()) write(6,2001) n,comp,smin,smax
                end if
             enddo
             call multifab_copy_c(snew_old(n),1,snew(n),1,nscal)
          enddo
          if (k .eq. sdc_iters) then
             if (parallel_IOProcessor()) write(*,*)
          end if
       end if

    end do  ! sdc_iters loop

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
         call multifab_destroy(adv_s(n,i-1))
       end do
       do i = 1,n_diff
         call multifab_destroy(D_s(n,i-1))
       end do
    enddo

   if (mass_fractions) then
      do n = 1,nlevs
         call multifab_destroy(adv_rho(n))
      end do
   end if

2000 format('... level ', i2,' new min/max : density      ',e17.10,2x,e17.10)
2001 format('... level ', i2,' new min/max :  tracer ',i2,'   ',e17.10,2x,e17.10)

  end subroutine scalar_advance_sdc



end module scalar_advance_sdc_module
