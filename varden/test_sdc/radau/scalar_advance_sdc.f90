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
    type(multifab) :: source(mla%nlevel)
    type(multifab) :: I_AD(mla%nlevel,n_interp_pts+1)
    type(multifab), save, allocatable :: I_R(:,:)

    logical        :: is_conservative(nscal)
    logical        :: umac_nodal_flag(mla%dim)
    ! advection term at the various times/quadrtature points:  
    ! adv_s(level,time) = - div(su) 
    type(multifab)  :: adv_s(mla%nlevel,n_adv)
    type(multifab), allocatable  :: adv_rho(:)
    ! diffusion at various times:  D_s(level,time)
    type(multifab)  :: D_s(mla%nlevel,n_diff)
    type(multifab)  :: snew_temp(mla%nlevel)

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

    ! for using time-lagged I_R as a source term in Godunov
    if (iter < 1) then
       allocate(I_R(mla%nlevel,n_interp_pts+1))
       do n = 1, nlevs
          do j = 1, n_interp_pts+1
             call multifab_build( I_R(n,j),mla%la(n),nspec,0)
             call setval(I_R(n,j),zero,all=.true.)
          end do
       end do
    endif

    do n = 1, nlevs
       call multifab_build( scal_force(n),ext_scal_force(n)%la,nscal,1)
       call multifab_build( divu(n),mla%la(n),    1,1)
       call multifab_build( source(n),mla%la(n),nspec,0)
       do j = 1, n_diff
          call multifab_build( D_s(n,j),mla%la(n),nspec,0)
          call setval(D_s(n,j),zero,all=.true.)
       enddo
       do j = 1, n_adv
          call multifab_build( adv_s(n,j),mla%la(n),nspec,0)
          call setval(adv_s(n,j),zero,all=.true.)
       enddo

       do j = 1, n_interp_pts+1
          call multifab_build( I_AD(n,j),mla%la(n),nspec,0)
          call setval(I_AD(n,j),zero,all=.true.)         
       enddo

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

       call multifab_build_copy(snew_temp(n),sold(n))

       call multifab_build(snew_old(n),mla%la(n),  nscal,0)
       call multifab_build(difference(n),mla%la(n),nscal,0)
       call copy(snew_old(n), sold(n), 0)

       call multifab_copy(snew(n),sold(n))

    enddo
  

!-------------------------------------------------------------------
! Compute provisional solution using Godunov scheme for advection 
! and implicit eulers for diffusion
!-------------------------------------------------------------------
    !************************************************************
    ! Create scalar force at time n.
    !************************************************************
    ! could use D_n and/or time lagged I_R as source
    diff_fac = zero     !in mksource: diff_fac*source
    call mksource(nlevs,scal_force,ext_scal_force,source,diff_fac) 

    !***********************************
    ! Create edge state scalars/fluxes
    !***********************************
    call mkflux(mla,sold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                the_bc_tower%bc_tower_array,is_vel,is_conservative)



    do j = 1, n_interp_pts  !loop over quadrature points

       ! For computing the integral of the reactions
       do n = 1, nlevs
          call multifab_copy_c(I_R(n,j),1,snew(n),2,nspec)        
       enddo

       !*****************************************************************
       ! Update the scalars with conservative or convective differencing.
       !*****************************************************************
       call update(mla,snew_temp,umac,sedge,sflux,scal_force,snew_temp,&
                   adv_s(:,1),dx,&
                   dt*sdc_dt_fac(j),t,is_vel,is_conservative,&
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

       !***************************************
       ! Solve u_t = A + D using implicit euler
       !***************************************       
       if (diff_coef > ZERO) then
          do n = 1,nlevs
             call setval(source(n),0.0_dp_t,all=.true.)
          end do

          do  comp = 2, nscal
             bc_comp = dm+comp
             visc_mu = dt*sdc_dt_fac(j)*diff_coef
             call diff_scalar_solve(mla,snew_temp,source,dx,t,visc_mu,&
                             the_bc_tower,comp,bc_comp,adj_index=.true.)
          end do
       end if

       !***************************************
       ! Compute a provisional D_j+1 for reactions
       !***************************************
       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(:,j),snew_temp,comp,&
                  bc_comp,dx,the_bc_tower,adj_index=.true.,is_vel=.false.)
          end do
          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(n,j),diff_coef)
          enddo
       endif

       !************************
       ! Integrate entire system
       !************************
       call react(mla,the_bc_tower,snew,dx,dt,t,j,adv_s,adv_rho,D_s,&
                  sdc_flag=1)
       ! use sdc_flag=1 for provisional, =2 for SDC

       !**************************************
       ! Compute the integral of the reactions
       !**************************************      
       call mk_provis_I_R(j,snew,dt*sdc_dt_fac(j))

       !******************************
       ! Recompute D_j+1 from snew_j+1
       !******************************
       if (diff_coef .gt. ZERO) then
          do comp = 2, nscal
             bc_comp = dm+comp
             call get_explicit_diffusive_term(mla,D_s(:,j),snew,comp,&
                                bc_comp,dx,the_bc_tower,adj_index=.true.,is_vel=.false.)
          end do
          do n = 1, nlevs
             call multifab_mult_mult_s(D_s(n,j),diff_coef)
          enddo
       endif

       if (mass_fractions) then 
          do n =1, nlevs
             call multifab_copy_c(snew_temp(n),2,snew(n),2,nspec)
          end do
       else 
          do n = 1, nlevs !rho not integrated w/ rxns, give snew updated rho
             call multifab_copy_c(snew_temp(n),2,snew(n),2,nspec)
             call multifab_copy_c(snew(n),1,snew_temp(n),1,1)
          end do
       endif

    enddo  ! do j = 1, n_nterp_pts
    
    !*****************************
    ! Compute aprrox to int(A+D)
    !*****************************
    call mk_I_AD(I_AD,adv_s,D_s,nlevs)

    !******************************
    ! Check changes over a timestep
    !******************************
    if (sdc_verbose > 0) then
       if (parallel_IOProcessor()) write(*,*)  'Initial/Provisional change in s'
       if (parallel_IOProcessor()) write(*,*) '(s_new^0 - s_old) min & max'
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

!----------------------------------------------------------------
! SDC iterations
!----------------------------------------------------------------
    do k = 1,sdc_iters
       
       !**********************************
       ! Create new forcing term at time n
       !********************************** 
       call mk_sdc_adv_source(scal_force,n_interp_pts+1)

       !*********************************
       ! Create edge state scalars/fluxes
       !*********************************
       call mkflux(mla,sold,sedge,sflux,umac,scal_force,divu,dx,dt, &
                   the_bc_tower%bc_tower_array,is_vel,is_conservative)
              
       do n = 1, nlevs
          call multifab_copy(snew_temp(n),sold(n))
          call multifab_copy(snew(n),sold(n))
          call setval(I_R(n,n_interp_pts+1),ZERO,all=.true.)
       end do


       do j = 1, n_interp_pts ! loop over quadrature points        

          !*****************************
          ! Compute aprrox to int(A+D)
          !*****************************
          call mk_sdc_adv_source(scal_force,j)

          ! For computing the integral of the reactions
          do n = 1, nlevs
             call multifab_copy_c(I_R(n,j),1,snew(n),2,nspec)
          enddo

          !****************************************************************
          ! Update the scalars with conservative or convective differencing
          !****************************************************************
          call update(mla,snew_temp,umac,sedge,sflux,scal_force,snew_temp,&
                      adv_s(:,1),dx,dt*sdc_dt_fac(j),t,&
                      is_vel,is_conservative,the_bc_tower%bc_tower_array,&
                      adv_rho)

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
          ! Update Diffusion:
          ! Solve snew = sold +dt*(delta A(s) + delta D(s)) + I_ADR
          !********************************************************
          do n = 1, nlevs
             call multifab_copy_c(source(n),1,D_s(n,j),1,nspec)
             call multifab_mult_mult_s(source(n),-ONE/diff_coef)
          enddo

          if (diff_coef > ZERO) then
             do  comp = 2, nscal
                bc_comp = dm+comp
                visc_mu = sdc_dt_fac(j)*dt*diff_coef
                call diff_scalar_solve(mla,snew_temp,source,dx,t,visc_mu,&
                                       the_bc_tower,comp,bc_comp,adj_index=.true.)
             end do
          end if

          !******************************************
          ! Compute a provisional D_j+1 for reactions
          !******************************************
          if (diff_coef .gt. ZERO) then
             do comp = 2, nscal
                bc_comp = dm+comp
                call get_explicit_diffusive_term(mla,D_s(:,j+n_interp_pts),&
                     snew_temp,comp,bc_comp,dx,the_bc_tower,adj_index=.true.,is_vel=.false.)
             end do
             do n = 1, nlevs
                call multifab_mult_mult_s(D_s(n,j+n_interp_pts),diff_coef)
             enddo
          endif

          !************************
          ! Integrate entire system
          !************************
          call react(mla,the_bc_tower,snew,dx,dt,t,j,adv_s,adv_rho,D_s,&
                     sdc_flag=2)
          ! use sdc_flag=1 for provisional, =2 for SDC

          !**************************************
          ! Compute the integral of the reactions
          !**************************************      
          call mk_sdc_I_R(j,snew,dt*sdc_dt_fac(j))

          !****************************
          ! Recompute a D_j+1 from snew
          !****************************
          if (diff_coef .gt. ZERO) then
             do comp = 2, nscal
                bc_comp = dm+comp
                call get_explicit_diffusive_term(mla,D_s(:,j+n_interp_pts),&
                     snew,comp,bc_comp,dx,the_bc_tower,adj_index=.true.,is_vel=.false.)
             end do
             do n = 1, nlevs
                call multifab_mult_mult_s(D_s(n,j+n_interp_pts),diff_coef)
             enddo
          endif

          if (mass_fractions) then
             do n = 1, nlevs
                call multifab_copy_c(snew_temp(n),2,snew(n),2,nspec)
             end do
          else 
             do n = 1, nlevs !rho not integrated w/ rxns, give snew updated rho
                call multifab_copy_c(snew_temp(n),2,snew(n),2,nspec)
                call multifab_copy_c(snew(n),1,snew_temp(n),1,1)
             end do
          endif

       enddo  ! do j = 1, n_nterp_pts

       if (k .NE. sdc_iters) then
          ! move new values into old 
          do j = 1, n_interp_pts
             do n = 1, nlevs
                call multifab_copy(D_s(n,j),D_s(n,j+n_interp_pts))
             enddo
          end do
          
          !*****************************
          ! Compute aprrox to int(A+D)
          !*****************************
          call mk_I_AD(I_AD,adv_s,D_s,nlevs)         
       end if
       
       !*******************************
       ! Check convergence of SDC iters
       !*******************************
       if (sdc_verbose > 0) then
          if (parallel_IOProcessor()) write(*,*) 'SDC corrections:  k=',k
          if (parallel_IOProcessor()) write(*,*) '(s_k - s_k+1) min & max'
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
       end if
          
    end do  ! sdc_iters loop

1000 FORMAT(5(E15.8,1X))
1001 FORMAT(I5,4(E15.8,1X))

    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(source(n))
       call multifab_destroy(snew_temp(n))
       call multifab_destroy(snew_old(n))
       call multifab_destroy(difference(n))
       do i = 1,dm
         call multifab_destroy(sflux(n,i))
         call multifab_destroy(sedge(n,i))
       end do
       do i = 1,n_adv
         call multifab_destroy(adv_s(n,i))
       end do
       do i = 1,n_diff
         call multifab_destroy(D_s(n,i))
       end do
       do i = 1, n_interp_pts+1
         call multifab_destroy(I_AD(n,i))          
!         call multifab_destroy(I_R(n,i))          
      enddo
   end do

   if (mass_fractions) then
      do n = 1,nlevs
         call multifab_destroy(adv_rho(n))
      end do
   end if

2000 format('... level ', i2,' min/max diff : density      ',e17.10,2x,e17.10)
2001 format('... level ', i2,' min/max diff :  tracer ',i2,'   ',e17.10,2x,e17.10)    

  contains
  
   subroutine write_plotfile(istep_to_write, n_plot_comps, mf)

    integer,         intent(in   ) :: istep_to_write
    integer,         intent(in   ) :: n_plot_comps
    type(multifab),  intent(in   ) :: mf(:)


    type(multifab), allocatable :: plotdata(:) 
    integer        :: n

 
    allocate(plotdata(nlevs))
 
    do n = 1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy(plotdata(n), mf(n))
    end do
        
    write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), t, dx(1,:))

    do n = 1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)
    
  end subroutine write_plotfile

  subroutine mk_provis_I_R(j,s,dt_sdc)
    ! 
    ! This is really I_R/(sdc_fac*dt) plus Adv terms from I_AD, 
    ! and -Aold from delta(A) that appears in the s update
    !
    integer,        intent(in   ) :: j
    type(multifab), intent(in   ) :: s(:)
    real(dp_t),     intent(in   ) :: dt_sdc

     ! local
    integer                 :: n,i,ng,comp
    integer                 :: ix,iy,iz
    integer                 :: lo(s(1)%dim),hi(s(1)%dim)
    real(dp_t), pointer     :: sop(:,:,:,:)
    type(dpt_pntr) :: A(n_adv), D(n_diff), IR(n_interp_pts+1) 


    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle      
          call get_interp_pts_array(adv_s,A,n,i,n_adv)
          call get_interp_pts_array(D_s,D,n,i,n_diff)
          call get_interp_pts_array(I_R,IR,n,i,n_interp_pts+1)
          sop => dataptr(s(n), i)         
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)

                   do comp = 1, nspec
                      IR(j)%p(ix,iy,1,comp) = &
                        (sop(ix,iy,1,comp+1) - IR(j)%p(ix,iy,1,comp))/dt_sdc &
                        - (D(j)%p(ix,iy,1,comp) + A(1)%p(ix,iy,1,comp))
                      IR(n_interp_pts+1)%p(ix,iy,1,comp) = &
                           IR(n_interp_pts+1)%p(ix,iy,1,comp) + &
                           sdc_dt_fac(j)*IR(j)%p(ix,iy,1,comp) 
                   enddo
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 1, nspec
                         IR(j)%p(ix,iy,iz,comp) = &
                          (sop(ix,iy,iz,comp+1)-IR(j)%p(ix,iy,iz,comp))/dt_sdc &
                          - (D(j)%p(ix,iy,iz,comp) + A(1)%p(ix,iy,iz,comp))
                         IR(n_interp_pts+1)%p(ix,iy,iz,comp) = &
                              IR(n_interp_pts+1)%p(ix,iy,iz,comp) + &
                              sdc_dt_fac(j)*IR(j)%p(ix,iy,iz,comp) 
                      enddo

                   end do
                end do
             end do

          end select
       end do
    end do

  end subroutine mk_provis_I_R

  subroutine mk_sdc_I_R(j,s,dt_sdc)
    ! 
    ! This is really I_R/(sdc_fac*dt) plus Adv terms from I_AD, 
    ! and -Aold from delta(A) appearing in the s update
    !
    integer,        intent(in   ) :: j
    type(multifab), intent(in   ) :: s(:)
    real(dp_t),     intent(in   ) :: dt_sdc        

     ! local
    integer                 :: n,i,ng,comp
    integer                 :: ix,iy,iz
    integer                 :: lo(s(1)%dim),hi(s(1)%dim)
    real(dp_t), pointer     :: sop(:,:,:,:)
    type(dpt_pntr) :: A(n_adv),D(n_diff),IR(n_interp_pts+1),IAD(n_interp_pts+1) 


    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle      
          call get_interp_pts_array(adv_s,A,n,i,n_adv)
          call get_interp_pts_array(D_s,D,n,i,n_diff)
          call get_interp_pts_array(I_R,IR,n,i,n_interp_pts+1)
          call get_interp_pts_array(I_AD,IAD,n,i,n_interp_pts+1)
          sop => dataptr(s(n), i)         
          lo = lwb(get_box(s(n), i))
          hi = upb(get_box(s(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)

                   do comp = 1, nspec
                      IR(j)%p(ix,iy,1,comp) = &
                           (sop(ix,iy,1,comp+1)-IR(j)%p(ix,iy,1,comp))/dt_sdc &
                           - IAD(j)%p(ix,iy,1,comp) - (A(1)%p(ix,iy,1,comp) &
                           + D(j+n_interp_pts)%p(ix,iy,1,comp) &
                           - D(j)%p(ix,iy,1,comp))
                      IR(n_interp_pts+1)%p(ix,iy,1,comp) = &
                           IR(n_interp_pts+1)%p(ix,iy,1,comp) + &
                           sdc_dt_fac(j)*IR(j)%p(ix,iy,1,comp) 
                   enddo
                end do
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 1, nspec
                         IR(j)%p(ix,iy,iz,comp) = &
                           (sop(ix,iy,iz,comp+1)-IR(j)%p(ix,iy,iz,comp))/dt_sdc&
                            - IAD(j)%p(ix,iy,iz,comp) - (A(1)%p(ix,iy,iz,comp)&
                            + D(j+n_interp_pts)%p(ix,iy,iz,comp) &
                            - D(j)%p(ix,iy,iz,comp))
                         IR(n_interp_pts+1)%p(ix,iy,iz,comp) = &
                              IR(n_interp_pts+1)%p(ix,iy,iz,comp) + &
                              sdc_dt_fac(j)*IR(j)%p(ix,iy,iz,comp) 
                      enddo

                   end do
                end do

             end do

          end select
       end do
    end do   

  end subroutine mk_sdc_I_R

  subroutine mk_sdc_adv_source(src,j)
    type(multifab), intent(inout) :: src(:)
    integer,        intent(in   ) :: j      

    ! local
    integer             :: n,i,ng,comp
    integer             :: ix,iy,iz
    integer             :: lo(src(1)%dim),hi(src(1)%dim)
    real(dp_t), pointer :: sop(:,:,:,:),esfop(:,:,:,:)
    type(dpt_pntr)      :: IR(n_interp_pts+1),IAD(n_interp_pts+1),A(n_adv) 


    ng = src(1)%ng

    do n=1,nlevs
       do i = 1, src(n)%nboxes
          if ( multifab_remote(src(n), i) ) cycle      
          call get_interp_pts_array(I_R,IR,n,i,n_interp_pts+1)
          call get_interp_pts_array(I_AD,IAD,n,i,n_interp_pts+1)
          call get_interp_pts_array(adv_s,A,n,i,n_adv)
          sop    => dataptr(src(n), i)         
          esfop  => dataptr(ext_scal_force(n), i)         
          lo = lwb(get_box(src(n), i))
          hi = upb(get_box(src(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)
                   ! scal-force was set = ext_scal-force in the beginning
                   ! so sop(comp = 1) remains unchanged
                   do comp = 1, nspec
                      sop(ix,iy,1,comp+1) = IAD(j)%p(ix,iy,1,comp) + &
                           IR(j)%p(ix,iy,1,comp)&! - A(1)%p(ix,iy,1,comp) &
                           + esfop(ix,iy,1,comp+1)
                   enddo
                end do
             end do

             do comp = 2,nscal
             ! use 0th order extrapolation for ghost cells
                do ix = lo(1), hi(1)
                   iy = lo(2)-1
                   sop(ix,iy,1,comp) = sop(ix,iy+1,1,comp)
                   iy = hi(2)+1
                   sop(ix,iy,1,comp) = sop(ix,iy-1,1,comp)
                enddo            
                do iy = lo(2), hi(2)
                   ix = lo(1)-1
                   sop(ix,iy,1,comp) = sop(ix+1,iy,1,comp)
                   ix = hi(1)+1
                   sop(ix,iy,1,comp) = sop(ix-1,iy,1,comp)
                enddo
             end do

          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 1, nspec
                         sop(ix,iy,iz,comp+1) = IAD(j)%p(ix,iy,iz,comp) + &
                           IR(j)%p(ix,iy,iz,comp) + esfop(ix,iy,iz,comp+1)
                      enddo

                   end do
                end do
             end do

             do comp = 2,nscal
             ! use 0th order extrapolation for ghost cells
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
                      iz = lo(3)-1
                      sop(ix,iy,iz,comp) = sop(ix,iy,iz+1,comp)
                      iz = hi(3)+1
                      sop(ix,iy,iz,comp) = sop(ix,iy,iz-1,comp)
                   enddo
                enddo
                do iy = lo(2), hi(2)
                   do iz = lo(3),hi(3)
                      ix = lo(1)-1
                      sop(ix,iy,iz,comp) = sop(ix+1,iy,iz,comp)
                      ix = hi(1)+1
                      sop(ix,iy,iz,comp) = sop(ix-1,iy,iz,comp)
                   enddo
                enddo
                do ix = lo(1), hi(1)
                   do iz = lo(3), hi(3)
                      iy = lo(2)-1
                      sop(ix,iy,iz,comp) = sop(ix,iy+1,iz,comp)
                      iy = hi(2)+1
                      sop(ix,iy,iz,comp) = sop(ix,iy-1,iz,comp)
                   enddo
                enddo
             end do

          end select
       end do

       call multifab_fill_boundary(src(n))
    end do ! do n = 1, nlevs
 
  end subroutine mk_sdc_adv_source

  subroutine mk_sdc_diff_source(src,j)
    type(multifab), intent(inout) :: src(:)
    integer,        intent(in   ) :: j

    ! local
    integer             :: n,i,ng,comp
    integer             :: ix,iy,iz
    integer             :: lo(src(1)%dim),hi(src(1)%dim)
    real(dp_t), pointer :: sop(:,:,:,:),esfop(:,:,:,:)
    type(dpt_pntr) :: A(n_adv), D(n_diff), IR(n_interp_pts),IAD(n_interp_pts) 


    ng = src(1)%ng

    do n=1,nlevs
       do i = 1, src(n)%nboxes
          if ( multifab_remote(src(n), i) ) cycle      
          call get_interp_pts_array(adv_s,A,n,i,n_adv)
          call get_interp_pts_array(D_s,D,n,i,n_diff)
          call get_interp_pts_array(I_R,IR,n,i,n_interp_pts)
          call get_interp_pts_array(I_AD,IAD,n,i,n_interp_pts)
          sop    => dataptr(src(n), i)         
          esfop  => dataptr(ext_scal_force(n), i)         
          lo = lwb(get_box(src(n), i))
          hi = upb(get_box(src(n), i))
          select case (dm)
          case (2)
             do iy = lo(2), hi(2)          
                do ix = lo(1),hi(1)

                   do comp = 1, nspec
                      sop(ix,iy,1,comp) = (A(1)%p(ix,iy,1,comp) &
                           - D(j)%p(ix,ix,1,comp) + IAD(j)%p(ix,iy,1,comp) &
                           + IR(j)%p(ix,iy,1,comp) &
                           + esfop(ix,iy,1,comp+1))/diff_coef
                      ! dividing by diff_coef here because diff_scalar_solve()
                      ! will mult src by dt*diff_coef
                   enddo
                end do
             end do
             
          case (3)
             do iz = lo(3), hi(3)
                do iy = lo(2), hi(2)
                   do ix = lo(1), hi(1)
    
                      do comp = 1, nspec
                         sop(ix,iy,iz,comp) = (A(1)%p(ix,iy,iz,comp) &
                              - D(j)%p(ix,ix,iz,comp)+IAD(j)%p(ix,iy,iz,comp)&
                              + IR(j)%p(ix,iy,iz,comp) &
                              + esfop(ix,iy,iz,comp+1))/diff_coef
                         ! dividing by diff_coef here because diff_scalar_solve()
                         ! will mult src by dt*diff_coef
                      enddo

                   end do
                end do
             end do

          end select
       end do
    end do
  end subroutine mk_sdc_diff_source

  end subroutine scalar_advance_sdc


end module scalar_advance_sdc_module
