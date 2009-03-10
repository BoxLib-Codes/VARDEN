module explicit_diffusive_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use viscous_module
  use macproject_module
  use probin_module, only : stencil_order, verbose, mg_verbose, cg_verbose, &
                            sdc_iters, nspec

  implicit none

  private

  public :: get_explicit_diffusive_term

contains

  subroutine get_explicit_diffusive_term(mla,lap_data,data,data_comp,&
                                         bc_comp,dx,the_bc_tower,adj_index,&
                                         is_vel)

    use bl_constants_module
    use probin_module, only: mass_fractions
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: lap_data(:)
    type(multifab) , intent(in   ) :: data(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: data_comp
    integer        , intent(in   ) :: bc_comp
    logical, optional, intent(in ) :: adj_index
    logical, optional, intent(in ) :: is_vel

    ! local variables
    type(multifab) :: alpha(mla%nlevel), beta(mla%nlevel)
    type(multifab) :: phi(mla%nlevel), Lphi(mla%nlevel)
    integer         :: n, comp,i
    integer         :: nlevs,dm
    logical         :: ladj_index
    logical         :: lis_vel


    if(present(adj_index)) then; ladj_index = adj_index
    else; ladj_index = .false.
    endif
    if(present(is_vel)) then; lis_vel = is_vel
    else; lis_vel = .true.
    endif

    nlevs = mla%nlevel
    dm    = mla%dim

    !***********************************
    ! build temporary arrays
    !***********************************
       
    do n = 1, nlevs
       call multifab_build(  phi(n),mla%la(n),    1,1)
       call multifab_build( Lphi(n),mla%la(n),    1,1)
       call multifab_build(alpha(n),mla%la(n),    1,1)
       call multifab_build( beta(n),mla%la(n),   dm,1)
       call setval( phi(n),0.0_dp_t,all=.true.)
       call setval(Lphi(n),0.0_dp_t,all=.true.)
       call setval(alpha(n),ZERO, all=.true.)
       call setval(beta(n),-ONE, all=.true.)

       if (mass_fractions .AND. (.NOT. lis_vel)) then
          do i = 1,dm
             ! mult by density
             call multifab_mult_mult_c(beta(n),i,data(n),1,1,1)
          end do
       end if
    enddo
       
    !***********************************
    ! Compute explicit diffusive_term
    !***********************************
     do n = 1, nlevs
        call multifab_copy_c(phi(n),1,data(n),data_comp,1,1)
        if (mass_fractions .AND. (.NOT.lis_vel)) then
           ! divide out density
           call multifab_div_div_c(phi(n),1,data(n),1,1,1)
        end if
     enddo

     call mac_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,&
                      bc_comp,stencil_order,mla%mba%rr,mg_verbose,cg_verbose)

     if (ladj_index) then
        comp = data_comp-1
     else 
        comp = data_comp
     endif

     do n = 1, nlevs
        call multifab_copy_c(lap_data(n),comp,Lphi(n),1)
     enddo

     do n = 1,nlevs
        call multifab_destroy(alpha(n))
        call multifab_destroy( beta(n))
        call multifab_destroy(  phi(n))
        call multifab_destroy( Lphi(n))
     enddo
 
   end subroutine get_explicit_diffusive_term

end module explicit_diffusive_module
