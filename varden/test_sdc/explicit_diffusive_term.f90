module explicit_diffusive_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use viscous_module
  use macproject_module
  use probin_module, only : stencil_order, verbose, mg_verbose, cg_verbose, &
                            sdc_iters

  implicit none

  private

  public :: get_explicit_diffusive_term

contains

  subroutine get_explicit_diffusive_term(mla,lap_data,data,data_comp,&
                                         bc_comp,dx,the_bc_tower,adj_index)

    use bl_constants_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: lap_data(:)
    type(multifab) , intent(in   ) :: data(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: data_comp
    integer        , intent(in   ) :: bc_comp
    logical, optional, intent(in ) :: adj_index

    ! local variables
    type(multifab) :: alpha(mla%nlevel), beta(mla%nlevel)
    type(multifab) :: phi(mla%nlevel), Lphi(mla%nlevel)
    integer         :: n, comp
    integer         :: nlevs,dm
    logical         :: ladj_index


    if(present(adj_index)) then; ladj_index = adj_index
    else; ladj_index = .false.
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
    enddo

    !***********************************
    ! Compute explicit diffusive_term
    !***********************************
     do n = 1, nlevs
        call multifab_copy_c(phi(n),1,data(n),data_comp,1,1)
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
