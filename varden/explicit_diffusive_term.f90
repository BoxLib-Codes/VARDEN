module explicit_diffusive_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: get_explicit_diffusive_term

contains

  subroutine get_explicit_diffusive_term(mla,lap_data,data,data_comp,bc_comp,dx,&
                                         the_bc_tower)

    use cc_applyop_module
    use probin_module     , only : stencil_order, verbose, mg_verbose, cg_verbose

    use bl_constants_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: lap_data(:)
    type(multifab) , intent(in   ) :: data(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: data_comp
    integer        , intent(in   ) :: bc_comp

    ! local variables
    integer        :: n,d,nlevs,dm
    type(multifab) :: alpha(mla%nlevel), beta(mla%nlevel,mla%dim)
    type(multifab) :: phi(mla%nlevel), Lphi(mla%nlevel)

    nlevs = mla%nlevel
    dm    = mla%dim

    !***********************************
    ! Allocate and build temporary arrays
    !***********************************

    do n = 1, nlevs
       call multifab_build(  phi(n),mla%la(n),    1,1)
       call multifab_build( Lphi(n),mla%la(n),    1,1)
       call multifab_build(alpha(n),mla%la(n),    1,1)
       do d = 1,dm
          call multifab_build_edge(beta(n,d),mla%la(n),1,1,d)
       end do
       call setval( phi(n),0.0_dp_t,all=.true.)
       call setval(Lphi(n),0.0_dp_t,all=.true.)
       call setval(alpha(n),ZERO, all=.true.)
       do d = 1,dm
          call setval(beta(n,d),-ONE, all=.true.)
       end do
    enddo

    !***********************************
    ! Compute explicit diffusive_term
    !***********************************
     do n = 1, nlevs
        call multifab_copy_c(phi(n),1,data(n),data_comp,1,1)
     enddo

     call cc_applyop(mla,Lphi,phi,alpha,beta,dx,the_bc_tower,&
                     bc_comp,stencil_order)

     do n = 1, nlevs
        call multifab_copy_c(lap_data(n),data_comp,Lphi(n),1)
     enddo

     do n = 1,nlevs
        call multifab_destroy(alpha(n))
        do d = 1,dm
           call multifab_destroy( beta(n,d))
        end do
        call multifab_destroy(  phi(n))
        call multifab_destroy( Lphi(n))
     enddo

  end subroutine get_explicit_diffusive_term

end module explicit_diffusive_module
