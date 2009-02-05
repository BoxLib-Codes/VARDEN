!
! Computes error for convergence study
!

module convergence

  use bl_types
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use average
  use fabio_module

  implicit none
 
contains
 
  subroutine conv(loc,files)
    integer, intent(in   ) :: loc(:)
    integer, intent(in   ) :: files(:)

    !local variables   
    integer, parameter :: dm        = 2, &
                          levs      = 1, &
                          nscal     = 4, &
                          nspec     = 3, &
                          n_grids   = 4, &
                          ng_cell   = 3, &
                          ng_grow   = 1, &
                          exact_loc = 2048

    integer                   :: nlevs,n,i
    integer                   :: ir(n_grids,dm)
    type(ml_layout)           :: mla(n_grids+1)
    real(dp_t)                :: time,dt
    type(multifab)            :: s_exact,s(n_grids)
    type(multifab)            :: uold,gp,p
    type(ml_boxarray)         :: mba
    logical                   :: nodal(dm)
    logical                   :: pmask(dm)
    type(multifab)            :: avg(n_grids)
    real(dp_t)                :: l2_error(n_grids,nspec+1), l1_error(n_grids,nspec+1)
    real(dp_t)                :: dx(n_grids,dm)
    character(len=20)         :: plot_names(nscal)
    character(len=20)         :: sd_name
    type(multifab), allocatable  :: plotdata(:)


    nlevs = 1
    nodal(:) = .true.
    pmask(:) = .true.
    
    ir(1,:) = 64
    dx(1,:) = 1.d0/dble(ir(1,:))

    do n = 2, n_grids
       ir(n,:) = ir(n-1,:)/2.d0
       dx(n,:) = dx(n-1,:)/2.d0
    end do

    plot_names(1) = 'density'
    plot_names(2) = 'species A'
    plot_names(3) = 'species B'
    plot_names(4) = 'species C'



    write(*,*) loc(1),loc(2),loc(3),loc(4)
    write(*,*) files(1),files(2)
! read  in data
    call initialize_from_restart(mla(n_grids+1),exact_loc,time,dt,pmask,uold,&
                                 s_exact,gp,p)
    call destroy(uold)
    call destroy(gp)
    call destroy(p)

    do n = 1, n_grids
       call initialize_from_restart(mla(n),loc(n),time,dt,pmask,uold,&
                                    s(n),gp,p)
       call destroy(uold)
       call destroy(gp)
       call destroy(p)
    enddo

!call print(s(1))
!call print(s(2))
!call print(s(3))
!call print(s(4))

! avg the exact data down onto course grids
    do n = 1,n_grids
       call multifab_build(avg(n),s(n)%la,nscal)
       call ml_cc_restriction(avg(n),s_exact,ir(n,:))
       call write_plotfile(n, mla(n),avg(n))
       call multifab_sub_sub(avg(n),s(n))
       call write_plotfile(n+nscal, mla(n),avg(n))
    enddo

! calculate l1 error
    do i = 1,n_grids
       l1_error(i,nspec+1) = 0.d0
       do n = 1,nspec
          l1_error(i,n) = multifab_norm_l1_c(avg(i),n+1,1,all=.false.)
          l1_error(i,nspec+1) = l1_error(i,nspec+1) + l1_error(i,n)
       enddo
       l1_error(i,:) = l1_error(i,:)*(dx(i,1)*dx(i,2))
       l1_error(i,nspec+1) = l1_error(i,nspec+1)/dble(nspec)
    end do


! calculate l2 error
    do i = 1,n_grids
       l2_error(i,nspec+1) = 0.d0
       do n = 1,nspec
          l2_error(i,n) = multifab_norm_l2_c(avg(i),n+1,1,all=.false.)
          l2_error(i,nspec+1) = l2_error(i,nspec+1) + l2_error(i,n)
       enddo
       l2_error(i,:) = l2_error(i,:)*sqrt(dx(i,1)*dx(i,2))
       l2_error(i,nspec+1) = l2_error(i,nspec+1)/dble(nspec)
    end do

    write(files(1),*)'#Convergence study data -- L1 norm'
    write(files(1),*)'#    dx     SpeciesA     SpeciesB    SpeciesC     Combined'
    do n = 1, n_grids
       write(files(1),1000) dx(n,1),l1_error(n,1),l1_error(n,2),l1_error(n,3),&
                   l1_error(n,4)
    end do

    write(files(2),*)'#Convergence study data -- L2 norm'
    write(files(2),*)'#    dx     SpeciesA     SpeciesB    SpeciesC    Combined'
    do n = 1, n_grids
       write(files(2),1000) dx(n,1),l2_error(n,1),l2_error(n,2),l2_error(n,3),&
                   l2_error(n,4)
    end do

    do n = 1,n_grids+1
       call destroy(mla(n))
    end do
    
1000 FORMAT(11(E15.8,1X)) 

contains

   subroutine initialize_from_restart(mla,restart,time,dt,pmask,uold,sold,gp,p)
 
     type(ml_layout),intent(out)   :: mla
     integer       , intent(in   ) :: restart
     real(dp_t)    , intent(  out) :: time,dt
     logical       , intent(in   ) :: pmask(:)
     type(multifab)       :: sold,uold,gp,p

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(multifab), pointer   :: chk_p(:)

     integer :: n


     call fill_restart_data(restart,mba,chkdata,chk_p,time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     call multifab_build(   uold, mla%la(1),    dm, ng_cell)
     call multifab_build(   sold, mla%la(1), nscal, ng_cell)
     call multifab_build(     gp, mla%la(1),    dm, ng_grow)
     call multifab_build(      p, mla%la(1),     1, ng_grow, nodal)

     call multifab_copy_c(uold,1,chkdata(1),1         ,dm)
     call multifab_copy_c(sold,1,chkdata(1),1+dm      ,nscal)
     call multifab_copy_c(  gp,1,chkdata(1),1+dm+nscal,dm)
     call multifab_copy_c(   p,1,  chk_p(1),1         ,1)
     !
     ! The layouts for chkdata and chk_p are built standalone, level
     ! by level, and need to be destroy()d as such as well.
     !
     do n = 1,nlevs     
        call destroy(chkdata(n)%la)
        call destroy(chk_p(n)%la)
        call multifab_destroy(chkdata(n))
        call multifab_destroy(chk_p(n))
     end do
     deallocate(chkdata,chk_p)
     call destroy(mba)

  end subroutine initialize_from_restart

  subroutine fill_restart_data(restart_int,mba,chkdata,chk_p,time,dt)

    integer          , intent(in   ) :: restart_int
    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba
    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chk_p(:)

    character(len=7)                  :: sd_name
    integer                           :: n
    integer                           :: rrs(nlevs)

    write(unit=sd_name,fmt='("chk",i4.4)') restart_int
    print *,'Reading ',sd_name,' to get state data for restart'
    call checkpoint_read(chkdata, chk_p, sd_name, rrs, time, dt, nlevs)

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = rrs(n-1)
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

  end subroutine fill_restart_data

  subroutine write_plotfile(istep_to_write, mla, mf)

    integer,         intent(in   ) :: istep_to_write
    type(ml_layout), intent(in   ) :: mla
    type(multifab),  intent(in   ) :: mf
  
    integer                        :: n,n_plot_comps

    allocate(plotdata(levs))
    n_plot_comps = nscal

    do n = 1,levs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy_c(plotdata(n),1        ,mf,1,nscal)

    end do

    write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), time, dx(1,:))

    do n = 1,levs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)

  end subroutine write_plotfile

end subroutine conv

end module convergence
