module update_module
  
  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: update

contains

  subroutine update(mla,sold,umac,sedge,flux,force,snew,dx,dt,is_vel,is_cons, &
                    the_bc_level)

    use bl_constants_module
    use ml_restrict_fill_module

    type(ml_layout)   , intent(in   ) :: mla
    type(multifab)    , intent(in   ) :: sold(:)
    type(multifab)    , intent(in   ) :: umac(:,:)
    type(multifab)    , intent(in   ) :: sedge(:,:)
    type(multifab)    , intent(in   ) :: flux(:,:)
    type(multifab)    , intent(in   ) :: force(:)
    type(multifab)    , intent(inout) :: snew(:)
    real(kind = dp_t) , intent(in   ) :: dx(:,:),dt
    logical           , intent(in   ) :: is_vel,is_cons(:)
    type(bc_level)    , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpx(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpy(:,:,:,:)
    real(kind=dp_t), pointer :: fluxpz(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    integer :: lo(get_dim(sold(1))),hi(get_dim(sold(1)))
    integer :: i,dm,nscal,n,nlevs
    integer :: ng_s,ng_u,ng_e,ng_f,ng_o
    
    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"update")

    nlevs = mla%nlevel
    dm    = mla%dim

    nscal = multifab_ncomp(sold(1))

    ng_s = sold(1)%ng
    ng_u = umac(1,1)%ng
    ng_e = sedge(1,1)%ng
    ng_f = flux(1,1)%ng
    ng_o = force(1)%ng

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(sop,snp,ump,vmp,wmp,sepx,sepy,sepz) &
    !$omp private(fluxpx,fluxpy,fluxpz,fp,lo,hi)

    do n=1,nlevs
       call mfiter_build(mfi, sold(n), tiling=.true.)

          do while (more_tile(mfi))
             i = get_fab_index(mfi)

             tilebox = get_tilebox(mfi)
             tlo = lwb(tilebox)
             thi = upb(tilebox)

!       do i = 1, nfabs(sold(n))
          sop    => dataptr(sold(n),i)
          snp    => dataptr(snew(n),i)
          ump    => dataptr(umac(n,1),i)
          vmp    => dataptr(umac(n,2),i)
          sepx   => dataptr(sedge(n,1),i)
          sepy   => dataptr(sedge(n,2),i)
          fluxpx => dataptr(flux(n,1),i)
          fluxpy => dataptr(flux(n,2),i)
          fp     => dataptr(force(n),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call update_2d(sop(:,:,1,:), ump(:,:,1,1), vmp(:,:,1,1), &
                            sepx(:,:,1,:), sepy(:,:,1,:), &
                            fluxpx(:,:,1,:), fluxpy(:,:,1,:), &
                            fp(:,:,1,:) , snp(:,:,1,:), &
                            lo, hi, ng_s, ng_u, ng_e, ng_f, ng_o, &
                            dx(n,:), dt, is_vel, is_cons,tlo,thi)
          case (3)
             wmp    => dataptr( umac(n,3),i)
             sepz   => dataptr(sedge(n,3),i)
             fluxpz => dataptr( flux(n,3),i)
             call update_3d(sop(:,:,:,:), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                            sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                            fluxpx(:,:,:,:), fluxpy(:,:,:,:), fluxpz(:,:,:,:), &
                            fp(:,:,:,:) , snp(:,:,:,:), &
                            lo, hi, ng_s, ng_u, ng_e, ng_f, ng_o, &
                            dx(n,:), dt, is_vel, is_cons,tlo,thi)
          end select
       end do

    enddo ! end loop over levels
    !$omp end parallel

    if (is_vel) then
       call ml_restrict_and_fill(nlevs, snew, mla%mba%rr, the_bc_level, bcomp=1)
    else
       call ml_restrict_and_fill(nlevs, snew, mla%mba%rr, the_bc_level, bcomp=dm+1)
    end if

    call destroy(bpt)

  end subroutine update

  subroutine update_2d(sold,umac,vmac,sedgex,sedgey,fluxx,fluxy,force,snew,&
                       glo,ghi,ng_s,ng_u, ng_e, ng_f, ng_o, dx,dt,is_vel,is_cons,tlo,thi)

    use bl_constants_module

    integer           , intent(in   ) :: glo(:), ghi(:), ng_s, ng_u, ng_e, ng_f, ng_o
    integer           , intent(in   ) :: tlo(:), thi(:)
    real (kind = dp_t), intent(in   ) ::    sold(glo(1)-ng_s:,glo(2)-ng_s:,:)  
    real (kind = dp_t), intent(  out) ::    snew(glo(1)-ng_s:,glo(2)-ng_s:,:)  
    real (kind = dp_t), intent(in   ) ::    umac(glo(1)-ng_u:,glo(2)-ng_u:)  
    real (kind = dp_t), intent(in   ) ::    vmac(glo(1)-ng_u:,glo(2)-ng_u:)  
    real (kind = dp_t), intent(in   ) ::  sedgex(glo(1)-ng_e:,glo(2)-ng_e:,:)  
    real (kind = dp_t), intent(in   ) ::  sedgey(glo(1)-ng_e:,glo(2)-ng_e:,:)  
    real (kind = dp_t), intent(in   ) ::   fluxx(glo(1)-ng_f:,glo(2)-ng_f:,:)  
    real (kind = dp_t), intent(in   ) ::   fluxy(glo(1)-ng_f:,glo(2)-ng_f:,:) 
    real (kind = dp_t), intent(in   ) ::   force(glo(1)-ng_o:,glo(2)-ng_o:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: is_vel
    logical           , intent(in   ) :: is_cons(:)

    integer :: i, j, comp
    real (kind = dp_t) ubar,vbar
    real (kind = dp_t) ugradu,ugradv,ugrads
    real (kind = dp_t) :: divsu

    if (is_vel) then

       do j = tlo(2), thi(2)
          do i = tlo(1), thi(1)

             ubar = HALF*(umac(i,j) + umac(i+1,j))
             vbar = HALF*(vmac(i,j) + vmac(i,j+1))

             ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
                      vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

             ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
                      vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

             snew(i,j,1) = sold(i,j,1) - dt * ugradu + dt * force(i,j,1)
             snew(i,j,2) = sold(i,j,2) - dt * ugradv + dt * force(i,j,2)

          enddo
       enddo

    else

       do comp = 1,size(sold,dim=3)
          if (is_cons(comp)) then
             do j = tlo(2), thi(2)
                do i = tlo(1), thi(1)
                   divsu = (fluxx(i+1,j,comp)-fluxx(i,j,comp))/dx(1) &
                         + (fluxy(i,j+1,comp)-fluxy(i,j,comp))/dx(2)
                   snew(i,j,comp) = sold(i,j,comp) - dt * divsu + dt * force(i,j,comp)
                enddo
             enddo
          else
             do j = tlo(2), thi(2)
                do i = tlo(1), thi(1)
                   ubar = HALF*(umac(i,j) + umac(i+1,j))
                   vbar = HALF*(vmac(i,j) + vmac(i,j+1))
                   ugrads = ubar*(sedgex(i+1,j,comp) - sedgex(i,j,comp))/dx(1) + &
                            vbar*(sedgey(i,j+1,comp) - sedgey(i,j,comp))/dx(2)
                   snew(i,j,comp) = sold(i,j,comp) - dt * ugrads + dt * force(i,j,comp)
                enddo
             enddo
          end if
       end do

    end if

  end subroutine update_2d

  subroutine update_3d(sold,umac,vmac,wmac,sedgex,sedgey,sedgez,fluxx,fluxy,fluxz, &
                       force,snew,glo,ghi,ng_s,ng_u,ng_e,ng_f,ng_o,dx,dt,is_vel,is_cons,tlo,thi)

    use bl_constants_module

    integer           , intent(in   ) :: glo(:), ghi(:), ng_s, ng_u, ng_e, ng_f, ng_o
    integer           , intent(in   ) :: tlo(:), thi(:)
    real (kind = dp_t), intent(in   ) ::    sold(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)  
    real (kind = dp_t), intent(  out) ::    snew(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)  
    real (kind = dp_t), intent(in   ) ::    umac(glo(1)-ng_u:,glo(2)-ng_u:,glo(3)-ng_u:)  
    real (kind = dp_t), intent(in   ) ::    vmac(glo(1)-ng_u:,glo(2)-ng_u:,glo(3)-ng_u:)  
    real (kind = dp_t), intent(in   ) ::    wmac(glo(1)-ng_u:,glo(2)-ng_u:,glo(3)-ng_u:)  
    real (kind = dp_t), intent(in   ) ::  sedgex(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)  
    real (kind = dp_t), intent(in   ) ::  sedgey(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)  
    real (kind = dp_t), intent(in   ) ::  sedgez(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:,:)  
    real (kind = dp_t), intent(in   ) ::   fluxx(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)  
    real (kind = dp_t), intent(in   ) ::   fluxy(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)  
    real (kind = dp_t), intent(in   ) ::   fluxz(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:) 
    real (kind = dp_t), intent(in   ) ::   force(glo(1)-ng_o:,glo(2)-ng_o:,glo(3)-ng_o:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: dt
    logical           , intent(in   ) :: is_vel
    logical           , intent(in   ) :: is_cons(:)

    !     Local variables
    integer :: i, j, k, comp
    real (kind = dp_t) ubar,vbar,wbar
    real (kind = dp_t) :: ugradu,ugradv,ugradw,ugrads
    real (kind = dp_t) :: divsu

    if (is_vel) then

       do k = tlo(3), thi(3)
          do j = tlo(2), thi(2)
             do i = tlo(1), thi(1)
                ubar = half*(umac(i,j,k) + umac(i+1,j,k))
                vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
                wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))

                ugradu = ubar*(sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1) + &
                         vbar*(sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2) + &
                         wbar*(sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)

                ugradv = ubar*(sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1) + &
                         vbar*(sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2) + &
                         wbar*(sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)

                ugradw = ubar*(sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1) + &
                         vbar*(sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2) + &
                         wbar*(sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)

                snew(i,j,k,1) = sold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
                snew(i,j,k,2) = sold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
                snew(i,j,k,3) = sold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)
             enddo
          enddo
       enddo

    else

       do comp = 1,size(sold,dim=4)
          if (is_cons(comp)) then
             do k = tlo(3), thi(3)
                do j = tlo(2), thi(2)
                   do i = tlo(1), thi(1)
                      divsu = (fluxx(i+1,j,k,comp)-fluxx(i,j,k,comp))/dx(1) &
                            + (fluxy(i,j+1,k,comp)-fluxy(i,j,k,comp))/dx(2) &
                            + (fluxz(i,j,k+1,comp)-fluxz(i,j,k,comp))/dx(3)
                      snew(i,j,k,comp) = sold(i,j,k,comp) - dt * divsu + dt * force(i,j,k,comp)
                   enddo
                enddo
             enddo

          else 

             do k = tlo(3), thi(3)
                do j = tlo(2), thi(2)
                   do i = tlo(1), thi(1)
                      ubar = half*(umac(i,j,k) + umac(i+1,j,k))
                      vbar = half*(vmac(i,j,k) + vmac(i,j+1,k))
                      wbar = half*(wmac(i,j,k) + wmac(i,j,k+1))
                      ugrads = ubar*(sedgex(i+1,j,k,comp) - sedgex(i,j,k,comp))/dx(1) + &
                               vbar*(sedgey(i,j+1,k,comp) - sedgey(i,j,k,comp))/dx(2) + &
                               wbar*(sedgez(i,j,k+1,comp) - sedgez(i,j,k,comp))/dx(3)
                      snew(i,j,k,comp) = sold(i,j,k,comp) - dt * ugrads + dt * force(i,j,k,comp)
                   enddo
                enddo
             enddo
          end if
       enddo

    end if

  end subroutine update_3d

end module update_module
