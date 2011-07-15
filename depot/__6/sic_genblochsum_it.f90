subroutine sic_genblochsum_it
use modmain
use mod_sic
use mod_nrkp
implicit none
integer ikloc,ik,j,igg,ig,ig1,ir,i,ispn,l,lm,n
integer vgl(3)
real(8) vg(3),t1,tp(2),vgk(3)
complex(8) zt1,zt2
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: wgk(:,:,:)
complex(8), allocatable :: wvgk(:,:,:)
complex(8), allocatable :: zf(:)
complex(8), allocatable :: ylmgk(:)
complex(8), allocatable :: zf1(:,:)
!
s_wkit=zzero
s_wvkit=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
allocate(wgk(ngvec,nspinor,sic_wantran%nwan))
allocate(wvgk(ngvec,nspinor,sic_wantran%nwan))
allocate(zf1(s_nr_min,lmmaxwan))
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  wgk=zzero
  wvgk=zzero
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jl,ylmgk,zf,t1,tp,ir,j,n,ispn,zt1,zt2,l,lm,ig1,vgl,igg,vgk)
  allocate(jl(s_nr_min,0:lmaxwan))
  allocate(ylmgk(lmmaxwan))
  allocate(zf(s_nr_min))
!$OMP DO
  do ig=1,ngvec
    vgk(:)=vgc(:,ig)+vkc(:,ik)
    call sphcrd(vgk,t1,tp)
    write(*,*)"ik,ig,|k+G|=",ik,ig,t1
! generate Bessel functions j_l(|G+k|x)
    do ir=1,s_nr_min
      call sbessel(lmaxwan,t1*s_r(ir),jl(ir,:))
    enddo
    call genylm(lmaxwan,tp,ylmgk)
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        do lm=1,lmmaxwan
          zf1(:,lm)=s_wlm(lm,:,ispn,j)
        enddo
! compute <G+k|W_n>
        zt1=zzero
        do l=0,lmaxwan
          zf=zzero
          do lm=l**2+1,(l+1)**2
            zf(:)=zf(:)+ylmgk(lm)*zf1(:,lm)
          enddo
          !do ir=1,s_nr_min
          !    do lm=l**2+1,(l+1)**2
          !      zf(ir)=zf(ir)+s_wlm(lm,ir,ispn,j)*ylmgk(lm)
          !    enddo
          !enddo
          zt2=zzero
          do ir=1,s_nr_min
            zt2=zt2+jl(ir,l)*zf(ir)*s_rw(ir)
          enddo
          zt1=zt1+zt2*fourpi*((-zi)**l)/sqrt(omega)
        enddo !l
        wgk(ig,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vgk(:)))
!! compute <G+k|(WV)_n>
!        zt1=zzero
!        do l=0,lmaxwan
!          zf=zzero
!          do ir=1,s_nr
!            do lm=l**2+1,(l+1)**2
!              zf(ir)=zf(ir)+s_wvlm(lm,ir,ispn,j)*ylmgk(lm)
!            enddo
!          enddo
!          zt2=zzero
!          do ir=1,s_nr
!            zt2=zt2+jl(ir,l)*zf(ir)*s_rw(ir)
!          enddo
!          zt1=zt1+zt2*fourpi*((-zi)**l)/sqrt(omega)
!        enddo !l
!        wvgk(ig,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vgkc(:,ig,1,ikloc)))
      enddo !ispn
    enddo !j
  enddo !ig
!$OMP END DO
  deallocate(jl,ylmgk,zf)
!$OMP DO 
  do ig=1,ngk(1,ik)
    do ig1=1,ngvec
      vgl(:)=ivg(:,igkig(ig,1,ikloc))-ivg(:,ig1)
      if (vgl(1).ge.intgv(1,1).and.vgl(1).le.intgv(1,2).and.&
          vgl(2).ge.intgv(2,1).and.vgl(2).le.intgv(2,2).and.&
          vgl(3).ge.intgv(3,1).and.vgl(3).le.intgv(3,2)) then
        igg=ivgig(vgl(1),vgl(2),vgl(3))
        s_wkit(ig,:,:,ikloc)=s_wkit(ig,:,:,ikloc)+&
          wgk(ig1,:,:)*cfunig(igg)
        !s_wvkit(ig,:,:,ikloc)=s_wvkit(ig,:,:,ikloc)+&
        !  wvgk(ig1,:,:)*cfunig(igg)
      endif
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo !ikloc
deallocate(wgk)
deallocate(wvgk)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_genblochsum_it] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end subroutine

