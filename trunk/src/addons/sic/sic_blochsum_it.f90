subroutine sic_blochsum_it
use modmain
use mod_sic
use mod_nrkp
implicit none
integer ikloc,ik,j,igg,ig,ig1,ir,i,ispn,l,lm,n
integer vgl(3)
real(8) vg(3),t1,tp(2)
complex(8) zt1,zt2
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: wgk(:,:,:)
complex(8), allocatable :: wvgk(:,:,:)
complex(8), allocatable :: zf(:)
complex(8), allocatable :: ylmgk(:)
!
s_wankir=zzero
s_wvkir=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
!allocate(jl(s_nr,0:lmaxwan))
!allocate(zf(s_nr))
!allocate(ylmgk(lmmaxwan))
allocate(wgk(ngkmax,nspinor,sic_wantran%nwan))
allocate(wvgk(ngkmax,nspinor,sic_wantran%nwan))
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  wgk=zzero
  wvgk=zzero
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(jl,ylmgk,zf,t1,tp,ir,j,n,ispn,zt1,zt2,l,lm,ig1,vgl,igg)
  allocate(jl(s_nr,0:lmaxwan))
  allocate(ylmgk(lmmaxwan))
  allocate(zf(s_nr))
!$OMP DO
  do ig=1,ngk(1,ik)
    call sphcrd(vgkc(1,ig,1,ikloc),t1,tp)
! generate Bessel functions j_l(|G+G'+k|x)
    do ir=1,s_nr
      call sbessel(lmaxwan,t1*s_r(ir),jl(ir,:))
    enddo
    call genylm(lmaxwan,tp,ylmgk)
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
! compute <G+k|W_n>
        zt1=zzero
        do l=0,lmaxwan
          zf=zzero
          do ir=1,s_nr
            do lm=l**2+1,(l+1)**2
              zf(ir)=zf(ir)+s_wanlm(lm,ir,ispn,j)*ylmgk(lm)
            enddo
          enddo
          zt2=zzero
          do ir=1,s_nr
            zt2=zt2+jl(ir,l)*zf(ir)*s_rw(ir)
          enddo
          zt1=zt1+zt2*fourpi*((-zi)**l)/sqrt(omega)
        enddo !l
        wgk(ig,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vgkc(:,ig,1,ikloc)))
! compute <G+k|(WV)_n>
        zt1=zzero
        do l=0,lmaxwan
          zf=zzero
          do ir=1,s_nr
            do lm=l**2+1,(l+1)**2
              zf(ir)=zf(ir)+s_wvlm(lm,ir,ispn,j)*ylmgk(lm)
            enddo
          enddo
          zt2=zzero
          do ir=1,s_nr
            zt2=zt2+jl(ir,l)*zf(ir)*s_rw(ir)
          enddo
          zt1=zt1+zt2*fourpi*((-zi)**l)/sqrt(omega)
        enddo !l
        wvgk(ig,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vgkc(:,ig,1,ikloc)))
      enddo !ispn
    enddo !j
  enddo !ig
!$OMP END DO
  deallocate(jl,ylmgk,zf)
!$OMP DO 
  do ig=1,ngk(1,ik)
    do ig1=1,ngk(1,ik)
      vgl(:)=ivg(:,igkig(ig,1,ikloc))-ivg(:,igkig(ig1,1,ikloc))
      igg=ivgig(vgl(1),vgl(2),vgl(3))
      s_wankir(ig,:,:,ikloc)=s_wankir(ig,:,:,ikloc)+&
        wgk(ig1,:,:)*cfunig(igg)
      s_wvkir(ig,:,:,ikloc)=s_wvkir(ig,:,:,ikloc)+&
        wvgk(ig1,:,:)*cfunig(igg)
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
enddo !ikloc
deallocate(wgk)
deallocate(wvgk)
!deallocate(jl)
!deallocate(zf)
!deallocate(ylmgk)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum_it] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end

