subroutine sic_blochsum_it
use modmain
use mod_sic
use mod_nrkp
implicit none
integer ikloc,ik,j,igg,ig,ig1,ir,i,ispn,l,lm,n
integer vgl(3),glim(2,3),gglim(2,3),ngg
real(8) vg(3),t1,tp(2)
complex(8) zt1,zt2
integer, allocatable :: ggdone(:)
integer, allocatable :: ggidx(:,:,:)
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: wggk(:,:,:)
complex(8), allocatable :: wvggk(:,:,:)
complex(8), allocatable :: zf(:)
complex(8), allocatable :: ylmgk(:)
!
s_wankir=zzero
s_wvkir=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
allocate(jl(s_nr,0:lmaxwan))
allocate(zf(s_nr))
allocate(ylmgk(lmmaxwan))
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  glim=0
  do ig=1,ngk(1,ik)
    do i=1,3
      glim(1,i)=min(glim(1,i),ivg(i,igkig(ig,1,ikloc)))
      glim(2,i)=max(glim(2,i),ivg(i,igkig(ig,1,ikloc)))
    enddo
  enddo
  do i=1,3
    gglim(1,i)=glim(1,i)+intgv(i,1)
    gglim(2,i)=glim(2,i)+intgv(i,2)
  enddo
  allocate(ggidx(gglim(1,1):gglim(2,1),gglim(1,2):gglim(2,2),gglim(1,3):gglim(2,3)))
  ngg=0
  ggidx=0
  do ig=1,ngk(1,ik)
    do ig1=1,ngrtot
      vgl(:)=ivg(:,igkig(ig,1,ikloc))+ivg(:,ig1)
      if (ggidx(vgl(1),vgl(2),vgl(3)).eq.0) then
        ngg=ngg+1
        ggidx(vgl(1),vgl(2),vgl(3))=ngg
      endif
    enddo
  enddo
  allocate(wggk(ngg,nspinor,sic_wantran%nwan))
  allocate(wvggk(ngg,nspinor,sic_wantran%nwan))
  allocate(ggdone(ngg))
  wggk=zzero
  wvggk=zzero
  ggdone=0
  do ig=1,ngk(1,ik)
    do ig1=1,ngrtot
      vgl(:)=ivg(:,igkig(ig,1,ikloc))+ivg(:,ig1)
      igg=ggidx(vgl(1),vgl(2),vgl(3))
      if (igg.ne.0.and.ggdone(igg).eq.0) then
        ggdone(igg)=1
        vg(:)=vgkc(:,ig,1,ikloc)+vgc(:,ig1)
        call sphcrd(vg,t1,tp)
! generate Bessel functions j_l(|G+G'+k|x)
        do ir=1,s_nr
          call sbessel(lmaxwan,t1*s_r(ir),jl(ir,:))
        enddo
        call genylm(lmaxwan,tp,ylmgk)
        do j=1,sic_wantran%nwan
          n=sic_wantran%iwan(j)
          do ispn=1,nspinor
! compute <G+G'+k|W_n>
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
            wggk(igg,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vg(:)))
! compute <G+G'+k|(WV)_n>
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
            wvggk(igg,ispn,j)=zt1*exp(-zi*dot_product(wanpos(:,n),vg(:)))
           enddo !ispn
        enddo !j
      endif
    enddo !ig1
  enddo !ig
  do ig=1,ngk(1,ik)
    do ig1=1,ngrtot
      vgl(:)=ivg(:,igkig(ig,1,ikloc))+ivg(:,ig1)
      igg=ggidx(vgl(1),vgl(2),vgl(3))
      if (igg.ne.0) then
        s_wankir(ig,:,:,ikloc)=s_wankir(ig,:,:,ikloc)+&
          wggk(igg,:,:)*cfunig(ig1)
        s_wvkir(ig,:,:,ikloc)=s_wvkir(ig,:,:,ikloc)+&
          wvggk(igg,:,:)*cfunig(ig1)
      endif
    enddo
  enddo
  deallocate(ggidx)
  deallocate(ggdone)
  deallocate(wggk)
  deallocate(wvggk)
enddo !ikloc
deallocate(jl)
deallocate(zf)
deallocate(ylmgk)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum_it] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end

