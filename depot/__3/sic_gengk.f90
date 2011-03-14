subroutine sic_gengk
use modmain
use mod_sic
use mod_nrkp
implicit none
integer ikloc,ik,j,igg,ig,ig1,ir,i,ispn,l,lm
complex(8), allocatable :: wgk(:,:,:)
complex(8), allocatable :: wggk(:,:,:)
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: zf(:)
complex(8) zt1,zt2
integer vgl(3),glim(2,3),gglim(2,3)
real(8) vg(3),t1,tp(2)
complex(8), allocatable :: ylmgk(:)
complex(8), allocatable :: zfir(:)
integer, allocatable :: ggidx(:,:,:)
integer, allocatable :: ggdone(:)
integer ngg
write(*,*)"in sic_gengk"

allocate(wgk(ngkmax,nspinor,sic_wantran%nwan))
allocate(jl(s_nr,0:lmaxwan))
allocate(zf(s_nr))
allocate(ylmgk(lmmaxwan))
allocate(zfir(ngrtot))

do ikloc=1,1 !nkptloc
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
  write(*,*)"ngg:",ngg
  allocate(wggk(ngg,nspinor,sic_wantran%nwan))
  allocate(ggdone(ngg))
  wggk=zzero
  ggdone=0
  do ig=1,ngk(1,ik)
    write(*,*)"ig=",ig
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
          do ispn=1,nspinor
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
            wggk(igg,ispn,j)=zt1
          enddo !ispn
        enddo !j
      endif
    enddo !ig1
  enddo !ig
  wgk=zzero
  do ig=1,ngk(1,ik)
    do ig1=1,ngrtot
      vgl(:)=ivg(:,igkig(ig,1,ikloc))+ivg(:,ig1)
      igg=ggidx(vgl(1),vgl(2),vgl(3))
      if (igg.ne.0) wgk(ig,:,:)=wgk(ig,:,:)+wggk(igg,:,:)*cfunig(ig1)
    enddo
  enddo
  zt1=zzero
  do ig=1,ngk(1,ik)
    zt1=zt1+dconjg(wgk(ig,1,1))*wann_unkit(ig,1,1,1)
  enddo
  write(*,*)zt1
  deallocate(ggidx)
  deallocate(ggdone)
  deallocate(wggk)
enddo !ikloc

open(220,file="cfunir.dat",form="formatted",status="replace")
do ir=1,ngrtot
  if (abs(vgrc(1,ir)-vgrc(2,ir)).lt.1d-10.and.&
      abs(vgrc(2,ir)-vgrc(3,ir)).lt.1d-10) then
    write(220,'(2G18.10)')sqrt(sum(vgrc(:,ir)**2)),cfunir(ir)
  endif
enddo
close(220)

zfir=zzero
do ig=1,ngknr(1)
  zfir(igfft(igkignr(ig,1)))=wann_unkit(ig,1,1,1)/sqrt(omega)
enddo
call zfftifc(3,ngrid,1,zfir)
open(220,file="wan_unk.dat",form="formatted",status="replace")
do ir=1,ngrtot
  if (abs(vgrc(1,ir)-vgrc(2,ir)).lt.1d-10.and.&
      abs(vgrc(2,ir)-vgrc(3,ir)).lt.1d-10) then
    write(220,'(3G18.10)')sqrt(sum(vgrc(:,ir)**2)),dreal(zfir(ir)),dimag(zfir(ir))
  endif
enddo
close(220)

open(220,file="wan_unk_dot_cfunir.dat",form="formatted",status="replace")
do ir=1,ngrtot
  if (abs(vgrc(1,ir)-vgrc(2,ir)).lt.1d-10.and.&
      abs(vgrc(2,ir)-vgrc(3,ir)).lt.1d-10) then
    write(220,'(3G18.10)')sqrt(sum(vgrc(:,ir)**2)),dreal(cfunir(ir)*zfir(ir)),&
      dimag(cfunir(ir)*zfir(ir))
  endif
enddo
close(220)

zfir=zzero
do ig=1,ngknr(1)
  zfir(igfft(igkignr(ig,1)))=wgk(ig,1,1)/sqrt(omega)
enddo
call zfftifc(3,ngrid,1,zfir)
open(220,file="wan_unk_dot_cfunig.dat",form="formatted",status="replace")
do ir=1,ngrtot
  if (abs(vgrc(1,ir)-vgrc(2,ir)).lt.1d-10.and.&
      abs(vgrc(2,ir)-vgrc(3,ir)).lt.1d-10) then
    write(220,'(3G18.10)')sqrt(sum(vgrc(:,ir)**2)),dreal(zfir(ir)),dimag(zfir(ir))
  endif
enddo
close(220)



deallocate(wgk)
deallocate(jl)
deallocate(zf)
deallocate(ylmgk)

return
end

