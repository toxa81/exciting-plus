subroutine genixc(ixcft)
use modmain
implicit none
complex(8), intent(out) :: ixcft(ngvec)

integer is,ia,ias,l,m,lm,ir,ig,nr,itp
real(8) t1,t2,rt1
complex(8) zsum1,zsum2
real(8), allocatable :: rftp1(:)
real(8), allocatable :: rftp2(:)
real(8), allocatable :: jl(:,:)
real(8), allocatable :: fr1(:)
real(8), allocatable :: fr2(:)
real(8), allocatable :: gr(:)
real(8), allocatable :: cf(:,:)
complex(8), allocatable :: zt1(:)
complex(8), allocatable :: zt2(:,:,:)
complex(8), allocatable :: zt3(:)

allocate(rftp1(lmmaxvr))
allocate(rftp2(lmmaxvr))
allocate(jl(0:lmaxvr,nrmtmax))
allocate(fr1(nrmtmax))
allocate(fr2(nrmtmax))
allocate(gr(nrmtmax))
allocate(cf(3,nrmtmax))
allocate(zt1(lmmaxvr))
allocate(zt2(lmmaxvr,nrmtmax,natmtot))
allocate(zt3(ngrtot))

ixcft=dcmplx(0.d0,0.d0)

! calculate Ixc(r)=Bxc(r)/m(r) inside muffin-tins
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
! transform Bxc and m from real spherical harmonics to spherical coordinates 
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,bxcmt(1,ir,ias,1),1, &
        0.d0,rftp1,1)
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,magmt(1,ir,ias,1),1, &
        0.d0,rftp2,1)
! calculate I(r)
      do itp=1,lmmaxvr
        if (abs(rftp2(itp)).lt.1d-10) rftp2(itp)=1d10
        zt1(itp)=dcmplx(rftp1(itp)/rftp2(itp),0.d0)
      enddo
! transform I(r) from spherical coordinates to complex spherical harmonics
      call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zt1,1,zzero, &
        zt2(1,ir,ias),1)
    enddo
  enddo
enddo
! calculate muffin-tin part of Fourier transform of Ixc(r) 
do ig=1,ngvec  
  do is=1,nspecies
    nr=nrmt(is)
! generate Bessel functions
    do ir=1,nr
      rt1=gc(ig)*spr(ir,is)
      call sbessel(lmaxvr,rt1,jl(0,ir))
    enddo
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nr
        zsum1=dcmplx(0.d0,0.d0)
        do l=0,lmaxvr
          zsum2=dcmplx(0.d0,0.d0)
          do m=-l,l
            lm=idxlm(l,m)
            zsum2=zsum2+zt2(lm,ir,ias)*ylmg(lm,ig)
          enddo !m
          zsum1=zsum1+jl(l,ir)*dconjg(zil(l))*zsum2
        enddo !l
        rt1=spr(ir,is)**2
        fr1(ir)=dreal(zsum1)*rt1
        fr2(ir)=dimag(zsum1)*rt1
      enddo !ir
      call fderiv(-1,nr,spr(1,is),fr1,gr,cf)
      t1=gr(nr)
      call fderiv(-1,nr,spr(1,is),fr2,gr,cf)
      t2=gr(nr)
      ixcft(ig)=ixcft(ig)+fourpi*dconjg(sfacg(ig,ias))*dcmplx(t1,t2)
    enddo !ia
  enddo !is
enddo !ig

! calculate Ixc(r)=Bxc(r)/m(r) in interstitial
do ir=1,ngrtot
  rt1=magir(ir,1)
  if (abs(rt1).lt.1d-10) rt1=1d10
  zt3(ir)=dcmplx(bxcir(ir,1)/rt1,0.d0)*cfunir(ir)*omega
enddo       
call zfftifc(3,ngrid,-1,zt3)
do ig=1,ngvec
  ixcft(ig)=ixcft(ig)+zt3(igfft(ig))
enddo 
ixcft=ixcft/omega         

deallocate(rftp1)
deallocate(rftp2)
deallocate(jl)
deallocate(fr1)
deallocate(fr2)
deallocate(gr)
deallocate(cf)
deallocate(zt1)
deallocate(zt2)
deallocate(zt3)
return
end