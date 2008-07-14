subroutine zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_resp,zrhofc)
use modmain
implicit none
integer    ,intent(in   ) :: ngvec_resp
complex(8) ,intent(in)    :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
complex(8) ,intent(inout) :: zrhoir(ngrtot)
real(8)    ,intent(in) :: jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec)
complex(8) ,intent(in) :: ylmgq0(lmmaxvr,ngvec)
complex(8) ,intent(in) :: sfacgq0(ngvec,natmtot)
complex(8) ,intent(out) :: zrhofc(ngvec_resp)

integer                :: ig,is,ia,ias,nr,ir,l,m,lm
real(8)                :: t1,t2
complex(8)             :: zsum1,zsum2
real(8) ,allocatable   :: fr1(:),fr2(:),gr(:),cf(:,:)
  
allocate(fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax))
  
zrhofc=dcmplx(0.d0,0.d0)

! muffin-tin part
do ig=1,ngvec_resp
  do is=1,nspecies
    nr=nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nr
        zsum1=dcmplx(0.d0,0.d0)  
        do l=0,lmaxvr
          zsum2=dcmplx(0.d0,0.d0)
          do m=-l,l
            lm=idxlm(l,m)
            zsum2=zsum2+zrhomt(lm,ir,ias)*ylmgq0(lm,ig)
          enddo !m
!         i^l*j_l(|G+q0|x)*\sum_m...
          zsum1=zsum1+jlgq0r(ir,l,is,ig)*dconjg(zil(l))*zsum2
        enddo !l
        t1=rcmt(ir,is)**2
        fr1(ir)=dreal(zsum1)*t1
        fr2(ir)=dimag(zsum1)*t1
      enddo !ir
      call fderiv(-1,nr,rcmt(:,is),fr1,gr,cf)
      t1=gr(nr)
      call fderiv(-1,nr,rcmt(:,is),fr2,gr,cf)
      t2=gr(nr)
      zrhofc(ig)=zrhofc(ig)+(fourpi)*dconjg(sfacgq0(ig,ias))*dcmplx(t1,t2)
    enddo !ia
  enddo !is
enddo !ig
! interstitial part
do ir=1,ngrtot
  zrhoir(ir)=zrhoir(ir)*cfunir(ir)*omega
enddo
call zfftifc(3,ngrid,-1,zrhoir)
do ig=1,ngvec_resp
  zrhofc(ig)=zrhofc(ig)+zrhoir(igfft(ig))
enddo

!do ig=1,ngvec_resp
!  if (abs(zrhofc(ig)).lt.1d-3) zrhofc(ig)=dcmplx(0.d0,0.d0)
!enddo

deallocate(fr1,fr2,gr,cf)
return
end
