module   response
implicit none

contains

  !calculates integral <psi_{nk}|e^{i(G+q)r|psi_{n'k'}> = int_{BZ} rho_{nk,n'k'}(r) e^{i(G+q)r dr
  !where k' = k+q 
  !subroutine zrhogqint(zrhomt,zrhoir)
  subroutine zrhogqint
  use        modmain
  implicit   none
  !complex(8) ,intent(in) :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
  !complex(8) ,intent(in) :: zrhoir(ngrtot)
  complex(8) ,allocatable :: zrhomt(:,:,:)
  complex(8) ,allocatable :: zrhoir(:)
  
  real(8)                :: vq0l(3),vq0c(3)
  integer                :: ngq0
  !integer ,allocatable   :: igq0ig(:)
  !real(8) ,allocatable   :: vgq0l(:,:)
  real(8) ,allocatable   :: vgq0c(:,:)
  real(8) ,allocatable   :: gq0(:)
  real(8) ,allocatable   :: tpgq0(:,:)
  complex(8) ,allocatable   :: sfacgq0(:,:)
  integer                :: ig,is,ia,ias,nr,ir,l,m,lm
  complex(8) ,allocatable :: ylmgq0(:,:)
  real(8) ,allocatable   :: jlgq0r(:,:)
  real(8)                :: t1,t2
  complex(8)             :: zsum1,zsum2,zsum
  real(8) ,allocatable   :: fr1(:),fr2(:),gr(:),cf(:,:)
  
  call init0
  call init1
  
  !allocate(igq0ig(ngkmax))
  !allocate(vgq0l(3,ngkmax))
  allocate(vgq0c(3,ngvec))
  allocate(gq0(ngvec))
  allocate(tpgq0(2,ngvec))
  allocate(sfacgq0(ngvec,natmtot))
  allocate(ylmgq0(lmmaxvr,ngvec))
  allocate(jlgq0r(0:lmaxvr,nrcmtmax))
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
  allocate(zrhoir(ngrtot))
  allocate(fr1(nrcmtmax),fr2(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax))
  
  vq0l = (/0.5d0,0.5d0,0.5d0/)
  
  !--- get q0 in lattice coordinates
  call r3mv(bvec,vq0l,vq0c)
  
  !--- generate G+q0 vectors
  !call gengpvec(vq0l,vq0c,ngq0,igq0ig,vgq0l,vgq0c,gq0,tpgq0)
  do ig = 1, ngvec
    vgq0c(:,ig) = vgc(:,ig) + vq0c(:)
    !--- get spherical coordinates and length of G+q0
    call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
    !--- generate spherical harmonics for G+q0
    call genylm(lmaxvr,tpgq0(:,ig),ylmgq0(:,ig))
  enddo
  
  !--- generate structure factor for G+q0 vectors
  call gensfacgp(ngvec,vgq0c,ngvec,sfacgq0)
  
  do ig = 1, 20
  
  zsum = dcmplx(0.d0,0.d0)
  do is = 1, nspecies
    nr = nrcmt(is)
    do ir = 1, nr
      !--- |G+q0|*x
      t1 = gq0(ig)*rcmt(ir,is)
      call sbessel(lmaxvr,t1,jlgq0r(:,ir))
    enddo
    do ia = 1, natoms(is)
      ias = idxas(ia,is)
      do ir = 1, nr
        zsum1 = dcmplx(0.d0,0.d0)  
        do l = 0, lmaxvr
          zsum2 = dcmplx(0.d0,0.d0)
          do m = -l,l
            lm = idxlm(l,m)
            zsum2 = zsum2 + zrhomt(lm,ir,ias)*ylmgq0(lm,ig)
          enddo !m
          !--- i^l*j_l(|G+q0|x)*\sum_m...
          zsum1 = zsum1 + jlgq0r(l,ir)*dconjg(zil(l))*zsum2
        enddo !l
        t1 = rcmt(ir,is)**2
        fr1(ir) = dreal(zsum1)*t1
        fr2(ir) = dimag(zsum1)*t1
      enddo !ir
      call fderiv(-1,nr,rcmt(:,is),fr1,gr,cf)
      t1 = gr(nr)
      call fderiv(-1,nr,rcmt(:,is),fr2,gr,cf)
      t2 = gr(nr)
      zsum = zsum + (fourpi/omega)*dconjg(sfacgq0(ig,ias))*dcmplx(t1,t2)
    enddo !ia
  enddo !is
  
  enddo
  
  end subroutine 

end module
