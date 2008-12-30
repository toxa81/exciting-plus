subroutine wann_init1
use modmain
implicit none

integer ik,n,i
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)

call init0
call init1

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr

call geturf

allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(wfsvit(nmatmax,nstsv,nspinor))

allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,wf_dim,nspinor,nkpt))
allocate(wann_unkit(nmatmax,wf_dim,nspinor,nkpt))
wann_unkmt=dcmplx(0.d0,0.d0)
wann_unkit=dcmplx(0.d0,0.d0)

! read and transform eigen-vectors
do ik=1,nkpt
  call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
  call getwfc(ik,wfc(1,1,1,ik))
  call match(ngk(1,ik),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmt)
  call genwfsvit(ngk(1,ik),evecfv,evecsv,wfsvit)
  do n=1,wf_dim
    do i=1,nstfv
      wann_unkmt(:,:,:,n,1,ik)=wann_unkmt(:,:,:,n,1,ik)+wfsvmt(:,:,:,i,1)*wfc(n,i,1,ik)
      wann_unkit(:,n,1,ik)=wann_unkit(:,n,1,ik)+wfsvit(:,i,1)*wfc(n,i,1,ik)
    enddo
  enddo
enddo !ik

return
end

subroutine wann_unk(n,vpl,vrc,val)
use modmain
implicit none
integer, intent(in) :: n
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)

integer isym,ik,i1,i2,i3,np2,is,ia,ias,ir,ir0,io,l,j,m,i,lm,ig
integer ntr(3)
real(8) vpc(3),rmt2,pos(3),v1(3),t1,tr(3),tp(2)
complex(8) zt1,zt2,ylm(lmmaxvr)
logical l1
real(8) ya(nprad),c(nprad),t2
real(8), external :: polynom

call findkpt(vpl,isym,ik)
call r3mv(bvec,vpl,vpc)
call getntr(vrc,ntr)

np2=nprad/2
zt1=dcmplx(0.d0,0.d0)
zt2=dcmplx(0.d0,0.d0)

l1=.false.
! check if point is in a muffin-tin
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do i1=ntr(1)-1,ntr(1)+1
    do i2=ntr(2)-1,ntr(2)+1
    do i3=ntr(3)-1,ntr(3)+1
      tr(:)=i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      pos(:)=atposc(:,ia,is)+tr(:)
      v1(:)=vrc(:)-pos(:)
      t1=v1(1)**2+v1(2)**2+v1(3)**2
      if (t1.lt.rmt2) then
        call sphcrd(v1,t1,tp)
        call genylm(lmaxvr,tp,ylm)
        do ir=1,nrmt(is)
          if (spr(ir,is).ge.t1) then
            if (ir.le.np2) then
              ir0=1
            else if (ir.gt.nrmt(is)-np2) then
              ir0=nrmt(is)-nprad+1
            else
              ir0=ir-np2
            end if
            t1=max(t1,spr(1,is))
            do io=1,nrfmax
              do l=0,lmaxvr
                do j=1,nprad
                  i=ir0+j-1
                  ya(j)=urf(i,l,io,ias)
                end do
                t2=polynom(0,nprad,spr(ir0,is),ya,c,t1)
                do m=-l,l
                  lm=idxlm(l,m)
                  zt1=zt1+t2*ylm(lm)*wann_unkmt(lm,io,ias,n,1,ik)
                 end do !m
              end do  !l
            enddo !io
            zt1=zt1*exp(dcmplx(0.d0,dot_product(vpc,tr(:))))
            l1=.true.
            goto 10
          end if
        end do !ir
      end if
    end do
    end do
    end do
  end do
end do !is
10 continue
! otherwise use interstitial function
if (.not.l1) then
  do ig=1,ngk(1,ik)
!    t1=vgc(1,ig)*vrc(1)+vgc(2,ig)*vrc(2)+vgc(3,ig)*vrc(3)
    t1=vgkc(1,ig,1,ik)*vrc(1)+vgkc(2,ig,1,ik)*vrc(2)+vgkc(3,ig,1,ik)*vrc(3)
    zt2=zt2+cmplx(cos(t1),sin(t1),8)*wann_unkit(ig,n,1,ik)/sqrt(omega)
  enddo
endif
val(1)=dreal(zt1+zt2)
val(2)=dimag(zt1+zt2)

return
end
