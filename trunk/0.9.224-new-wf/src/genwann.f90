subroutine genwann(ik,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)

complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: prjao(:,:,:)
complex(8), allocatable :: s(:,:)
integer ispn,i,j,n,m1,m2,io1,io2,ias,lm1,lm2,ierr,l
complex(8), allocatable :: zt2(:,:)
complex(8), allocatable :: wfcnew(:,:)
integer, external :: ikglob

allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

call match(ngk(1,ikglob(ik)),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ikglob(ik)),evecfv,evecsv,apwalm,wfsvmt)
call genwann2(evalsv(1,ikglob(ik)),wfsvmt,wann_c(1,1,1,ik))
do ispn=1,wann_nspin
  allocate(zt2(nwann(ispn),nwann(ispn)))
! compute H(k) in WF basis
  zt2=dcmplx(0.d0,0.d0)
  do m1=1,nwann(ispn)
    do m2=1,nwann(ispn)
      do j=1,nstfv
        zt2(m1,m2)=zt2(m1,m2)+dconjg(wann_c(m1,j,ispn,ik))*wann_c(m2,j,ispn,ik) * &
	      evalsv(j+(ispn-1)*nstfv,ikglob(ik))
      enddo
    enddo
  enddo
  wf_h(1:nwann(ispn),1:nwann(ispn),ispn,ikglob(ik))=zt2(:,:)
  call diagzhe(nwann(ispn),zt2,wf_e(1,ispn,ikglob(ik)))
  deallocate(zt2)
enddo !ispn

deallocate(wfsvmt,apwalm)

return
end

subroutine genwann2(e,wfsvmt,wf)
use modmain
implicit none
! arguments
real(8), intent(in) :: e(nstsv)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(out) :: wf(wann_nmax,nstfv,wann_nspin)

complex(8), allocatable :: prjao(:,:,:)
complex(8), allocatable :: s(:,:)
integer ispn,i,j,n,m1,m2,io1,io2,ias,lm1,lm2,ierr,l,itype
integer n1n2(2,2,wann_ntype),n1,n2

n1n2=0
! find bands for a given energy interval 
if (wann_use_eint) then
  do ispn=1,wann_nspin
    do j=1,wann_ntype
      n1n2(1,ispn,j)=1
      do i=1,nstfv
        if (e(i+(ispn-1)*nstfv).lt.wann_eint(1,j)) n1n2(1,ispn,j)=i+1
        if (e(i+(ispn-1)*nstfv).le.wann_eint(2,j)) n1n2(2,ispn,j)=i
      enddo
    enddo    
  enddo
endif

! compute <\psi|g_n>
allocate(prjao(wann_nmax,nstfv,wann_nspin))
prjao=dcmplx(0.d0,0.d0)
do ispn=1,wann_nspin
  do n=1,nwann(ispn)
    ias=iwann(n,ispn,1)
    itype=iwann(n,ispn,4)
    if (wann_use_eint) then
      n1=n1n2(1,ispn,itype)
      n2=n1n2(2,ispn,itype)
    else
      n1=wann_nint(1,itype)
      n2=wann_nint(2,itype)
    endif
    do j=n1,n2
      l=iwann(n,ispn,3)
      do m1=-l,l
        lm1=idxlm(l,m1)
        lm2=iwann(n,ispn,2)
        do io1=1,nrfmax
          io2=2
          prjao(n,j,ispn)=prjao(n,j,ispn)+dconjg(wfsvmt(lm1,io1,ias,j+(ispn-1)*nstfv,ispn)) * &
            urfprod(l,io1,io2,ias)*rylm_lcs(lm2,lm1,ias)
        enddo !io1
      enddo !m
    enddo !j
  enddo !n
enddo !ispn


do ispn=1,wann_nspin
  allocate(s(nwann(ispn),nwann(ispn)))
! compute ovelap matrix
  s=dcmplx(0.d0,0.d0)
  do m1=1,nwann(ispn)
    do m2=1,nwann(ispn)
      do j=1,nstfv
        s(m1,m2)=s(m1,m2)+prjao(m1,j,ispn)*dconjg(prjao(m2,j,ispn))
      enddo
    enddo
  enddo
! compute S^{-1/2}
  call isqrtzhe(nwann(ispn),s,ierr)
  if (ierr.ne.0) then
    write(*,*)
    write(*,'("Error(genwann2): failed to calculate S^{-1/2} for spin ",I1)')ispn
    do n=1,nwann(ispn)
      write(*,*)'  prjao=',abs(prjao(n,:,ispn))
    enddo
    write(*,*)
  endif
! compute Wannier function expansion coefficients
  wf(:,:,ispn)=dcmplx(0.d0,0.d0)
  if (ierr.eq.0) then
    do m1=1,nwann(ispn)
      do m2=1,nwann(ispn)
        wf(m1,:,ispn)=wf(m1,:,ispn)+prjao(m2,:,ispn)*dconjg(s(m2,m1))
      enddo
    enddo
  else
    wf(:,:,ispn)=prjao(:,:,ispn)
  endif
  deallocate(s)
enddo !ispn

deallocate(prjao)

return
end
