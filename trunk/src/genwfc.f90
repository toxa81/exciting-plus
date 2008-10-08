subroutine genwfc(ik,lmax,lmmax,mtord,uu,evecfv,evecsv)
use modmain
use modwann
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: mtord
real(8), intent(in) :: uu(0:lmax,mtord,mtord,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)

complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: acoeff(:,:,:,:)
complex(8), allocatable :: prjao(:,:,:)
complex(8), allocatable :: s(:,:)
integer ispn,i,j,n,m1,m2,io1,io2,ias,lm1,lm2,ierr,l
complex(8) zt2(wf_dim,wf_dim)
integer, external :: ikglob

allocate(acoeff(lmmax,mtord,natmtot,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

call match(ngk(ikglob(ik),1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
call getacoeff(lmax,lmmax,ngk(ikglob(ik),1),mtord,apwalm,evecfv,evecsv,acoeff)

! find bands for a given energy interval 
if (wann_use_lhen) then
  do ispn=1,wann_nspins
    do n=1,wf_dim
      wf_lhbnd(1,ispn,n)=1
      do i=1,nstfv
        if (evalsv(i+(ispn-1)*nstfv,ikglob(ik)).lt.wf_lhen(1,ispn,n)) &
          wf_lhbnd(1,ispn,n)=i+1
        if (evalsv(i+(ispn-1)*nstfv,ikglob(ik)).le.wf_lhen(2,ispn,n)) &
          wf_lhbnd(2,ispn,n)=i
      enddo
    enddo    
  enddo
endif

! compute <\psi|g_n>
allocate(prjao(wf_dim,nstfv,wann_nspins))
prjao=dcmplx(0.d0,0.d0)
do n=1,wf_dim
  ias=wf_n(n,1)
  do ispn=1,wann_nspins
     do j=wf_lhbnd(1,ispn,n),wf_lhbnd(2,ispn,n)
      l=wf_n(n,3)
      do m1=-l,l
        lm1=idxlm(l,m1)
        lm2=wf_n(n,2)
        do io1=1,mtord
          io2=2
          prjao(n,j,ispn)=prjao(n,j,ispn)+dconjg(acoeff(lm1,io1,ias,j+(ispn-1)*nstfv)) * &
            uu(l,io1,io2,ias)*ylm2rlm(lm2,lm1)
        enddo !io1
      enddo !m
    enddo !j
  enddo !ispn
enddo !n

allocate(s(wf_dim,wf_dim))
do ispn=1,wann_nspins
! compute ovelap matrix
  s=dcmplx(0.d0,0.d0)
  do m1=1,wf_dim
    do m2=1,wf_dim
      do j=1,nstfv
        s(m1,m2)=s(m1,m2)+prjao(m1,j,ispn)*dconjg(prjao(m2,j,ispn))
      enddo
    enddo
  enddo
! compute S^{-1/2}
  call isqrtzhe(wf_dim,s,ierr)
  if (ierr.ne.0) then
    write(*,*)
    write(*,'("Error(wann_a_ort): faild to calculate S^{-1/2} for spin ",I1,&
      &" during iteration ",I4)')ispn,iscl
    write(*,'("  non-orthogonal WF will be used")')
    write(*,*)
  endif
! compute Wannier function expansion coefficients
  wfc(:,:,ispn,ik)=dcmplx(0.d0,0.d0)
  if (ierr.eq.0) then
    do m1=1,wf_dim
      do m2=1,wf_dim
        wfc(m1,:,ispn,ik)=wfc(m1,:,ispn,ik)+prjao(m2,:,ispn)*dconjg(s(m2,m1))
      enddo
    enddo
  else
    wfc(:,:,ispn,ik)=prjao(:,:,ispn)
  endif
! compute H(k) in WF basis
  zt2=dcmplx(0.d0,0.d0)
  do m1=1,wf_dim
    do m2=1,wf_dim
      do j=1,nstfv
        zt2(m1,m2)=zt2(m1,m2)+dconjg(wfc(m1,j,ispn,ik))*wfc(m2,j,ispn,ik) * &
	      evalsv(j+(ispn-1)*nstfv,ikglob(ik))
      enddo
    enddo
  enddo
  wf_h(:,:,ispn,ikglob(ik))=zt2(:,:)
  call diag_mtrx(wf_dim,zt2,wf_e(1,ispn,ikglob(ik)))
enddo !ispn

deallocate(acoeff,apwalm,prjao,s)

return
end
