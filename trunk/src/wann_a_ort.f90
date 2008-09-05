subroutine wann_a_ort(ik,lmax,lmmax,mtord,uu,evecfv,evecsv)
use modmain
use modwann
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: mtord
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
real(8), intent(in) :: uu(0:lmax,mtord,mtord,natmtot)

complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: acoeff(:,:,:,:,:)
complex(8), allocatable :: acoeff_r(:,:,:,:,:)
complex(8), allocatable :: a_tmp(:,:,:)
complex(8), allocatable :: s(:,:)
integer ispn,j,n,m1,m2,io1,io2,ias,lm1,lm2,ierr,l,ikglob
complex(8) zt1(16)

ikglob=ikptloc(iproc,1)+ik-1

allocate(acoeff(lmmax,mtord,natmtot,nspinor,nstsv))
allocate(acoeff_r(lmmax,mtord,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

call match(ngk(ikglob,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
call getacoeff(lmax,lmmax,ngk(ikglob,1),mtord,apwalm,evecfv,evecsv,acoeff)

do j=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      do io1=1,mtord
        zt1=dcmplx(0.d0,0.d0)
        do lm1=1,16
          do lm2=1,16
            zt1(lm1)=zt1(lm1)+rlm2ylm(lm2,lm1)*acoeff(lm2,io1,ias,ispn,j)
          enddo
        enddo
        acoeff_r(1:16,io1,ias,ispn,j)=zt1(1:16)
      enddo !io1
    enddo !ias
  enddo !ispn
enddo !j

allocate(a_tmp(wf_dim,nstfv,nspinor))
a_tmp=dcmplx(0.d0,0.d0)
do n=1,wf_dim
  ias=wf_n(n,1)
  do ispn=1,nspinor
     do j=wf_lhbnd(n,1),wf_lhbnd(n,2) ! <-spin for bands should be here
      l=wf_n(n,3)
      do m1=-l,l
        lm1=idxlm(l,m1)
        lm2=wf_n(n,2)
        do io1=1,mtord
          io2=1
!          a_tmp(n,j,ispn)=a_tmp(n,j,ispn)+dconjg(acoeff(lm1,io1,ias,ispn,j)) * &
!            acoeff_r(lm2,io2,ias,ispn,j)*uu(l,io1,io2,ias)*ylm2rlm(lm2,lm1)
          a_tmp(n,j,ispn)=a_tmp(n,j,ispn)+dconjg(acoeff(lm1,io1,ias,ispn,j+(ispn-1)*nstfv)) * &
            uu(l,io1,io2,ias)*ylm2rlm(lm2,lm1)
        enddo !io1
      enddo !m
    enddo !j
  enddo !ispn
!  write(*,*)'n=',n
!  do j=1,nstsv
!    write(*,*)'j=',j,'a=',abs(a_tmp(n,j))
!  enddo
enddo !n

a_ort(:,:,:,ikglob)=dcmplx(0.d0,0.d0)
wf_h(:,:,:,ikglob)=dcmplx(0.d0,0.d0)
allocate(s(wf_dim,wf_dim))
do ispn=1,nspinor
  s=dcmplx(0.d0,0.d0)
  do m1=1,wf_dim
    do m2=1,wf_dim
      do j=1,nstfv
        s(m1,m2)=s(m1,m2)+a_tmp(m1,j,ispn)*dconjg(a_tmp(m2,j,ispn))
      enddo
    enddo
!    write(*,'(255F12.6)')(abs(s(m1,m2)),m2=1,wf_dim)
  enddo
  call decomp_ovm(wf_dim,s,ierr)
  if (ierr.ne.0) then
    write(*,*)'Error caclulating S^{-1/2}'
  endif
  do m1=1,wf_dim
    do m2=1,wf_dim
      a_ort(m1,:,ispn,ikglob)=a_ort(m1,:,ispn,ikglob)+a_tmp(m2,:,ispn)*dconjg(s(m2,m1))
    enddo
  enddo
! check orthonormality
  s=dcmplx(0.d0,0.d0)
  do m1=1,wf_dim
    do m2=1,wf_dim
      do j=1,nstfv
        s(m1,m2)=s(m1,m2)+a_ort(m1,j,ispn,ikglob)*dconjg(a_ort(m2,j,ispn,ikglob))
      enddo
      if (m1.eq.m2) then
        if (abs(s(m1,m2)-1.d0).gt.1d-10) then
          write(*,*)'Bad norm'
        endif
      else
        if (abs(s(m1,m2)).gt.1d-10) then
          write(*,*)'Bad orth'
        endif
      endif
    enddo
  enddo
!  do n=1,wf_dim
!    write(*,*)'wf=',n
!    do j=1,nstsv
!      write(*,*)'  j=',j,'a=',abs(a_ort(n,j,1,ik))
!    enddo
!  enddo
  
  do m1=1,wf_dim
    do m2=1,wf_dim
      do j=1,nstfv
        wf_h(m1,m2,ispn,ikglob)=wf_h(m1,m2,ispn,ikglob)+ &
	  dconjg(a_ort(m1,j,ispn,ikglob))*a_ort(m2,j,ispn,ikglob)*evalsv(j+(ispn-1)*nstfv,ikglob)
      enddo
    enddo
!    write(*,'(255F12.6)')(abs(wf_h(m1,m2,ispn,ik)),m2=1,wf_dim)
  enddo
  call diag_mtrx(wf_dim,wf_h(1,1,ispn,ikglob),wf_e(1,ispn,ikglob))
  
    
  
enddo

deallocate(acoeff,acoeff_r,apwalm,a_tmp,s)

return
end
