subroutine gendmatrsh
use modmain
implicit none
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:)
complex(8), allocatable :: wf1(:,:,:)
complex(8), allocatable :: dmatrsh(:,:,:)
integer ik,is,ia,ias,l,lm1,lm2,m1,m2,io1,io2,ispn,j

integer, external :: ikglob


allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wf1(lmmaxlu,nrfmax,nstsv))
allocate(dmatrsh(lmmaxlu,lmmaxlu,natmtot))
dmatrsh=dcmplx(0.d0,0.d0)

! begin loop over k-points
do ik=1,nkptloc(iproc)
  call match(ngk(ikglob(ik),1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(ikglob(ik),1),evecfvloc(1,1,1,ik),evecsvloc(1,1,ik),apwalm,wfsvmt)
! begin loop over atoms and species
  do is=1,nspecies
    l=llu(is)
    if (l.lt.0) goto 10
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      wf1=dcmplx(0.d0,0.d0)
      do j=1,nstsv
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            do io1=1,nrfmax
              wf1(lm1,io1,j)=wf1(lm1,io1,j)+wfsvmt(lm2,io1,ias,j+(ispn-1)*nstfv)*rlm2ylm(lm2,lm1)
            enddo
          enddo
        enddo
      enddo
          
!      do j=1,nstfv !5,17
!        do m1=-l,l
!          lm1=idxlm(l,m1)
!          do m2=-l,l
!            lm2=idxlm(l,m2)
!            do io1=1,nrfmax
!            do io2=1,nrfmax
!              ispn=1
!              prjao(lm2,j,ispn)=prjao(lm2,j,ispn)+dconjg(wfsvmt(lm1,io1,ias,j+(ispn-1)*nstfv)) * &
!                urfprod(l,io1,io2,ias)*ylm2rlm(lm2,lm1)*wfsvmt(lm1,io2,ias,j+(ispn-1)*nstfv) 
!            enddo
!            enddo !io1
!          enddo !m2
!        enddo !m1
!        do m1=-l,l
!          lm1=idxlm(l,m1)
!          do io1=1,nrfmax
!            do io2=1,nrfmax
!              ispn=1
!              prjao(lm1,j,ispn)=prjao(lm1,j,ispn)+dconjg(wfsvmt(lm1,io1,ias,j+(ispn-1)*nstfv)) * &
!                urfprod(l,io1,io2,ias)*wfsvmt(lm1,io2,ias,j+(ispn-1)*nstfv) 
!            enddo
!          enddo !io1
!        enddo !m1
!      enddo !j
      do j=5,17 !1,nstfv
        do io1=1,nrfmax
        do io2=1,nrfmax
        do m1=-l,l
        do m2=-l,l
          lm1=idxlm(l,m1)
          lm2=idxlm(l,m2)
          ispn=1
          !dmatrsh(lm1,lm2,ias)=dmatrsh(lm1,lm2,ias)+dconjg(wfsvmt(lm1,io1,ias,j+(ispn-1)*nstfv))* &
          !    wfsvmt(lm2,io2,ias,j+(ispn-1)*nstfv)*urfprod(l,io1,io2,ias)*wkpt(ik)*occsv(j,ik)
          dmatrsh(lm1,lm2,ias)=dmatrsh(lm1,lm2,ias)+dconjg(wf1(lm1,io1,j))*wf1(lm2,io2,j)* &
            urfprod(l,io1,io2,ias)*wkpt(ik) !*occsv(j,ik)
        enddo
        enddo
        enddo
        enddo
      enddo !j
    enddo !ia
10 continue
  enddo !is
enddo !ik
call zsync(dmatrsh,lmmaxlu*lmmaxlu*natmtot,.true.,.true.)
if (iproc.eq.0) then
  open(50,file='DMATRSH.OUT',form='formatted',status='replace')
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        write(50,'("ias : ",I2)')ias
        write(50,'(" real part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dreal(dmatrsh(idxlm(l,m1),idxlm(l,m2),ias)),m2=-l,l)  
        enddo
        write(50,'(" imag part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dimag(dmatrsh(idxlm(l,m1),idxlm(l,m2),ias)),m2=-l,l)  
        enddo
      enddo
    endif
  enddo !is
  close(50)
endif  

deallocate(wfsvmt,apwalm,wf1,dmatrsh)

return
end
