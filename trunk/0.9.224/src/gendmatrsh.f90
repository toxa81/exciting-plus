subroutine gendmatrsh
use modmain
implicit none
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wf1(:,:,:)
complex(8), allocatable :: dmatrsh(:,:,:,:,:)
complex(8), allocatable :: dmatrshlcs(:,:,:,:,:)
integer ik,is,ia,ias,l,lm1,lm2,lm3,m1,m2,io1,io2,ispn,j,n1,n2,isym,lspl
real(8), allocatable :: mtrx(:,:)
real(8), allocatable :: eval(:)
complex(8), allocatable :: zflm(:,:)
complex(8), allocatable :: ulm(:,:,:)
complex(8), allocatable :: ulmrsh(:,:,:)
complex(8), allocatable :: dm1(:,:,:,:,:)
complex(8), allocatable :: dm2(:,:)
complex(8), allocatable :: dm3(:,:)
integer, external :: ikglob

allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wf1(lmmaxlu,nrfmax,nstsv))
allocate(dmatrsh(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dmatrshlcs(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dm1(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dm2(lmmaxlu,lmmaxlu))
dm1=dcmplx(0.d0,0.d0)

ispn=1

! begin loop over k-points
do ik=1,nkptloc(iproc)
  call match(ngk(1,ikglob(ik)),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ikglob(ik)),evecfvloc(1,1,1,ik),evecsvloc(1,1,ik),apwalm,wfsvmt)
! begin loop over atoms and species
  do is=1,nspecies
    l=llu(is)
    if (l.lt.0) goto 10
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      wf1(1:lmmaxlu,:,:)=wfsvmt(1:lmmaxlu,:,ias,:,ispn)
      do j=dmbnd1,dmbnd2
        do io1=1,nrfmax
        do io2=1,nrfmax
        do m1=-l,l
        do m2=-l,l
          lm1=idxlm(l,m1)
          lm2=idxlm(l,m2)
          dm1(lm1,lm2,ispn,ispn,ias)=dm1(lm1,lm2,ispn,ispn,ias)+ &
            wf1(lm1,io1,j)*dconjg(wf1(lm2,io2,j))*urfprod(l,io1,io2,ias)* &
            wkpt(ik)*occsv(j,ik)
        enddo
        enddo
        enddo
        enddo
      enddo !j
    enddo !ia
10 continue
  enddo !is
enddo !ik

call zsync(dm1,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
call symdmat(lmaxlu,lmmaxlu,dm1)

dmatrsh=dcmplx(0.d0,0.d0)
dmatrshlcs=dcmplx(0.d0,0.d0)
do ias=1,natmtot
! dens.mtrx. in GCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm(lm1,lm3))*dm1(lm3,lm2,ispn,ispn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dmatrsh(lm1,lm2,ispn,ispn,ias)=dmatrsh(lm1,lm2,ispn,ispn,ias) + &
          dm2(lm1,lm3)*rylm(lm2,lm3)
      enddo
    enddo
  enddo
! dens.mtrx. in LCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm_lcs(lm1,lm3,ias))*dm1(lm3,lm2,ispn,ispn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dmatrshlcs(lm1,lm2,ispn,ispn,ias)=dmatrshlcs(lm1,lm2,ispn,ispn,ias) + &
          dm2(lm1,lm3)*rylm_lcs(lm2,lm3,ias)
      enddo
    enddo
  enddo
enddo

if (iproc.eq.0) then
  open(50,file='DMATRSH.OUT',form='formatted',status='replace')
  write(50,'("bands for density matrix : ",2I4)')dmbnd1,dmbnd2
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      allocate(mtrx(2*l+1,2*l+1))
      allocate(eval(2*l+1))
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do m1=-l,l
          do m2=-l,l
            mtrx(m1+l+1,m2+l+1)=dreal(dmatrsh(idxlm(l,m1),idxlm(l,m2),1,1,ias))
          enddo
        enddo
        call diagdsy(2*l+1,mtrx,eval)
        write(50,'("ias : ",I2)')ias
	write(50,'(" GCS :")')
        write(50,'(" real part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dreal(dmatrsh(idxlm(l,m1),idxlm(l,m2),1,1,ias)),m2=-l,l)  
        enddo
        write(50,'(" imag part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dimag(dmatrsh(idxlm(l,m1),idxlm(l,m2),1,1,ias)),m2=-l,l)  
        enddo
        write(50,'(" eigen-vectors : ")')
        do m1=1,2*l+1
          write(50,'(14G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
        enddo
        write(50,'(" eigen-values : ")')
        write(50,'(14F12.6)')(eval(m1),m1=1,2*l+1)
	write(50,'(" LCS :")')
	write(50,'(" real part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dreal(dmatrshlcs(idxlm(l,m1),idxlm(l,m2),1,1,ias)),m2=-l,l)  
        enddo
        write(50,'(" imag part : ")')
        do m1=-l,l
          write(50,'(14F12.6)')(dimag(dmatrshlcs(idxlm(l,m1),idxlm(l,m2),1,1,ias)),m2=-l,l)  
        enddo

      enddo !ia
      deallocate(mtrx,eval)
    endif
  enddo !is
  close(50)
endif  

deallocate(wfsvmt,apwalm,wf1,dmatrsh)
deallocate(dmatrshlcs,dm1,dm2)

return
end
