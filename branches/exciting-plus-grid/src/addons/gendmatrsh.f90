subroutine gendmatrsh
use modmain
implicit none
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wf1(:,:,:)
complex(8), allocatable :: dmatrlm(:,:,:,:,:)
complex(8), allocatable :: dmatylm(:,:,:,:,:)
complex(8), allocatable :: dmatrlmlcs(:,:,:,:,:)
complex(8), allocatable :: ematrlm(:,:,:,:,:)
complex(8), allocatable :: ematylm(:,:,:,:,:)
complex(8), allocatable :: ematrlmlcs(:,:,:,:,:)
integer is,ia,ias,l,lm1,lm2,lm3,m1,m2,io1,io2,ispn,jspn,j,nlcs
real(8), allocatable :: mtrx(:,:)
real(8), allocatable :: eval(:)
complex(8), allocatable :: dm2(:,:)
complex(8) z1
integer i1
real(8) d1(nspinor)
integer ikloc
character*2 :: c2
character*20 :: fmt

allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wf1(lmmaxlu,nrfmax,nstsv))
allocate(dmatrlm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dmatylm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dmatrlmlcs(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(ematrlm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(ematylm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(ematrlmlcs(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dm2(lmmaxlu,lmmaxlu))

dmatylm=dcmplx(0.d0,0.d0)
ematylm=dcmplx(0.d0,0.d0)
! compute matrix in Ylm
! begin loop over k-points
do ikloc=1,nkptloc
  call match(ngk(1,ikloc),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc), &
    sfacgk(1,1,1,ikloc),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ikloc),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),apwalm,wfsvmt)
! begin loop over atoms and species
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do j=1,nstsv
          do ispn=1,nspinor
          do jspn=1,nspinor
            do io1=1,nrfmax
            do io2=1,nrfmax
              do m1=-l,l
              do m2=-l,l
                lm1=idxlm(l,m1)
                lm2=idxlm(l,m2)
                z1=wfsvmt(lm1,io1,ias,ispn,j)*dconjg(wfsvmt(lm2,io2,ias,jspn,j))*&
                  urfprod(l,io1,io2,ias)*wkpt(ikloc)
                if (ldensmtrx) then
                  if (evalsv(j,ikloc).ge.dm_e1.and. &
                      evalsv(j,ikloc).le.dm_e2) then
                    dmatylm(lm1,lm2,ispn,jspn,ias)=dmatylm(lm1,lm2,ispn,jspn,ias)+z1
                  endif                 
                else
                  dmatylm(lm1,lm2,ispn,jspn,ias)=dmatylm(lm1,lm2,ispn,jspn,ias)+&
                    z1*occsv(j,ikloc)
                endif
                ematylm(lm1,lm2,ispn,jspn,ias)=ematylm(lm1,lm2,ispn,jspn,ias)+&
                    z1*evalsv(j,ikloc)
              enddo
              enddo
            enddo
            enddo
          enddo !ispn
          enddo !jspn
        enddo !j
      enddo !ia
    endif
  enddo !is
enddo !ikloc

!call zsync(dmatylm,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
!call zsync(ematylm,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
call symdmat(lmaxlu,lmmaxlu,dmatylm)
call symdmat(lmaxlu,lmmaxlu,ematylm)

! compute matrix in Rlm
dmatrlm=dcmplx(0.d0,0.d0)
dmatrlmlcs=dcmplx(0.d0,0.d0)
ematrlm=zzero
ematrlmlcs=zzero
do ispn=1,nspinor
do jspn=1,nspinor
do ias=1,natmtot
! dens.mtrx. in GCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm(lm1,lm3))*dmatylm(lm3,lm2,ispn,jspn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dmatrlm(lm1,lm2,ispn,jspn,ias)=dmatrlm(lm1,lm2,ispn,jspn,ias) + &
          dm2(lm1,lm3)*rylm(lm2,lm3)
      enddo
    enddo
  enddo
! e.mtrx. in GCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm(lm1,lm3))*ematylm(lm3,lm2,ispn,jspn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        ematrlm(lm1,lm2,ispn,jspn,ias)=ematrlm(lm1,lm2,ispn,jspn,ias) + &
          dm2(lm1,lm3)*rylm(lm2,lm3)
      enddo
    enddo
  enddo
! dens.mtrx. in LCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm_lcs(lm1,lm3,ias))*dmatylm(lm3,lm2,ispn,jspn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dmatrlmlcs(lm1,lm2,ispn,jspn,ias)=dmatrlmlcs(lm1,lm2,ispn,jspn,ias) + &
          dm2(lm1,lm3)*rylm_lcs(lm2,lm3,ias)
      enddo
    enddo
  enddo
! e.mtrx. in LCS in RSH
  dm2=dcmplx(0.d0,0.d0)
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        dm2(lm1,lm2)=dm2(lm1,lm2)+dconjg(rylm_lcs(lm1,lm3,ias))*ematylm(lm3,lm2,ispn,jspn,ias)
      enddo
    enddo
  enddo
  do lm1=1,lmmaxlu
    do lm2=1,lmmaxlu
      do lm3=1,lmmaxlu
        ematrlmlcs(lm1,lm2,ispn,jspn,ias)=ematrlmlcs(lm1,lm2,ispn,jspn,ias) + &
          dm2(lm1,lm3)*rylm_lcs(lm2,lm3,ias)
      enddo
    enddo
  enddo
enddo
enddo
enddo

! cut numerical noise
do ias=1,natmtot
  do ispn=1,nspinor
    do jspn=1,nspinor
      do lm1=1,lmmaxlu
        do lm2=1,lmmaxlu
          i1=int(dreal(dmatrlm(lm1,lm2,ispn,jspn,ias))*1.0d8)
          dmatrlm(lm1,lm2,ispn,jspn,ias)=dcmplx(i1*1.0d-8,0.d0)
        enddo
      enddo
    enddo
  enddo
enddo
    

if (iproc.eq.0) then
  open(50,file='DMATRSH.OUT',form='formatted',status='replace')
  if (ldensmtrx) then
    write(50,'("Density matrix")')
    write(50,'("  band interval (Ha) : ",2F8.2)')dm_e1,dm_e2
  else
    write(50,'("Occupancy matrix")')
  endif
  nlcs=0
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      allocate(mtrx(2*l+1,2*l+1))
      allocate(eval(2*l+1))
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        nlcs=nlcs+1
        write(c2,'(I2)')2*l+1
        fmt=trim(c2)//"F12.6"
        if (spinpol) then
          fmt=trim(fmt)//",4X,"//trim(fmt)
        endif  
        fmt="("//trim(fmt)//")"
        write(50,'("ias : ",I2)')ias
!        write(50,'("  in complex harmonics (global basis) ")')
!        write(50,'("  real part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dreal(dmatylm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
!        write(50,'("  imag part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dimag(dmatylm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
        write(50,'("  in real harmonics (global basis) ")')
        write(50,'("  real part : ")')
        do ispn=1,nspinor
          do m1=-l,l
          write(50,trim(fmt)) &
            ((dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
              m2=-l,l),jspn=1,nspinor)
          enddo
          write(50,*)
        enddo
!        write(50,'("  imag part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dimag(dmatrlm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
        write(50,'("  in real harmonics (local basis) ")')
        write(50,'("  real part : ")')
        do ispn=1,nspinor
          do m1=-l,l
          write(50,trim(fmt)) &
            ((dreal(dmatrlmlcs(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
              m2=-l,l),jspn=1,nspinor)
          enddo
          write(50,*)
        enddo
!        write(50,'("  imag part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dimag(dmatrlmlcs(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
        write(50,'("  eigen vectors and values : ")')
        do ispn=1,nspinor
        do jspn=1,nspinor
        write(50,'("  ispn jspn : ",2I2)')ispn,jspn
        do m1=-l,l
          do m2=-l,l
            mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias))
          enddo
        enddo
        call diagdsy(2*l+1,mtrx,eval)
        do m1=1,2*l+1
          write(50,'(2X,7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
        enddo
        write(50,*)
        write(50,'(2X,7G18.10)')(eval(m1),m1=1,2*l+1)
        write(50,*)
        enddo
        enddo
        if (.not.ldensmtrx) then
          d1=0.d0
          do ispn=1,nspinor
            do m1=-l,l
              d1(ispn)=d1(ispn)+dreal(dmatrlm(idxlm(l,m1),idxlm(l,m1),ispn,ispn,ias))
            enddo
            write(50,'("  spin : ",I1,"  occupancy : ",F12.6)')ispn,d1(ispn)
          enddo
          write(50,'("     total  occupancy : ",F12.6)')sum(d1)
          if (spinpol) write(50,'("     spin moment : ",F12.6)')d1(1)-d1(2)
          write(50,*)
        endif

!        write(50,'("Energy matrix in real harmonics (local basis) ")')
!        write(50,'("  real part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dreal(ematrlmlcs(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
!        write(50,'("  imag part : ")')
!        do ispn=1,nspinor
!          do m1=-l,l
!          write(50,trim(fmt)) &
!            ((dimag(ematrlmlcs(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
!              m2=-l,l),jspn=1,nspinor)
!          enddo
!          write(50,*)
!        enddo
!       d1=0.d0
!        do ispn=1,nspinor
!          do m1=-l,l
!            d1(ispn)=d1(ispn)+dreal(ematrlmlcs(idxlm(l,m1),idxlm(l,m1),ispn,ispn,ias))
!          enddo
!          write(50,'("  spin : ",I1,"  E_l: ",F12.6)')ispn,d1(ispn)/(2*l+1)
!        enddo
        write(50,*)
      enddo !ia
      deallocate(mtrx,eval)
    endif
  enddo !is
  if (.not.spinpol) then
    write(50,*)
    write(50,'("add to input file : ")')
    write(50,*)
    write(50,'("lcs")')
    write(50,'(I4)')nlcs
    do is=1,nspecies
      l=llu(is)
      if (l.gt.0) then
        allocate(mtrx(2*l+1,2*l+1))
        allocate(eval(2*l+1))
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          write(50,'(2I4)')ias,l
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
          enddo
        enddo !ia
        deallocate(mtrx,eval)
      endif
    enddo !is
    write(50,*)
    write(50,'("transponse of eigen-vector matrix (for easy copy-paste) : ")')
    write(50,*)
    do is=1,nspecies
      l=llu(is)
      if (l.gt.0) then
        allocate(mtrx(2*l+1,2*l+1))
        allocate(eval(2*l+1))
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          write(50,'(2I4)')ias,l
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
          enddo
        enddo !ia
        deallocate(mtrx,eval)
      endif
    enddo !is
  else
    write(50,*)
    write(50,'("spin-polarized case : ")')
    do is=1,nspecies
      l=llu(is)
      if (l.gt.0) then
        allocate(mtrx(2*l+1,2*l+1))
        allocate(eval(2*l+1))
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          write(50,'("ias : ",I4,"     l : ",I1)')ias,l
! uu
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          write(50,'("eigen vectors and values of uu matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
          enddo
          write(50,*)
          write(50,'(2X,7G18.10)')(eval(m1),m1=1,2*l+1)
          write(50,*)
          write(50,'("transponse eigen vectors of uu matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
          enddo
          write(50,*)
! dd
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),2,2,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          write(50,'("eigen vectors and values of dd matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
          enddo
          write(50,*)
          write(50,'(2X,7G18.10)')(eval(m1),m1=1,2*l+1)
          write(50,*)
          write(50,'("transponse eigen vectors of dd matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
          enddo
          write(50,*)
! uu+dd
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))+&
                dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),2,2,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          write(50,'("eigen vectors and values of uu+dd matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
          enddo
          write(50,*)
          write(50,'(2X,7G18.10)')(eval(m1),m1=1,2*l+1)
          write(50,*)
          write(50,'("transponse eigen vectors of uu+dd matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
          enddo
          write(50,*)
! dd-uu
          do m1=-l,l
            do m2=-l,l
              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),2,2,ias))-&
                dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))
            enddo
          enddo
          call diagdsy(2*l+1,mtrx,eval)
          write(50,'("eigen vectors and values of dd-uu matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
          enddo
          write(50,*)
          write(50,'(2X,7G18.10)')(eval(m1),m1=1,2*l+1)
          write(50,*)
          write(50,'("transponse eigen vectors of dd-uu matrix : ")')
          do m1=1,2*l+1
            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
          enddo
          write(50,*)
        enddo
      endif
    enddo
  endif
  close(50)
endif  

deallocate(wfsvmt,apwalm,wf1,dmatrlm)
deallocate(dmatrlmlcs,dmatylm,dm2)

return
end
