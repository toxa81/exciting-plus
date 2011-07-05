subroutine gendmatrsh
use modmain
use modldapu
implicit none
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: dmatrlm(:,:,:,:,:)
complex(8), allocatable :: dmatylm(:,:,:,:,:)
complex(8), allocatable :: dmatrlmlps(:,:,:,:,:)
integer is,ia,ias,l,lm1,lm2,m1,m2,io1,io2,ispn,jspn,j,nlps
real(8), allocatable :: mtrx(:,:)
real(8), allocatable :: eval(:)
complex(8) z1
integer i1
real(8) d1(nspinor)
integer ikloc,ik
character*2 :: c2
character*20 :: fmt
logical, external :: bndint

allocate(wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(dmatrlm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dmatylm(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(dmatrlmlps(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))

dmatylm=zzero
! compute matrix in Ylm
! begin loop over k-points
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),sfacgk(1,1,1,ikloc),&
    apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),apwalm,wfsvmt)
! begin loop over atoms and species
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do j=1,nstsv
          if (ldensmtrx.or.(occsv(j,ik).gt.epsocc)) then
            do ispn=1,nspinor
            do jspn=1,nspinor
              do io1=1,nufrmax
              do io2=1,nufrmax
                do m1=-l,l
                do m2=-l,l
                  lm1=idxlm(l,m1)
                  lm2=idxlm(l,m2)
                  z1=wfsvmt(lm1,io1,ias,ispn,j)*&
                    dconjg(wfsvmt(lm2,io2,ias,jspn,j))*&
                    ufrp(l,io1,io2,ias2ic(ias))*wkpt(ik)
                  if (ldensmtrx.and.bndint(j,evalsv(j,ik),dm_e1,dm_e2)) then
                    dmatylm(lm1,lm2,ispn,jspn,ias)=dmatylm(lm1,lm2,ispn,jspn,ias)+z1
                  else
                    dmatylm(lm1,lm2,ispn,jspn,ias)=dmatylm(lm1,lm2,ispn,jspn,ias)+&
                      z1*occsv(j,ik)
                  endif
                enddo
                enddo
              enddo
              enddo
            enddo !ispn
            enddo !jspn
          endif
        enddo !j
      enddo !ia
    endif
  enddo !is
enddo !ikloc
call mpi_grid_reduce(dmatylm(1,1,1,1,1),lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,&
  dims=(/dim_k/),all=.true.)
call symdmat(lmaxlu,lmmaxlu,dmatylm)
dmatlu=dmatylm
! compute matrix in Rlm
do ias=1,natmtot
  do ispn=1,nspinor
    do jspn=1,nspinor
! density matrix in real spherical harmonics
      dmatrlm(:,:,ispn,jspn,ias)=dmatylm(:,:,ispn,jspn,ias)
      call unimtrxt(lmmaxlu,rylm,dmatrlm(1,1,ispn,jspn,ias))
! cut numerical noise
      do lm1=1,lmmaxlu
        do lm2=1,lmmaxlu
          i1=int(dreal(dmatrlm(lm1,lm2,ispn,jspn,ias))*1.0d8)
          dmatrlm(lm1,lm2,ispn,jspn,ias)=dcmplx(i1*1.0d-8,0.d0)
        enddo
      enddo
! density matrix is real spherical harmonics (local point symmetry)
      dmatrlmlps(:,:,ispn,jspn,ias)=dmatylm(:,:,ispn,jspn,ias)
      call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),dmatrlmlps(1,1,ispn,jspn,ias))  
    enddo
  enddo
enddo

if (mpi_grid_root()) then
  open(50,file="DMATRSH.OUT",form="FORMATTED",status="REPLACE")
  if (ldensmtrx) then
    write(50,'("Density matrix")')
    write(50,'("  band interval (Ha) : ",2F8.2)')dm_e1,dm_e2
  else
    write(50,'("Occupancy matrix")')
  endif
  nlps=0
  do is=1,nspecies
    l=llu(is)
    if (l.gt.0) then
      allocate(mtrx(2*l+1,2*l+1))
      allocate(eval(2*l+1))
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        nlps=nlps+1
! construct format for output
        write(c2,'(I2)')2*l+1
        fmt=trim(c2)//"F12.6"
        if (spinpol) then
          fmt=trim(fmt)//",4X,"//trim(fmt)
        endif  
        fmt="("//trim(fmt)//")"
! write
        write(50,'("ias : ",I2)')ias
        write(50,'("  in real harmonics (global basis) ")')
        do ispn=1,nspinor
          do m1=-l,l
            write(50,trim(fmt)) &
              ((dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)),m2=-l,l),&
              jspn=1,nspinor)
          enddo
          write(50,*)
        enddo
        write(50,'("  in real harmonics (local basis) ")')
        do ispn=1,nspinor
          do m1=-l,l
            write(50,trim(fmt)) &
              ((dreal(dmatrlmlps(idxlm(l,m1),idxlm(l,m2),ispn,jspn,ias)), &
              m2=-l,l),jspn=1,nspinor)
          enddo
          write(50,*)
        enddo
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
        write(50,*)
      enddo !ia
      deallocate(mtrx,eval)
    endif
  enddo !is
  if (.not.spinpol) then
    write(50,*)
    write(50,'("add to input file : ")')
    write(50,*)
    write(50,'("lps")')
    write(50,'(I4)')nlps
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
!    write(50,*)
!    write(50,'("transponse of eigen-vector matrix (for easy copy-paste) : ")')
!    write(50,*)
!    do is=1,nspecies
!      l=llu(is)
!      if (l.gt.0) then
!        allocate(mtrx(2*l+1,2*l+1))
!        allocate(eval(2*l+1))
!        do ia=1,natoms(is)
!          ias=idxas(ia,is)
!          write(50,'(2I4)')ias,l
!          do m1=-l,l
!            do m2=-l,l
!              mtrx(m1+l+1,m2+l+1)=dreal(dmatrlm(idxlm(l,m1),idxlm(l,m2),1,1,ias))
!            enddo
!          enddo
!          call diagdsy(2*l+1,mtrx,eval)
!          do m1=1,2*l+1
!            write(50,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
!          enddo
!        enddo !ia
!        deallocate(mtrx,eval)
!      endif
!    enddo !is
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
        deallocate(mtrx,eval)
      endif
    enddo
  endif
  close(50)
endif  

deallocate(wfsvmt,apwalm,dmatrlm)
deallocate(dmatrlmlps,dmatylm)

return
end
