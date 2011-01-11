subroutine bandchar(tsym,ikloc,lmax,lmmax,evecfv,evecsv,bndchr)
use modmain
implicit none
! arguments
logical, intent(in) :: tsym
integer, intent(in) :: ikloc
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
real(4), intent(out) :: bndchr(lmmax,natmtot,nspinor,nstsv)
! local variables
integer io1,io2,ispn1,ispn2,lm1,lm2,m1,m2
integer l,lm,j,ik,ispn,ias,is,ic
real(8) d1
complex(8) z1
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: dmatylm(:,:,:,:,:)
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
allocate(wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(dmatylm(lmmax,lmmax,nspinor,nspinor,natmtot))
call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),sfacgk(1,1,1,ikloc),apwalm)
call genwfsvmt(lmax,lmmax,ngk(1,ik),evecfv(:,:,1),evecsv,apwalm,wfsvmt)
do j=1,nstsv
  dmatylm=zzero  
  do ias=1,natmtot
    is=ias2is(ias)
    ic=ias2ic(ias)
    do l=0,lmax
      do ispn1=1,nspinor; do ispn2=1,nspinor
        do io1=1,nufr(l,is); do io2=1,nufr(l,is)
          do m1=-l,l; do m2=-l,l
            lm1=idxlm(l,m1)
            lm2=idxlm(l,m2)
            z1=wfsvmt(lm1,io1,ias,ispn1,j)*dconjg(wfsvmt(lm2,io2,ias,ispn2,j))*&
              ufrp(l,io1,io2,ic)
            dmatylm(lm1,lm2,ispn1,ispn2,ias)=dmatylm(lm1,lm2,ispn1,ispn2,ias)+z1
          enddo; enddo !m
        enddo; enddo !io
      enddo; enddo !ispn
    enddo !l
  enddo !ias
  if (tsym) then
    call symdmat(lmax,lmmax,dmatylm)
  endif
! compute matrix in Rlm (in local point symmetry)
  do ias=1,natmtot
    do ispn=1,nspinor
! density matrix is real spherical harmonics (local point symmetry)
      call unimtrxt(lmmax,rylm_lps(1,1,ias),dmatylm(1,1,ispn,ispn,ias))  
      do lm=1,lmmax
        d1=dreal(dmatylm(lm,lm,ispn,ispn,ias))
        bndchr(lm,ias,ispn,j)=real(d1)
      enddo
    enddo !ispn
  enddo !ias
enddo !j
deallocate(wfsvmt,apwalm,dmatylm)
return
end subroutine
!EOC
