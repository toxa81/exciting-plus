subroutine checknorm
use modmain
implicit none

integer ik,ist1,ist2,j
integer ngknr
real(8) t1
complex(8) zt1
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: tpgkcnr(:,:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: wfsvit(:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:)
real(8), allocatable :: wfnrmdev(:)

complex(8), external :: zfint


call init0
call init1
call readstate
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr

allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))

allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv))
allocate(wfsvit(nmatmax,nstsv))
allocate(wfnrmdev(nstsv*(nstsv+1)/2))
open(60,file='NORM.OUT',form='FORMATTED',status='REPLACE')

do ik=1,nkptnr
! generate the G+k vectors
  call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(vklnr(:,ik),vgklnr,evecfv)
  call getevecsv(vklnr(:,ik),evecsv)
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! calculate the wavefunctions for all states
  call genwfsv(.false.,ngknr,igkignr,evalsv,apwalm,evecfv,evecsv,wfmt,wfir)
  
  call genwfsvmt(lmaxvr,lmmaxvr,ngknr,evecfv,evecsv,apwalm,wfsvmt)
  call genwfsvit(ngknr,evecfv,evecsv,wfsvit)
  wfnrmdev=0.d0
  call wfsvprodk(ngknr,igkignr,wfsvmt,wfsvit,wfnrmdev)
    
  j=0
  do ist1=1,nstsv
    do ist2=ist1,nstsv
      j=j+1
      call vnlrho(.true.,wfmt(:,:,:,:,ist1),wfmt(:,:,:,:,ist2), &
        wfir(:,:,ist1),wfir(:,:,ist2),zrhomt,zrhoir)
        zt1=zfint(zrhomt,zrhoir)
      t1=0.d0
      if (ist1.eq.ist2) t1=1.d0
      t1=abs(zt1-t1)
      if (t1.gt.1.d-4) then
        write(60,'("ik : ",I4,4x,"n,n'' : ",2I4,4x,"dev : ",F18.10)')ik,ist1,ist2,t1
      end if
      if (t1.gt.1d-1) then
        write(*,*)
        write(*,'("Warning(checknorm) : very big deviation from norm")')
        write(*,'("ik : ",I4,4x,"n,n'' : ",2I4,4x,"dev : ",F18.10)')ik,ist1,ist2,t1
        write(*,*)
      endif
      write(*,*)wfnrmdev(j),t1
    end do
  end do
enddo
close(60)


return
end
