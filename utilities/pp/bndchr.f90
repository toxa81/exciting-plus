program main
implicit none
      
integer lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
integer lm,ias,ispn,ist,ik
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: dpp1d(:),e(:,:),wt(:,:),lines(:,:)
integer, allocatable :: orb(:,:,:)
integer i,j,ist1,n,iin
real(8) scale,emin,emax,wtmax,efermi
real(8), parameter :: ha2ev = 27.21138386d0
logical wannier,l1
integer nwann
real(8), allocatable :: wann_c(:,:,:)
character*2 e_units
character*20 strin
integer nspecies,ias1,is1,is
character(256), allocatable :: spsymb(:)
integer, allocatable :: natoms(:)
integer, allocatable :: iasis(:)

      
open(50,file="BNDCHR.OUT",form="FORMATTED",status="OLD")
read(50,*)lmmax,nspecies,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
read(50,*)efermi
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(dpp1d(nkpt))
allocate(e(nstsv,nkpt))
allocate(spsymb(nspecies))
allocate(natoms(nspecies))
allocate(iasis(natmtot))
do is=1,nspecies
  read(50,*)spsymb(is)
  read(50,*)natoms(is)
enddo
do ias=1,natmtot
  read(50,*)ias1,is1
  iasis(ias1)=is1
enddo      
do ik=1,nkpt
  read(50,*)dpp1d(ik)
  read(50,*)(e(ist,ik),ist=1,nstsv)
  read(50,*)((((bndchr(lm,ias,ispn,ist,ik),lm=1,lmmax), &
                ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
enddo
!read(50,*)wannier
!if (wannier) then
!  read(50,*)nwann
!  allocate(wann_c(nwann,nstsv,nkpt))
!  do ik=1,nkpt
!    read(50,*)((wann_c(n,i,ik),n=1,nwann),i=1,nstsv)
!  enddo
!endif
close(50)

35 continue
write(*,'("Put Fermi level to zero? (1-Yes,2-No)")')
read(*,*)iin
if (.not.(iin.eq.1.or.iin.eq.2)) goto 35
if (iin.eq.1) then
  e=e-efermi
  efermi=0.d0
endif

40 continue
write(*,'("Energy units? (1-eV,2-Ha)")')
read(*,*)iin
if (.not.(iin.eq.1.or.iin.eq.2)) goto 40
e_units="Ha"
if (iin.eq.1) then
  e=e*ha2ev
  e_units="eV"
endif
emin=minval(e(:,:))
emax=maxval(e(:,:))
write(*,'("Band energy range : ",2F8.2," [",A,"]")')emin,emax,e_units
write(*,'("Input energy interval (emin emax) [",A,"]")')e_units
read(*,*)emin,emax
     

allocate(wt(nstsv,nkpt))
wt=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    do lm=1,lmmax
      do ias=1,natmtot
        do ispn=1,nspinor
          wt(ist,ik)=wt(ist,ik)+bndchr(lm,ias,ispn,ist,ik)
        enddo
      enddo
    enddo
  enddo
enddo
wtmax=maxval(wt(:,:))
write(*,*)
!write(*,'("Maximum weight of all orbitals : ",G18.10)')wtmax
bndchr(:,:,:,:,:)=bndchr(:,:,:,:,:)/wtmax
!write(*,'("Band-character weight is normalized")')


allocate(orb(lmmax,nspinor,natmtot))
orb=0
22 continue
write(*,'("Input orbitals for band character")')
read(*,'(A)')strin
if (trim(strin).eq.'q') goto 23
call addorb(strin,orb,lmmax,nspinor,natmtot,iasis,spsymb)
goto 22
23 continue 

wt=0.d0
do ispn=1,nspinor
  write(*,'("spin : ",I1)')ispn
  do ias=1,natmtot
    write(*,'("  atom : ",I4," orb : ",16I2)')ias,(orb(lm,ispn,ias),lm=1,lmmax)
    do ik=1,nkpt
      do ist=1,nstsv
        do lm=1,lmmax
          wt(ist,ik)=wt(ist,ik)+orb(lm,ispn,ias)*bndchr(lm,ias,ispn,ist,ik)
        enddo
      enddo
    enddo
  enddo
enddo
      
scale=1.d0
write(*,'("Input maximum line width [",A,"]")')e_units
read(*,*)scale
      
open(50,file='BNDS.DAT',form='formatted',status='replace')
open(51,file='BNDS1.DAT',form='formatted',status='replace')
open(52,file='BNDS2.DAT',form='formatted',status='replace')
open(53,file='BNDS3.DAT',form='formatted',status='replace')
open(54,file='BNDS_xydy.DAT',form='formatted',status='replace')
do ist1=1,nstsv
  do ik=1,nkpt
    write(50,*)dpp1d(ik),e(ist1,ik)
    write(51,*)dpp1d(ik),e(ist1,ik)+scale*wt(ist1,ik)/2
    write(52,*)dpp1d(ik),e(ist1,ik)-scale*wt(ist1,ik)/2
    write(53,*)dpp1d(ik),e(ist1,ik)+scale*wt(ist1,ik)/2
    write(53,*)dpp1d(ik),e(ist1,ik)-scale*wt(ist1,ik)/2
    write(53,*)
    write(54,*)dpp1d(ik),e(ist1,ik),scale*wt(ist1,ik)
  enddo
  write(50,*)
  write(51,*)
  write(52,*)
  write(54,*)
enddo
close(50) 
close(51) 
close(52) 
close(53) 
close(54)
      
open(50,file='BNDS.GNU',form='formatted',status='replace')
write(50,*)'set term postscript portrait'
write(50,*)'set noxzeroaxis'
write(50,*)'set tics out'
write(50,*)'set noxtics'
write(50,*)'set nokey'
write(50,'(" set yrange [",F12.6,":",F12.6,"]")')emin,emax      
write(50,'(" set grid ytics")')
write(50,'(" set ytics ",F12.6)')(emax-emin)/40
write(50,'(" set xrange [",F12.6,":",F12.6,"]")')0.d0,dpp1d(nkpt)
write(50,*)"set ylabel 'Energy ("//e_units//")'"
write(50,*)"set output 'bnds.ps'"      
write(50,'("plot ''BNDS.DAT'' with lines lt 1, ''BNDS_xydy.DAT'' &
  &using 1:2:($3*4) with points lt 1 lc 1 pt 6 ps variable")')
!write(50,'("plot ''BNDS.DAT'' with line lt 1, ''BNDS1.DAT'' with line lt 1, &
!           & ''BNDS2.DAT'' with line lt 1, ''BNDS3.DAT'' with line lt 1, &
!           & ''BNDS4.DAT'' with line lt 2")')
close(50)      
allocate(lines(4,nlines))
open(50,file='BANDLINES.OUT',form='formatted',status='old')
do i=1,nlines
  read(50,*)lines(1,i),lines(2,i)
  read(50,*)lines(3,i),lines(4,i)
  read(50,*)
enddo
close(50)      
open(50,file='BNDS4.DAT',form='formatted',status='replace')
do i=1,nlines
  write(50,*)lines(1,i),emin+1.d-4
  write(50,*)lines(3,i),emax-1.d-4
  write(50,*)
enddo
write(50,*)0.d0,efermi
write(50,*)dpp1d(nkpt),efermi
close(50)
end
      