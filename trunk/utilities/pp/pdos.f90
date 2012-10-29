program main
implicit none

integer nspecies,natmtot,nspinor,lmax,lmmax,nwdos,n,nwann
character(256), allocatable :: spsymb(:)
integer, allocatable :: natoms(:)
real(8), allocatable :: w(:),pdos(:,:,:,:),doswan(:,:)
integer, allocatable :: iasis(:)
integer iw,is,ias,ispn,lm,ias1,is1
integer, allocatable :: orb(:,:,:)
character(20) :: strin
real(8), allocatable :: dos(:)
logical wannier,l1
integer iwann(100)

open(160,file='PDOS.OUT',form='UNFORMATTED',status='OLD')
read(160)nspecies,natmtot,nspinor,lmax,lmmax,nwdos
allocate(spsymb(nspecies))
allocate(natoms(nspecies))
allocate(iasis(natmtot))
allocate(w(nwdos))
allocate(pdos(nwdos,lmmax,nspinor,natmtot))
allocate(dos(nwdos))
dos=0.d0
do is=1,nspecies
  read(160)spsymb(is)
  read(160)natoms(is)
enddo
do ias=1,natmtot
  read(160)ias1,is1
  iasis(ias1)=is1
enddo
do iw=1,nwdos
  read(160)w(iw)
enddo
do ias=1,natmtot
  do ispn=1,nspinor
    do lm=1,lmmax
      do iw=1,nwdos
        read(160)pdos(iw,lm,ispn,ias)
      end do !iw
    end do !lm
  end do !ispn
end do !ias
read(160)wannier
if (wannier) then
  read(160)nwann
  allocate(doswan(nwdos,nwann))
  do n=1,nwann
    do iw=1,nwdos
      read(160)doswan(iw,n)
    enddo
  enddo
endif
close(160)

l1=.false.
if (wannier) then
  write(*,'("Wannier DOS? (T/F)")')
  read(*,*)l1
  if (l1) goto 40
endif

! regular pDOS
write(*,'("Number of species : ",I4)')nspecies
write(*,'("Number of atoms : ",I4)')natmtot
write(*,'("Number of spins : ",I4)')nspinor
allocate(orb(lmmax,nspinor,natmtot))
orb=0
10 continue
write(*,'("Input orbitals for pDOS")')
read(*,'(A)')strin
if (trim(strin).eq.'q') goto 20
call addorb(strin,orb,lmmax,nspinor,natmtot,iasis,spsymb)
goto 10
20 continue
do ispn=1,nspinor
  write(*,'("spin : ",I1)')ispn
  do ias=1,natmtot
    write(*,'("  atom : ",I4," orb : ",16I2)')ias,(orb(lm,ispn,ias),lm=1,lmmax)
    do lm=1,lmmax
      dos(:)=dos(:)+orb(lm,ispn,ias)*pdos(:,lm,ispn,ias)
    enddo
  enddo
enddo
goto 50
! Wannier pDOS
40 continue
iwann=-1
write(*,'("Input Wannier functions")')
read(*,*)iwann
do n=1,100
  if (iwann(n).ne.-1) dos(:)=dos(:)+doswan(:,iwann(n))
enddo
50 continue
open(50,file='dos.dat',form='formatted',status='replace')
do iw=1,nwdos
  write(50,'(2G18.10)')w(iw),dos(iw)
enddo
close(50)

end
