subroutine writegeom
use modmain
implicit none
! local variables
integer is,ia,ip
open(50,file='GEOMETRY.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("! Atomic positions generated by spacegroup version ",&
 &I1.1,".",I1.1,".",I2.2)') version
write(50,'("!  Hermann-Mauguin symbol : ",A)') trim(hrmg)
write(50,'("!  Hall symbol            : ",A)') trim(hall)
write(50,'("!  Schoenflies symbol     : ",A)') trim(schn)
write(50,'("!  space group number     : ",A)') trim(num)
write(50,'("!  lattice constants (a,b,c) : ",3G18.10)') a,b,c
write(50,'("!  angles in degrees (ab,ac,bc) : ",3G18.10)') ab,ac,bc
write(50,'("!  number of conventional unit cells : ",3I4)') ncell
write(50,'("!  reduction to primitive cell : ",L1)') primcell
write(50,'("!  Wyckoff positions :")')
do is=1,nspecies
  write(50,'("!   species : ",I4,", ",A)') is,trim(spfname(is))
  do ip=1,nwpos(is)
    write(50,'("!   ",3G18.10)') wpos(:,ip,is)
  end do
end do
write(50,*)
write(50,'("avec")')
write(50,'(3G18.10)') avec(:,1)
write(50,'(3G18.10)') avec(:,2)
write(50,'(3G18.10)') avec(:,3)
write(50,*)
write(50,'("atoms")')
write(50,'(I4,T40," : nspecies")') nspecies
do is=1,nspecies
  write(50,'("''",A,"''",T40," : spfname")') trim(spfname(is))
  write(50,'(I4,T40," : natoms; atposl, bfcmt below")') natoms(is)
  do ia=1,natoms(is)
    write(50,'(3F14.8,"  ",3F12.8)') atposl(:,ia,is),0.d0,0.d0,0.d0
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writegeom):")')
write(*,'(" EXCITING lattice vectors and atomic positions written to &
 &GEOMETRY.OUT")')
return
end subroutine
