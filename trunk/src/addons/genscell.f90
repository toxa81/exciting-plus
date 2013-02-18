subroutine genscell
use modmain
implicit none
real(8) sca(3,3),scainv(3,3)
integer i,i1,i2,i3,ia,is
real(8) v1(3),v2(3)
real(8), allocatable :: vecl(:,:)
integer iv(3),nsc,n,j
real(8), external :: r3mdet
logical l1

call init0
!call init1
call mpi_grid_initialize((/1/))

do i=1,3
  sca(:,i)=scvl(1,i)*avec(:,1)+scvl(2,i)*avec(:,2)+scvl(3,i)*avec(:,3)
enddo
nsc=int(abs(r3mdet(sca)/r3mdet(avec))+epslat)
write(*,*)
write(*,'("Lattice vectors")')
write(*,'(" a1 : ",3G18.10)')avec(:,1)
write(*,'(" a2 : ",3G18.10)')avec(:,2)
write(*,'(" a3 : ",3G18.10)')avec(:,3)
write(*,*)
write(*,'("Supercell lattice vectors")')
write(*,'(" a1 : ",3G18.10)')sca(:,1)
write(*,'(" a2 : ",3G18.10)')sca(:,2)
write(*,'(" a3 : ",3G18.10)')sca(:,3)
write(*,*)
write(*,'("Volume ratio : ",I3)')nsc

open(50,file='GEOMETRY_SCELL'//trim(filext),action='WRITE',form='FORMATTED')
write(50,'("avec")')
write(50,'(3G18.10)') sca(:,1)
write(50,'(3G18.10)') sca(:,2)
write(50,'(3G18.10)') sca(:,3)
write(50,*)
write(50,'("atoms")')
write(50,'(I4,T40," : nspecies")') nspecies

call r3minv(sca,scainv)

do is=1,nspecies
  write(50,'(" ''",A,"''",T40," : spfname")') trim(spfname(is))
  write(50,'(I4,T40," : natoms; atpos, bfcmt below")') natoms(is)*nsc
  allocate(vecl(3,natoms(is)*nsc))
  n=0
  do ia=1,natoms(is)
    do i1=-10,10
    do i2=-10,10
    do i3=-10,10
      v1(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      call r3mv(scainv,v1,v2)
      call r3frac(epslat,v2,iv)
      l1=.true.
      do j=1,n
        if (abs(v2(1)-vecl(1,j)).lt.epslat.and.&
            abs(v2(2)-vecl(2,j)).lt.epslat.and.&
            abs(v2(3)-vecl(3,j)).lt.epslat) l1=.false.
      enddo
      if (l1) then
        n=n+1
        vecl(:,n)=v2(:)
      endif
    enddo
    enddo
    enddo
  enddo
  do j=1,n
    write(50,'(3F14.8,"  ",3F12.8)') vecl(:,j),0.d0,0.d0,0.d0
  enddo  
  deallocate(vecl)
enddo
close(50)

open(50,file='input.xml',action='WRITE',form='FORMATTED')
write(50,'("<structure>")')
write(50,'("  <crystal>")')
write(50,'("    <basevect>",3G18.10,"</basevect>")') sca(:,1)
write(50,'("    <basevect>",3G18.10,"</basevect>")') sca(:,2)
write(50,'("    <basevect>",3G18.10,"</basevect>")') sca(:,3)
write(50,'("  </crystal>")')

call r3minv(sca,scainv)

do is=1,nspecies
  write(50,'("  <species speciesfile=""",A,""">")')trim(spfname(is))
  allocate(vecl(3,natoms(is)*nsc))
  n=0
  do ia=1,natoms(is)
    do i1=-2,2
    do i2=-2,2
    do i3=-2,2
      v1(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      call r3mv(scainv,v1,v2)
      call r3frac(epslat,v2,iv)
      l1=.true.
      do j=1,n
        if (abs(v2(1)-vecl(1,j)).lt.epslat.and.&
            abs(v2(2)-vecl(2,j)).lt.epslat.and.&
            abs(v2(3)-vecl(3,j)).lt.epslat) l1=.false.
      enddo
      if (l1) then
        n=n+1
        vecl(:,n)=v2(:)
      endif
    enddo
    enddo
    enddo
  enddo
  do j=1,n
    write(50,'("    <atom coord=""",3F14.8,"""/>")')vecl(:,j)
  enddo  
  deallocate(vecl)
  write(50,'("  </species>")')
enddo
write(50,'("</structure>")')
close(50)

return
end