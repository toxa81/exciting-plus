program main
implicit none

integer nspecies,natmtot,nspinor,lmax,lmmax,nwdos
character(256), allocatable :: spsymb(:)
integer, allocatable :: natoms(:)
real(8), allocatable :: w(:),pdos(:,:,:,:)
integer, allocatable :: iasis(:)
integer iw,is,ias,ispn,lm,ias1,is1
integer, allocatable :: orb(:,:,:)
character(20) :: strin
real(8), allocatable :: dos(:)

open(160,file='PDOS.OUT',form='UNFORMATTED',status='OLD')
read(160)nspecies,natmtot,nspinor,lmax,lmmax,nwdos
allocate(spsymb(nspecies))
allocate(natoms(nspecies))
allocate(iasis(natmtot))
allocate(w(nwdos))
allocate(pdos(nwdos,lmmax,nspinor,natmtot))
allocate(dos(nwdos))
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
close(160)

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
dos=0.d0
do ispn=1,nspinor
  write(*,'("spin : ",I1)')ispn
  do ias=1,natmtot
    write(*,'("  atom : ",I4," orb : ",16I2)')ias,(orb(lm,ispn,ias),lm=1,lmmax)
    do lm=1,lmmax
      dos(:)=dos(:)+orb(lm,ispn,ias)*pdos(:,lm,ispn,ias)
    enddo
  enddo
enddo
open(50,file='dos.dat',form='formatted',status='replace')
do iw=1,nwdos
  write(50,'(2G18.10)')w(iw),dos(iw)
enddo
close(50)

return
end

subroutine addorb(strin,orb,lmmax,nspinor,natmtot,iasis,spsymb)
implicit none
character(20), intent(in) :: strin
integer, intent(inout) :: orb(lmmax,nspinor,natmtot)
integer, intent(in) :: lmmax
integer, intent(in) :: nspinor
integer, intent(in) :: natmtot
integer, intent(in) :: iasis(natmtot)
character(256), intent(in) :: spsymb(*)
character(20) str1,str2,str3
integer i1,i2,ias,lm,ispn
integer orb_l(lmmax)
integer orb_s(nspinor)
integer orb_a(natmtot)

orb_l=0
orb_s=0
orb_a=0

! split input string to 2 or 3 substrings
i1=scan(strin,",",.false.)
if (i1.eq.0) then
  write(*,10)
  return
endif
str1=strin(1:i1-1)
i2=scan(strin(i1+1:),",",.false.)
if (i2.ne.0) then
  str2=strin(i1+1:i1+i2-1)
  str3=strin(i1+i2+1:)
else
  str2=strin(i1+1:)
  str3=''
endif

! try to read atom number
read(str1,'(I4)',err=20),ias
if (ias.gt.natmtot.or.ias.lt.1) then
  write(*,'("wrong atom number")')
  return
endif
orb_a(ias)=1
goto 21
20 continue
do ias=1,natmtot
  if (trim(str1).eq.trim(spsymb(iasis(ias)))) orb_a(ias)=1
enddo
21 continue

!try to read orbital
read(str2,'(I4)',err=30),lm
if (lm.gt.lmmax.or.lm.lt.1) then
  write(*,'("wrong orbital")')
  return
endif
orb_l(lm)=1
goto 31
30 continue
if (trim(str2).eq.'s') orb_l(1)=1
if (trim(str2).eq.'p') orb_l(2:4)=1
if (trim(str2).eq.'d') orb_l(5:9)=1
if (trim(str2).eq.'f') orb_l(10:16)=1
31 continue

!try to read spin
if (trim(str3).eq.'') then
  orb_s=1
else
  read(str3,'(I4)',err=40),ispn
  if (ispn.ge.1.and.ispn.le.nspinor) then
    orb_s(ispn)=1
    goto 41
  endif
40 continue
  write(*,'("wrong spin")')
  return
endif
41 continue

do ias=1,natmtot
  do lm=1,lmmax
    do ispn=1,nspinor
      if (orb(lm,ispn,ias).eq.0) orb(lm,ispn,ias)=orb_a(ias)*orb_l(lm)*orb_s(ispn)
    enddo
  enddo
enddo

10 format('wrong string')
return
end