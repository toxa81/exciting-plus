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

! TODO:
!  Input in the A,O,S (atom, orbital, spin) format
!  A is a number or a string
!  O is a number or a string
!  S is a number
!  ,,    : all states
!  ,,1   : spin up
!  ,,2   : spin down
!  ,d,1  : all d's of all atoms spin 1
!  ,2,   : y states of all atoms of all spins
!  2,,   : all states of atom No.2
!  Fe,d, : all d states of all Fe atoms

orb_l=0
orb_s=0
orb_a=0

! split input string to 2 or 3 substrings
i1=scan(strin,",",.false.)
if (i1.eq.0) then
  write(*,*)
  write(*,'("Warning(addorb): wrong string")')
  write(*,'(A)')trim(strin)
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
read(str1,'(I4)',err=20)ias
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
read(str2,'(I4)',err=30)lm
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
  read(str3,'(I4)',err=40)ispn
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