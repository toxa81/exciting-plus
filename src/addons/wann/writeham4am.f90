subroutine writeham4am
use modmain
use modldapu
use mod_nrkp
implicit none

integer values(8), seconds, hash, isp, hdim, ntype, itype, istart, l, norb, lnext
integer is, ias, ia, nkp
real(8)  reh, imh, sc, sc1, sc2, sc3, alat(3,3)
character*1 xx
character(256) block
integer, allocatable   :: basis_desc(:,:)
! basis_desc(12, ntype)
! ntype            : number of WF types
! basis_desc(1, )  : number of atom, functions of that are projected to itype
! basis_desc(2, )  : orbital quantum number of basis functions
! basis_desc(3, )  : number of the first orbital of this type in hamiltonian
! basis_desc(4, )  : amount of orbitals in this group
! basis_desc(5, )  : number of this orbital 1 for s- states,
!    ...           : (2,3,4) for p- states, (5, 6, 7, 8, 9) for d- states,
! basis_desc(12, ) : (10, 11, 12, 13, 14, 15, 16)  for f - states
integer i, j, ik
real*8,            parameter :: zero=0.d0

if (task.eq.829) then
  nkp = nkpt
else
  nkp = nkptnr
end if


call date_and_time(values=values)
! Inaccurate but easy calculation of timestamp.
hash = (values(1) - 1970) * 365 * 24 * 60 * 60  &
   + (values(2)-1) * 31 * 24 * 60 * 60  &
   + (values(3) + 11 ) * 24 * 60 * 60  &
   + (values(5)-values(4)) * 60 * 60  &
   + values(6) * 60  + values(7)  &
   + int(values(8)/1000)

open(200, file="hamilt.am", form="FORMATTED", status="REPLACE")

write(200,'("# This file was written on: ",I2.2,".",I2.2,".",I4.4,"  &
                                         ",I2.2,":",I2.2,":",I2.2)') &
           values(3),values(2),values(1),values(5),values(6),values(7)
write(200,*)

write(200,'(a5)') '&hash'
write(200,*) hash
write(200,*)

write(200,'(a10)') '&Codestamp'
write(200,*) 'exciting-plus'
write(200,*)


write(200,'(a6)') '&nspin'
write(200,*) nspinor
write(200,*)

if (task.eq.829) then
  write(200,'("# Hamiltonian is written along high-symmetry directions")')
endif
write(200,'(a4)') '&nkp'
write(200,*) nkp
write(200,*)

hdim = int(nwantot/nspinor)
write(200,'(a4)') '&dim'
write(200,*) hdim !nwantot
write(200,*)

open(50,file='TOTENERGY.OUT',action='READ',form='FORMATTED',status='OLD')
do while (.true.)
  read(50, '(G22.12)', iostat=i) engytot
  if (i /= 0) exit
end do
close(50)
write(200,'(a5)') '&etot'
write(200,*) engytot * ha2ev
write(200,*)

write(200,'(a6)') '&fermi'
write(200,*) efermi * ha2ev
write(200,*)


write(200,'(a11)') '&crystcoord'
write(200,*) 'true'
write(200,*)

write(200,'(a8)') '&kpoints'
if (task.eq.829)then
  do ik=1,nkp ! tnr
    write(200,'(f15.12,6f9.5)') 1.d0, vkl(:,ik)
  end do
else
  do ik=1,nkp ! tnr
    write(200,'(f15.12,6f9.5)') wkptnr(ik), vklnr(:,ik)
  end do
endif
write(200,*)

hdim = int(nwantot/nspinor)
write(200,'(a12)') '&hamiltonian'
do isp=1,nspinor
do ik = 1, nkp
  do i = 1, nwantot
  do j = i, nwantot
    if (wan_info(wi_spin,i) /= isp .or. wan_info(wi_spin,j) /= isp) cycle
    reh=dreal(wann_h(i,j,ik)) * ha2ev
    imh=dimag(wann_h(i,j,ik)) * ha2ev
    if(dabs(reh)<1.d-12 .and. dabs(imh)<1.d-12) then
    write(200,'(6x,f2.0,18x,f2.0)')zero,zero
    elseif(dabs(reh)<1.d-12) then
    write(200,'(6x,f2.0,12x,f20.12)')zero,imh
    elseif(dabs(imh)<1.d-12) then
    write(200,'(f20.12,6x,f2.0)')reh,zero
    else
    write(200,'(2f20.12)')reh,imh
    end if
  end do !n2
  end do   !n1
end do     !ikp
end do       !isp

close(200)

open(200, file="system.am", form="FORMATTED", status="REPLACE")

write(200,'("# This file was written on: ",I2.2,".",I2.2,".",I4.4,"  &
                                         ",I2.2,":",I2.2,":",I2.2)') &
           values(3),values(2),values(1),values(5),values(6),values(7)
write(200,*)

write(200,'(a5)') '&hash'
write(200,*) hash
write(200,*)

write(200,'(a10)') '&Codestamp'
write(200,*) 'exciting-plus'
write(200,*)

!print*,alat
sc  = 1.0
sc1 = 1.0
sc2 = 1.0
sc3 = 1.0
open(50,file='elk.in',action='READ',status='OLD',form='FORMATTED')
10 continue
read(50,*,end=30) block
if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
select case(trim(block))
case('avec')
  read(50,*,err=20) alat(:,1)
  read(50,*,err=20) alat(:,2)
  read(50,*,err=20) alat(:,3)
case('scale')
  read(50,*,err=20) sc
case('scale1')
  read(50,*,err=20) sc1
case('scale2')
  read(50,*,err=20) sc2
case('scale3')
  read(50,*,err=20) sc3
case('')
  goto 10
end select
goto 10
20 continue
write(*,*)
write(*,'("Error(readinput): error reading from elk.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
write(*,'("Check input convention in manual")')
write(*,*)
stop
30 continue
close(50)

alat(:,1)=sc1*alat(:,1)
alat(:,2)=sc2*alat(:,2)
alat(:,3)=sc3*alat(:,3)

write(200,'(a5)') '&cell'
write(200,'(f12.9)') sc !1.0
do i=1,3
  write(200,'(3f9.5)') alat(:,i) !avec(:,i) !/tmp
end do
write(200,*)

write(200,'(a6)') '&fermi'
write(200,*) efermi * ha2ev
write(200,*)

write(200,'(a11)') '&crystcoord'
write(200,*) 'true'
write(200,*)

write(200,'(a6)') '&atoms'
write(200,*) natmtot
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  write(200,'(a4,x,3f9.5)') trim(spsymb(is)), atposl(:,ia,is)
end do
write(200,*)

ntype = int(wann_natom/nspinor)
allocate (basis_desc(12,ntype))
basis_desc = 0
basis_desc(4,:) = 1
itype = 1
istart = 1
j = 0
do i = 1,hdim-1
  basis_desc(1,itype) = wan_info(wi_atom,i)
  l = wan_info(wi_lm,i)
  if (l < 2)                 basis_desc(2,itype) = 0
  if (l < 5 .and. l >= 2 )   basis_desc(2,itype) = 1
  if (l < 10 .and. l >= 5 )  basis_desc(2,itype) = 2
  if (l < 17 .and. l >= 10 ) basis_desc(2,itype) = 3
  basis_desc(3,itype) = istart
  l = basis_desc(2,itype)
  if (basis_desc(2,itype) == 0) norb = 1
  if (basis_desc(2,itype) == 1) norb = 3
  if (basis_desc(2,itype) == 2) norb = 5
  if (basis_desc(2,itype) == 3) norb = 7

  basis_desc(5+j,itype) = wan_info(wi_lm,i)

  lnext = wan_info(wi_lm,i+1)
  if (lnext < 2)                     lnext = 0
  if (lnext < 5 .and.  lnext >= 2 )  lnext = 1
  if (lnext < 10 .and. lnext >= 5 )  lnext = 2
  if (lnext < 17 .and. lnext >= 10 ) lnext = 3

  if ((wan_info(wi_atom,i) .ne. wan_info(wi_atom,i+1)) .or. &
        l .ne. lnext) then
    itype = itype +1
    istart = istart + norb
    j = 0
  else
    basis_desc(4, itype) = basis_desc(4,itype) + 1
    basis_desc(5+j+1,itype) = wan_info(wi_lm,i+1) !
    j = j + 1
  end if
end do

write(200,'(a20)') '# Basis description:'
write(200,'(a14)') '# dim, nblocks'
write(200,'(a74)') '# atom_sym, atom_num, l_sym, block_dim, block_start, orbitals(1:block_dim)'
write(200,'(a6)') '&basis'
write(200,'(i2,i4)') hdim, ntype
do i = 1,ntype
  select case(basis_desc(2,i))
  case(0)
    xx='s'
    norb = 1
  case(1)
    xx='p'
    norb = 3
  case(2)
    xx='d'
    norb = 5
  case(3)
    xx='f'
    norb = 7
  end select
  norb = basis_desc(4, i)
  is=ias2is(basis_desc(1,i))
  write(200,'(a3,i3,a2,i2,i4,4x,16i2)')trim(spsymb(is)), basis_desc(1,i), &
             xx, basis_desc(4, i), basis_desc(3,i), (basis_desc(5+j,i), j=0,norb-1)
end do

write(200,*)

deallocate(basis_desc)
close(200)

return
end subroutine writeham4am
