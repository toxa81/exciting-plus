! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genspecies(fnum)
implicit none
! arguments
integer, intent(in) :: fnum
! local variables
! maximum number of atomic states
integer, parameter :: maxst=40
! polynomial order for radial mesh integration
integer, parameter :: np=4
! use Perdew-Wang exchange-correlation functional
integer, parameter :: xctype=3,xcgrad=0
logical core(maxst),lorb(maxst)
integer nz,nmax,nmaxl,lmax,nst
integer ne,nrmt,nr,ir,i,j
integer nlorb,ist,jst
integer n(maxst),l(maxst),k(maxst)
real(8), parameter :: pi=3.1415926535897932385d0
! mass of 1/12 times carbon-12 in electron masses (CODATA 2006)
real(8), parameter :: amu=1822.88848426d0
! speed of light in atomic units (=1/alpha) (CODATA 2006)
real(8), parameter :: sol=137.035999679d0
! core-valence cut-off energy
real(8), parameter :: ecvcut=-3.5d0
! semi-core cut-off energy
real(8), parameter :: esccut=-0.5d0
! default APW band energy
real(8), parameter :: apwe0=0.15d0
real(8) mass,zn,t1,t2,t3
real(8) rmt,rmin,rmax
real(8) occ(maxst),eval(maxst)
logical l1(0:3,10),tlosc
character(256) symb,name
character*1 c1
character*100 str
integer, allocatable :: level(:,:)
character, parameter :: lname(4)=(/'s','p','d','f'/)
! allocatable arrays
real(8), allocatable :: r(:),rho(:),vr(:),rwf(:,:,:)
read(fnum,*,err=20) nz
if (nz.le.0) then
  write(*,*)
  write(*,'("Error(genspecies): atomic number negative : ",I8)') nz
  write(*,*)
  stop
end if
read(fnum,*,err=20) symb,name
read(fnum,*,err=20) mass
! convert from 'atomic mass units' to atomic units
mass=mass*amu
read(fnum,*,err=20) rmt
read(fnum,*,err=20) nst
if ((nst.le.0).or.(nst.gt.maxst)) then
  write(*,*)
  write(*,'("Error(genspecies): nst out of range : ",I8)') nst
  write(*,'(" for species ",A)') trim(name)
  write(*,*)
  stop
end if
ne=0
nmax=1
do ist=1,nst
  read(fnum,*,err=20) n(ist),l(ist),k(ist),i
  ne=ne+i
  occ(ist)=i
  nmax=max(nmax,n(ist))
end do
allocate(level(0:3,nmax))
level=0
do ist=1,nst
  level(l(ist),n(ist))=1
enddo
write(*,'("Info(genspecies): running Z = ",I4,", (",A,")")') nz,trim(name)
if (ne.ne.nz) then
  write(*,*)
  write(*,'("Warning(genspecies): atom not neutral, electron number : ",I4)') ne
end if
! nuclear charge in units of e
zn=-dble(nz)
! minimum radial mesh point proportional to 1/sqrt(Z)
rmin=2.d-6/sqrt(dble(nz))
! set the number of radial mesh points proportional to number of nodes
nrmt=150*(nmax+1)
! default effective infinity
rmax=100.d0
do i=1,2
! number of points to effective infinity
  t1=log(rmt/rmin)
  t2=log(rmax/rmin)
  t3=dble(nrmt)*t2/t1
  nr=int(t3)
  allocate(r(nr),rho(nr),vr(nr),rwf(nr,2,nst))
! generate logarithmic radial mesh
  t2=1.d0/dble(nrmt-1)
  do ir=1,nr
    r(ir)=rmin*exp(dble(ir-1)*t1*t2)
  end do
! solve the Kohn-Sham-Dirac equation for the atom
  call atom(sol,.true.,zn,nst,n,l,k,occ,xctype,xcgrad,np,nr,r,eval,rho,vr,rwf)
  do ir=nr,1,-1
    if (rho(ir).gt.1.d-20) then
      rmax=1.5d0*r(ir)
      goto 10
    end if
  end do
10 continue
  deallocate(r,rho,vr,rwf)
end do
! find which states belong to core
do ist=1,nst
  if (eval(ist).lt.ecvcut) then
    core(ist)=.true.
    level(l(ist),n(ist))=3
  else
    core(ist)=.false.
  end if
end do
! check that the state for same n and l but different k is also core
do ist=1,nst
  if (core(ist)) then
    do jst=1,nst
      if ((n(ist).eq.n(jst)).and.(l(ist).eq.l(jst))) core(jst)=.true.
    end do
  end if
end do
lmax=0
do ist=1,nst
  if (.not.core(ist)) lmax=max(lmax,l(ist))
end do
lmax=min(lmax+1,3)
! determine the local orbitals
nlorb=lmax+1
lorb(:)=.false.
do ist=1,nst
  if (.not.core(ist)) then
    if ((l(ist).eq.0).or.(l(ist).lt.k(ist))) then
      if (eval(ist).lt.esccut) then
        lorb(ist)=.true.
        nlorb=nlorb+1
        level(l(ist),n(ist))=2
      end if
    end if
  end if
end do
write(*,*)"levels :"
do i=1,nmax
  write(*,*)level(:,i)
enddo
! write the species file
open(60,file=trim(symb)//'.in',action='WRITE',form='FORMATTED')
write(60,'(" ''",A,"''",T45,": spsymb")') trim(symb)
write(60,'(" ''",A,"''",T45,": spname")') trim(name)
write(60,'(G14.6,T45,": spzn")') zn
write(60,'(G18.10,T45,": spmass")') mass
write(60,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') rmin,rmt,rmax, &
 nrmt
write(60,'(I4,T45,": spnst")') nst
write(60,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') n(1),l(1), &
 k(1),occ(1),core(1)
do ist=2,nst
  write(60,'(3I4,G14.6,L1)') n(ist),l(ist),k(ist),occ(ist),core(ist)
end do
write(60,'(I4,T45,": apword")') 1
write(60,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') apwe0,0,.false.
write(60,'(I4,T45,": nlx")') lmax+1
do i=0,lmax
  write(60,'(2I4,T45,": l, apword")')i,1
  write(60,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') apwe0,0,.true.
enddo
write(60,'(I4,T45,": nlorb")') nlorb
do i=0,lmax
  write(60,'(2I4,T45,": lorbl, lorbord")') i,2
  write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') apwe0,0,.true.
  write(60,'(F8.4,I4,"  ",L1)') apwe0,1,.true.
end do
do ist=1,nst
  if (lorb(ist)) then
    write(60,'(2I4,T45,": lorbl, lorbord")') l(ist),3
    write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') apwe0,0, &
     .false.
    write(60,'(F8.4,I4,"  ",L1)') apwe0,1,.false.
! subtract 0.25 Ha from eigenvalue because findband searches upwards
    write(60,'(F8.4,I4,"  ",L1)') eval(ist),0,.true.
  end if
end do
close(60)
! write to .json file
open(60,file=trim(symb)//'.json',action='WRITE',form='FORMATTED')
write(60,'("{")')
write(60,'("  ""name""   : """,A,""", ")') trim(name)
write(60,'("  ""symbol"" : """,A,""", ")') trim(symb)
write(str,'(I4)')nz
write(60,'("  ""number"" : ",A,", ")') trim(adjustl(str))
write(str,'(G18.10)')mass
write(60,'("  ""mass""   : ",A,", ")') trim(adjustl(str))
write(str,'(G18.10)')rmin
write(60,'("  ""rmin""   : ",A,", ")') trim(adjustl(str))
write(str,'(G18.10)')rmax
write(60,'("  ""rmax""   : ",A,", ")') trim(adjustl(str))
write(str,'(G18.10)')rmt
write(60,'("  ""rmt""    : ",A,", ")') trim(adjustl(str))
write(str,'(I6)')nrmt
write(60,'("  ""nrmt""   : ",A,", ")') trim(adjustl(str))
!write(60,'("  ""core""   : [")')
l1=.false.
do ist=1,nst
  if (core(ist)) then
    l1(l(ist),n(ist))=.true.
  endif
enddo
str=""
do i=1,10
  do j=0,3
    if (l1(j,i)) then
      write(c1,'(I1)')i
      str=trim(adjustl(str))//c1//lname(j+1)
    endif
  enddo
enddo
write(60,'("  ""core""   : """,A,""", ")') trim(str)
write(60,'("  ""apw"" : {")')
write(60,'("    ""default"" : [{""enu"" : 0.15, ""dme"" : 0}]")')
write(60,'("    ""explicit"" : [")')
do j=0,lmax
  write(60,'("    {")')
  write(60,'("      ""l""    : ",I1,",")')j
  nmaxl=0
  do ist=1,nst
    if (core(ist).and.j.eq.l(ist)) nmaxl=max(min(n(ist),j+1),nmaxl)
  enddo
  write(60,'("      ""conf"" : [{""n"" : ",I1,", ""enu"" : 0.15, ""dme"" : 0, ""auto"" : true}]")')nmax+1
  if (j.eq.lmax) then
    write(60,'("    }")')
  else
    write(60,'("    },")')
  endif
enddo
write(60,'("    ]")')
write(60,'("  }")')
tlosc=.false.
do ist=1,nst
  if (lorb(ist)) tlosc=.true.
enddo
write(60,'("  ""lo"" : [")')
do j=0,lmax
  write(60,'("  {")')
  write(60,'("    ""l""    : ",I1,",")')j
  nmax=0
  do ist=1,nst
    if (core(ist).and.j.eq.l(ist)) nmax=max(n(ist),nmax)
  enddo
  if (nmax.eq.0) nmax=j
  write(60,'("    ""conf"" : [")')
  write(60,'("      {""n"" : ",I1,", ""enu"" : 0.15, ""dme"" : 0, ""auto"" : true},")')nmax+1
  write(60,'("      {""n"" : ",I1,", ""enu"" : 0.15, ""dme"" : 1, ""auto"" : true}")')nmax+1
  write(60,'("    ]")')
  if (j.eq.lmax.and..not.tlosc) then
    write(60,'("  }")')
  else
    write(60,'("  },")')
  endif
enddo
do ist=1,nst
  if (lorb(ist)) then
    write(60,'("  {")')
    write(60,'("    ""l""    : ",I1,",")')j
    nmax=0
    !do ist=1,nst
    !  if (core(ist).and.j.eq.l(ist)) nmax=max(n(ist),nmax)
    !enddo
    if (nmax.eq.0) nmax=j
    write(60,'("    ""conf"" : [")')
    write(60,'("      {""n"" : ",I1,", ""enu"" : 0.15, ""dme"" : 0, ""auto"" : true},")')nmax+1
    write(60,'("      {""n"" : ",I1,", ""enu"" : 0.15, ""dme"" : 1, ""auto"" : true}")')nmax+1
    write(60,'("    ]")')
 
     write(60,'(2I4,T45,": lorbl, lorbord")') l(ist),3
    write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') apwe0,0, &
     .false.
    write(60,'(F8.4,I4,"  ",L1)') apwe0,1,.false.
! subtract 0.25 Ha from eigenvalue because findband searches upwards
    write(60,'(F8.4,I4,"  ",L1)') eval(ist),0,.true.
  end if
end do
write(60,'("  ]")')
!i=0
!do ist=1,nst
!  if (core(ist)) i=ist
!enddo
!do ist=1,nst
!  if (core(ist)) then
!    write(60,'("    {")')
!    write(60,'("      ""n"" : ",I1,",")')n(ist)
!    write(60,'("      ""l"" : ",I1,",")')l(ist)
!    write(60,'("      ""k"" : ",I1,",")')k(ist)
!    write(60,'("      ""occupancy"" : ",I1)')iocc(ist)
!    if (i.eq.ist) then
!      write(60,'("    }")')
!    else
!      write(60,'("    },")')
!    endif
!  endif
!enddo
!write(60,'("  ]")')

!write(60,'(G14.6,2F10.4,I6,T45,": sprmin, rmt, sprmax, nrmt")') rmin,rmt,rmax, &
! nrmt
!write(60,'(I4,T45,": spnst")') nst
!write(60,'(3I4,G14.6,L1,T45,": spn, spl, spk, spocc, spcore")') n(1),l(1), &
! k(1),occ(1),core(1)
!do ist=2,nst
!  write(60,'(3I4,G14.6,L1)') n(ist),l(ist),k(ist),occ(ist),core(ist)
!end do
!write(60,'(I4,T45,": apword")') 1
!write(60,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') apwe0,0,.false.
!write(60,'(I4,T45,": nlx")') lmax+1
!do i=0,lmax
!  write(60,'(2I4,T45,": l, apword")')i,1
!  write(60,'(F8.4,I4,"  ",L1,T45,": apwe0, apwdm, apwve")') apwe0,0,.true.
!enddo
!write(60,'(I4,T45,": nlorb")') nlorb
!do i=0,lmax
!  write(60,'(2I4,T45,": lorbl, lorbord")') i,2
!  write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') apwe0,0,.true.
!  write(60,'(F8.4,I4,"  ",L1)') apwe0,1,.true.
!end do
!do ist=1,nst
!  if (lorb(ist)) then
!    write(60,'(2I4,T45,": lorbl, lorbord")') l(ist),3
!    write(60,'(F8.4,I4,"  ",L1,T45,": lorbe0, lorbdm, lorbve")') apwe0,0, &
!     .false.
!    write(60,'(F8.4,I4,"  ",L1)') apwe0,1,.false.
!! subtract 0.25 Ha from eigenvalue because findband searches upwards
!    write(60,'(F8.4,I4,"  ",L1)') eval(ist),0,.true.
!  end if
!end do
close(60)
deallocate(level)
return
20 continue
write(*,*)
write(*,'("Error(genspecies): error reading species data")')
write(*,*)
stop
end subroutine

