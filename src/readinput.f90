
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
subroutine readinput
! !USES:
use modmain
use modldapu
use modrdm
use modtest
use modqpt
use mod_addons_q
use mod_sic
use mod_wannier
use mod_linresp
use mod_expigqr
! !DESCRIPTION:
!   Reads in the input parameters from the file {\tt elk.in}. Also sets default
!   values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,js,ia,ias,j,lm1,lm2
integer i,l,k,iv,iostat,n
real(8) sc,sc1,sc2,sc3
real(8) solscf,v(3)
character(256) block,str,str1
character*100 fname

!call getarg(1,fname)
!if (trim(adjustl(fname)).eq."") fname="elk.in"
fname="elk.in"

!------------------------!
!     default values     !
!------------------------!
ntasks=1
tasks(1)=-1
avec(:,:)=0.d0
avec(1,1)=1.d0
avec(2,2)=1.d0
avec(3,3)=1.d0
sc=1.d0
sc1=1.d0
sc2=1.d0
sc3=1.d0
epslat=1.d-6
primcell=.false.
tshift=.true.
ngridk(:)=1
vkloff(:)=0.d0
autokpt=.false.
radkpt=40.0
reducek=1
ngridq(:)=1
reduceq=1
rgkmax=7.d0
gmaxvr=12.d0
lmaxapw=8
lmaxvr=7
lmaxmat=5
lmaxinr=7
fracinr=0.d0
npsden=9
xctype(1)=3
xctype(2)=0
xctype(3)=0
stype=0
swidth=0.01d0
autoswidth=.false.
mstar=10.d0
epsocc=1.d-8
epschg=1.d-3
nempty=6
maxscl=200
mixtype=1
beta0=0.05d0
betamax=1.d0
epspot=1.d-6
epsengy=1.d-4
epsforce=5.d-4
molecule=.false.
nspecies=0
natoms(:)=0
atposl(:,:,:)=0.d0
atposc(:,:,:)=0.d0
bfcmt0(:,:,:)=0.d0
sppath='./'
scrpath='./'
nvp1d=2
if (allocated(vvlp1d)) deallocate(vvlp1d)
allocate(vvlp1d(3,nvp1d))
vvlp1d(:,1)=0.d0
vvlp1d(:,2)=1.d0
npp1d=200
vclp2d(:,:)=0.d0
vclp2d(1,2)=1.d0
vclp2d(2,3)=1.d0
np2d(:)=40
vclp3d(:,:)=0.d0
vclp3d(1,2)=1.d0
vclp3d(2,3)=1.d0
vclp3d(3,4)=1.d0
np3d(:)=20
nwdos=500
ngrdos=100
nsmdos=3
wdos(1)=-0.5d0
wdos(2)=0.5d0
dosocc=.false.
dosmsum=.false.
dosssum=.false.
lmirep=.true.
spinpol=.false.
spinorb=.false.
tau0atm=0.2d0
nstfsp=6
lradstp=4
chgexs=0.d0
nprad=4
scissor=0.d0
noptcomp=1
optcomp(:,1)=1
usegdft=.false.
intraband=.false.
evaltol=-1.d0
deband=0.0025d0
epsband=1.d-6
autolinengy=.false.
dlefe=-0.1d0
bfieldc0(:)=0.d0
efieldc(:)=0.d0
afieldc(:)=0.d0
fixspin=0
momfix(:)=0.d0
mommtfix(:,:,:)=0.d0
taufsm=0.01d0
autormt=.false.
rmtapm(1)=0.25d0
rmtapm(2)=0.95d0
isgkmax=-1
nosym=.false.
deltaph=0.03d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3,nphwrt))
vqlwrt(:,:)=0.d0
notelns=0
tforce=.false.
tfibs=.true.
maxitoep=120
tauoep(1)=1.d0
tauoep(2)=0.75d0
tauoep(3)=1.25d0
nkstlist=1
kstlist(:,1)=1
vklem(:)=0.d0
deltaem=0.025d0
ndspem=1
nosource=.false.
spinsprl=.false.
ssdph=.true.
vqlss(:)=0.d0
nwrite=2
tevecsv=.false.

! LDA+U defaults
ldapu=0
inptypelu=1
llu(:)=-1
ujlu(:,:)=0.d0
flu(:,:)=0.d0
elu(:,:)=0.d0
lambdalu(:)=0.d0
ulufix(:)=0.d0
lambdalu0(:)=0.d0
tmomlu=.false.
readalu=.false.

! reduced density matrix functional theory (RMDFT) defaults
rdmxctype=2
rdmmaxscl=1
maxitn=250
maxitc=10
taurdmn=1.d0
taurdmc=0.5d0
rdmalpha=0.7d0
rdmbeta=0.25d0
rdmtemp=0.d0
wrtvnlijji=.true.

reducebf=1.d0
ptnucl=.true.
tseqit=.false.
nseqit=40
tauseq=0.1d0
vecql(:)=0.d0
mustar=0.15d0
sqados(1:2)=0.d0
sqados(3)=1.d0
test=.false.
frozencr=.false.
spincore=.false.
solscf=1.d0
emaxelnes=-1.2d0

!--------------------------!
!     read from elk.in     !
!--------------------------!
open(50,file=trim(adjustl(fname)),action='READ',status='OLD',form='FORMATTED',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readinput): error opening elk.in")')
  write(*,*)
  stop
end if
10 continue
read(50,*,end=30) block
! check for a comment
if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
select case(trim(block))
case('tasks')
  do i=1,maxtasks
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): no tasks to perform")')
        write(*,*)
        stop
      end if
      ntasks=i-1
      goto 10
    end if
    read(str,*,iostat=iostat) tasks(i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading tasks")')
      write(*,'("(blank line required after tasks block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many tasks")')
  write(*,*)
  stop
case('species')
  call genspecies(50)
case('avec')
  read(50,*,err=20) avec(:,1)
  read(50,*,err=20) avec(:,2)
  read(50,*,err=20) avec(:,3)
case('scale')
  read(50,*,err=20) sc
case('scale1')
  read(50,*,err=20) sc1
case('scale2')
  read(50,*,err=20) sc2
case('scale3')
  read(50,*,err=20) sc3
case('epslat')
  read(50,*,err=20) epslat
  if (epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epslat <= 0 : ",G18.10)') epslat
    write(*,*)
    stop
  end if
case('primcell')
  read(50,*,err=20) primcell
case('tshift')
  read(50,*,err=20) tshift
case('autokpt')
  read(50,*,err=20) autokpt
case('radkpt')
  read(50,*,err=20) radkpt
  if (radkpt.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): radkpt <= 0 : ",G18.10)') radkpt
    write(*,*)
    stop
  end if
case('ngridk')
  read(50,*,err=20) ngridk(:)
  if ((ngridk(1).le.0).or.(ngridk(2).le.0).or.(ngridk(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
    write(*,*)
    stop
  end if
case('vkloff')
  read(50,*,err=20) vkloff(:)
case('reducek')
  read(50,*,err=20) reducek
case('ngridq')
  read(50,*,err=20) ngridq(:)
  if ((ngridq(1).le.0).or.(ngridq(2).le.0).or.(ngridq(3).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
    write(*,*)
    stop
  end if
case('reduceq')
  read(50,*,err=20) reduceq
case('rgkmax')
  read(50,*,err=20) rgkmax
  if (rgkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
    write(*,*)
    stop
  end if
case('gmaxvr')
  read(50,*,err=20) gmaxvr
case('lmaxapw')
  read(50,*,err=20) lmaxapw
  if (lmaxapw.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw < 0 : ",I8)') lmaxapw
    write(*,*)
    stop
  end if
  if (lmaxapw.ge.maxlapw) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw too large : ",I8)') lmaxapw
    write(*,'("Adjust maxlapw in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('lmaxvr')
  read(50,*,err=20) lmaxvr
  if (lmaxvr.lt.3) then
    write(*,*)
    write(*,'("Error(readinput): lmaxvr < 3 : ",I8)') lmaxvr
    write(*,*)
    stop
  end if
case('lmaxmat')
  read(50,*,err=20) lmaxmat
  if (lmaxmat.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxmat < 0 : ",I8)') lmaxmat
    write(*,*)
    stop
  end if
case('lmaxinr')
  read(50,*,err=20) lmaxinr
  if (lmaxinr.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxinr < 0 : ",I8)') lmaxinr
    write(*,*)
    stop
  end if
case('fracinr')
  read(50,*,err=20) fracinr
case('npsden')
  read(50,*,err=20) npsden
  if (npsden.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): npsden < 2 : ",I8)') npsden
    write(*,*)
    stop
  end if
case('spinpol')
  read(50,*,err=20) spinpol
case('spinorb')
  read(50,*,err=20) spinorb
case('xctype')
  read(50,'(A256)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) xctype
case('stype')
  read(50,*,err=20) stype
case('swidth')
  read(50,*,err=20) swidth
  if (swidth.lt.1.d-9) then
    write(*,*)
    write(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
     swidth
    write(*,*)
    stop
  end if
case('autoswidth')
  read(50,*,err=20) autoswidth
case('mstar')
  read(50,*,err=20) mstar
  if (mstar.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): mstar <= 0 : ",G18.10)') mstar
    write(*,*)
    stop
  end if
case('epsocc')
  read(50,*,err=20) epsocc
  if (epsocc.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
    write(*,*)
    stop
  end if
case('epschg')
  read(50,*,err=20) epschg
  if (epschg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
    write(*,*)
    stop
  end if
case('nempty')
  read(50,*,err=20) nempty
  if (nempty.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nempty <= 0 : ",I8)') nempty
    write(*,*)
    stop
  end if
case('mixtype')
  read(50,*,err=20) mixtype
case('beta0')
  read(50,*,err=20) beta0
  if (beta0.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): beta0 < 0 : ",G18.10)') beta0
    write(*,*)
    stop
  end if
case('betamax')
  read(50,*,err=20) betamax
  if ((betamax.lt.0.d0).or.(betamax.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): betmax not in [0,1] : ",G18.10)') betamax
    write(*,*)
    stop
  end if
case('maxscl')
  read(50,*,err=20) maxscl
  if (maxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    write(*,*)
    stop
  end if
case('epspot')
  read(50,*,err=20) epspot
case('epsengy')
  read(50,*,err=20) epsengy
case('epsforce')
  read(50,*,err=20) epsforce
case('sppath')
  read(50,*,err=20) sppath
  sppath=adjustl(sppath)
case('scrpath')
  read(50,*,err=20) scrpath
case('molecule')
  read(50,*,err=20) molecule
case('atoms')
  read(50,*,err=20) nspecies
  if (nspecies.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
    write(*,*)
    stop
  end if
  if (nspecies.gt.maxspecies) then
    write(*,*)
    write(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
    write(*,'("Adjust maxspecies in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do is=1,nspecies
    read(50,*,err=20) spfname(is)
    spfname(is)=adjustl(spfname(is))
    read(50,*,err=20) natoms(is)
    if (natoms(is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,*)
      stop
    end if
    if (natoms(is).gt.maxatoms) then
      write(*,*)
      write(*,'("Error(readinput): natoms too large : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,'("Adjust maxatoms in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do ia=1,natoms(is)
      read(50,'(A256)',err=20) str
      str=trim(str)//' 0.0 0.0 0.0'
      read(str,*,err=20) atposl(:,ia,is),bfcmt0(:,ia,is)
    end do
  end do
case('plot1d')
  read(50,*,err=20) nvp1d,npp1d
  if (nvp1d.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nvp1d < 1 : ",I8)') nvp1d
    write(*,*)
    stop
  end if
  if (npp1d.lt.nvp1d) then
    write(*,*)
    write(*,'("Error(readinput): npp1d < nvp1d : ",2I8)') npp1d,nvp1d
    write(*,*)
    stop
  end if
  if (allocated(vvlp1d)) deallocate(vvlp1d)
  allocate(vvlp1d(3,nvp1d))
  do iv=1,nvp1d
    read(50,*,err=20) vvlp1d(:,iv)
  end do
case('plot2d')
  read(50,*,err=20) vclp2d(:,1)
  read(50,*,err=20) vclp2d(:,2)
  read(50,*,err=20) vclp2d(:,3)
  read(50,*,err=20) np2d(:)
  if ((np2d(1).lt.1).or.(np2d(2).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np2d < 1 : ",2I8)') np2d
    write(*,*)
    stop
  end if
case('plot3d')
  read(50,*,err=20) vclp3d(:,1)
  read(50,*,err=20) vclp3d(:,2)
  read(50,*,err=20) vclp3d(:,3)
  read(50,*,err=20) vclp3d(:,4)
  read(50,*,err=20) np3d(:)
  if ((np3d(1).lt.1).or.(np3d(2).lt.1).or.(np3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np3d < 1 : ",3I8)') np3d
    write(*,*)
    stop
  end if
case('dos')
  read(50,*,err=20) nwdos,ngrdos,nsmdos
  if (nwdos.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nwdos < 2 : ",I8)') nwdos
    write(*,*)
    stop
  end if
  if (ngrdos.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ngrdos < 1 : ",I8)') ngrdos
    write(*,*)
    stop
  end if
  if (nsmdos.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nsmdos < 0 : ",I8)') nsmdos
    write(*,*)
    stop
  end if
  read(50,*,err=20) wdos(:)
  if (wdos(1).ge.wdos(2)) then
    write(*,*)
    write(*,'("Error(readinput): wdos(1) >= wdos(2) : ",2G18.10)') wdos
    write(*,*)
    stop
  end if
case('dosocc')
  read(50,*,err=20) dosocc
case('dosmsum')
  read(50,*,err=20) dosmsum
case('dosssum')
  read(50,*,err=20) dosmsum
case('lmirep')
  read(50,*,err=20) lmirep
case('tau0atm')
  read(50,*,err=20) tau0atm
case('nstfsp')
  read(50,*,err=20) nstfsp
  if (nstfsp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nstfsp <= 0 : ",I8)') nstfsp
    write(*,*)
    stop
  end if
case('lradstp')
  read(50,*,err=20) lradstp
  if (lradstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
    write(*,*)
    stop
  end if
case('chgexs')
  read(50,*,err=20) chgexs
case('nprad')
  read(50,*,err=20) nprad
  if (nprad.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nprad < 2 : ",I8)') nprad
    write(*,*)
    stop
  end if
case('scissor')
  read(50,*,err=20) scissor
case('optcomp')
  do i=1,27
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty optical component list")')
        write(*,*)
        stop
      end if
      noptcomp=i-1
      goto 10
    end if
    str=trim(str)//' 1'
    read(str,*,iostat=iostat) optcomp(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading optical component list")')
      write(*,'("(blank line required after optcomp block)")')
      write(*,*)
      stop
    end if
    if ((optcomp(1,i).lt.1).or.(optcomp(1,i).gt.3).or. &
        (optcomp(2,i).lt.1).or.(optcomp(2,i).gt.3).or. &
        (optcomp(3,i).lt.1).or.(optcomp(3,i).gt.3)) then
      write(*,*)
      write(*,'("Error(readinput): invalid optcomp : ",3I8)') optcomp
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): optical component list too long")')
  write(*,*)
  stop
case('usegdft')
  read(50,*,err=20) usegdft
case('intraband')
  read(50,*,err=20) intraband
case('evaltol')
  read(50,*,err=20) evaltol
case('deband')
  read(50,*,err=20) deband
  if (deband.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deband < 0 : ",G18.10)') deband
    write(*,*)
    stop
  end if
case('epsband')
  read(50,*,err=20) epsband
  if (epsband.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsband <= 0 : ",G18.10)') epsband
    write(*,*)
    stop
  end if
case('autolinengy')
  read(50,*,err=20) autolinengy
case('dlefe')
  read(50,*,err=20) dlefe
case('bfieldc')
  read(50,*,err=20) bfieldc0
case('efieldc')
  read(50,*,err=20) efieldc
case('afieldc')
  read(50,*,err=20) afieldc
case('fixspin')
  read(50,*,err=20) fixspin
case('momfix')
  read(50,*,err=20) momfix
case('mommtfix')
  do ias=1,maxspecies*maxatoms
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=iostat) is,ia,mommtfix(:,ia,is)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading muffin-tin fixed spin &
       &moments")')
      write(*,'("(blank line required after mommtfix block")')
      write(*,*)
      stop
    end if
  end do
case('taufsm')
  read(50,*,err=20) taufsm
  if (taufsm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
    write(*,*)
    stop
  end if
case('autormt')
  read(50,*,err=20) autormt
case('rmtapm')
  read(50,*,err=20) rmtapm(:)
  if (rmtapm(1).lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rmtapm(1) < 0 : ",G18.10)') rmtapm(1)
    write(*,*)
    stop
  end if
  if ((rmtapm(2).le.0.d0).or.(rmtapm(2).gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rmtapm(2) not in (0,1] : ",G18.10)') rmtapm(2)
    write(*,*)
    stop
  end if
case('isgkmax')
  read(50,*,err=20) isgkmax
case('nosym')
  read(50,*,err=20) nosym
case('deltaph')
  read(50,*,err=20) deltaph
case('phwrite')
  read(50,*,err=20) nphwrt
  if (nphwrt.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nphwrt <= 0 : ",I8)') nphwrt
    write(*,*)
    stop
  end if
  if (allocated(vqlwrt)) deallocate(vqlwrt)
  allocate(vqlwrt(3,nphwrt))
  do i=1,nphwrt
    read(50,*,err=20) vqlwrt(:,i)
  end do
case('notes')
  do i=1,maxnlns
    read(50,'(A80)') notes(i)
    if (trim(notes(i)).eq.'') then
      notelns=i-1
      goto 10
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many note lines")')
  write(*,*)
  stop
case('tforce')
  read(50,*,err=20) tforce
case('tfibs')
  read(50,*,err=20) tfibs
case('maxitoep')
  read(50,*,err=20) maxitoep
  if (maxitoep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
    write(*,*)
    stop
  end if
case('tauoep')
  read(50,*,err=20) tauoep(:)
  if ((tauoep(1).lt.0.d0).or.(tauoep(2).lt.0.d0).or.(tauoep(3).lt.0.d0)) then
    write(*,*)
    write(*,'("Error(readinput): tauoep < 0 : ",3G18.10)') tauoep
    write(*,*)
    stop
  end if
case('kstlist')
  do i=1,maxkst
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty k-point and state list")')
        write(*,*)
        stop
      end if
      nkstlist=i-1
      goto 10
    end if
    str=trim(str)//' 1 1'
    read(str,*,iostat=iostat) kstlist(:,i)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading k-point and state list")')
      write(*,'("(blank line required after kstlist block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): k-point and state list too long")')
  write(*,*)
  stop
case('vklem')
  read(50,*,err=20) vklem
case('deltaem')
  read(50,*,err=20) deltaem
case('ndspem')
  read(50,*,err=20) ndspem
  if ((ndspem.lt.1).or.(ndspem.gt.3)) then
    write(*,*)
    write(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
    write(*,*)
    stop
  end if
case('nosource')
  read(50,*,err=20) nosource
case('spinsprl')
  read(50,*,err=20) spinsprl
case('ssdph')
  read(50,*,err=20) ssdph
case('vqlss')
  read(50,*,err=20) vqlss
case('nwrite')
  read(50,*,err=20) nwrite
case('tevecsv')
  read(50,*,err=20) tevecsv
case('lda+u')
  read(50,*,err=20) ldapu,inptypelu
  do is=1,maxspecies
    read(50,'(A256)',err=20) str
    if (trim(str).eq.'') goto 10
    if (inptypelu.eq.1) then
      read(str,*,iostat=iostat) js,l,ujlu(1:2,js)
    else if (inptypelu.eq.2) then
      read(str,*,iostat=iostat) js,l,(flu(k,js),k=0,2*l,2)
    else if (inptypelu.eq.3) then
      read(str,*,iostat=iostat) js,l,(elu(k,js),k=0,l)
    else if (inptypelu.eq.4) then
      read(str,*,iostat=iostat) js,l,lambdalu(js)
    else if (inptypelu.eq.5) then
      read(str,*,iostat=iostat) js,l,ulufix(js)
    end if
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading LDA+U parameters")')
      write(*,'("(blank line required after lda+u block)")')
      write(*,*)
      stop
    end if
    if ((js.le.0).or.(js.ge.maxspecies)) then
      write(*,*)
      write(*,'("Error(readinput): invalid species number in lda+u block : ", &
       &I8)') js
      write(*,*)
      stop
    end if
    if (l.gt.lmaxlu) then
      write(*,*)
      write(*,'("Error(readinput): l > lmaxlu in lda+u block : ",2I8)') l,lmaxlu
      write(*,*)
      stop
    end if
    llu(js)=l
  end do
case('tmomlu')
  read(50,*,err=20) tmomlu
case('readalu')
  read(50,*,err=20) readalu
case('rdmxctype')
  read(50,*,err=20) rdmxctype
case('rdmmaxscl')
  read(50,*,err=20) rdmmaxscl
  if (rdmmaxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): rdmmaxscl < 0 : ",I8)') rdmmaxscl
    write(*,*)
  end if
case('maxitn')
  read(50,*,err=20) maxitn
case('maxitc')
  read(50,*,err=20) maxitc
case('taurdmn')
  read(50,*,err=20) taurdmn
  if (taurdmn.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmn < 0 : ",G18.10)') taurdmn
    write(*,*)
    stop
  end if
case('taurdmc')
  read(50,*,err=20) taurdmc
  if (taurdmc.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmc < 0 : ",G18.10)') taurdmc
    write(*,*)
    stop
  end if
case('rdmalpha')
  read(50,*,err=20) rdmalpha
  if ((rdmalpha.le.0.d0).or.(rdmalpha.ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rdmalpha not in (0,1) : ",G18.10)') rdmalpha
    write(*,*)
    stop
  end if
case('rdmbeta')
  read(50,*,err=20) rdmbeta
  if ((rdmbeta.le.0.d0).or.(rdmbeta.ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rdmbeta not in (0,1) : ",G18.10)') rdmbeta
    write(*,*)
    stop
  end if
case('rdmtemp')
  read(50,*,err=20) rdmtemp
  if (rdmtemp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rdmtemp < 0 : ",G18.10)') rdmtemp
    write(*,*)
    stop
  end if
case('wrtvnlijji')
  read(50,*,err=20) wrtvnlijji
case('reducebf')
  read(50,*,err=20) reducebf
  if ((reducebf.lt.0.d0).or.(reducebf.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): reducebf not in [0,1] : ",G18.10)') reducebf
    write(*,*)
    stop
  end if
case('ptnucl')
  read(50,*,err=20) ptnucl
case('tseqit')
  read(50,*,err=20) tseqit
case('nseqit')
  read(50,*,err=20) nseqit
  if (nseqit.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nseqit < 1 : ",I8)') nseqit
    write(*,*)
    stop
  end if
case('tauseq')
  read(50,*,err=20) tauseq
  if (tauseq.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tauseq < 0 : ",G18.10)') tauseq
    write(*,*)
    stop
  end if
case('vecql')
  read(50,*,err=20) vecql(:)
case('mustar')
  read(50,*,err=20) mustar
case('sqados')
  read(50,*,err=20) sqados(:)
case('test')
  read(50,*,err=20) test
case('frozencr')
  read(50,*,err=20) frozencr
case('spincore')
  read(50,*,err=20) spincore
case('solscf')
  read(50,*,err=20) solscf
  if (solscf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): solscf < 0 : ",G18.10)') solscf
    write(*,*)
    stop
  end if
case('emaxelnes')
  read(50,*,err=20) emaxelnes
case('lps')
  read(50,*,err=20) natlps
  allocate(lpsrsh(16,16,natlps))
  allocate(iatlps(natlps))
  lpsrsh=0.d0
  do i=1,natlps
    do j=1,16
      lpsrsh(j,j,i)=1.d0
    enddo
    read(50,*,err=20) iatlps(i),l
    do lm1=l**2+1,(l+1)**2
      read(50,*,err=20)(lpsrsh(lm1,lm2,i),lm2=l**2+1,(l+1)**2)
    enddo
  enddo 
case('vqm')
  read(50,*,err=20) nvq
  allocate(vqm(3,nvq))
  do i=1,nvq
    read(50,*,err=20) vqm(:,i)
  end do
case('tq0bz')
  read(50,*,err=20) tq0bz
case('vq0c')
  nvq0=1
  read(50,*,err=20) vq0c(:,1)
  aq0(1)=1.d0/sum(vq0c(:,1)**2)
case('gqmax')
  read(50,*,err=20) gqmax
case('gqsh')
  read(50,*,err=20) gqsh
case('tgqsh')
  read(50,*,err=20) tgqsh 
case('chi0_include_bands')
  read(50,*,err=20) chi0_include_bands(:)
case('chi0_exclude_bands')
  read(50,*,err=20) chi0_exclude_bands(:)  
case('response_w')
  read(50,*,err=20) lr_nw
  read(50,*,err=20) lr_w0,lr_w1,lr_eta  
case('response_iw')
  timgw=.true.
  read(50,*,err=20) lr_nw,lr_iw0,lr_iw1
case('response_fxc')
  read(50,*,err=20) fxctype,nfxca,fxca0,fxca1
case('response_wann')
  read(50,*,err=20) wannier_chi0_chi
  read(50,*,err=20) wannier_chi0_afm
case('wannier')
  read(50,*,err=20) wannier
  read(50,*,err=20) wann_add_poco
  read(50,*,err=20) wann_natom,wann_norbgrp,wann_ntype
  allocate(wann_norb(wann_norbgrp))
  allocate(wann_iorb(3,32,wann_norbgrp))
  allocate(wann_eint(2,wann_ntype))
  allocate(wann_v(wann_ntype))
  do i=1,wann_norbgrp
    read(50,*,err=20) wann_norb(i)
    read(50,*,err=20) (wann_iorb(1,l,i),l=1,wann_norb(i))
    read(50,*,err=20) (wann_iorb(2,l,i),l=1,wann_norb(i))
    read(50,*,err=20) (wann_iorb(3,l,i),l=1,wann_norb(i))
  enddo
  wann_eint=0.d0
  wann_v=0.d0
  do i=1,wann_ntype
      read(50,*,err=20) wann_eint(1,i),wann_eint(2,i),wann_v(i)
  enddo
  allocate(wann_iprj(2,wann_natom))
  do i=1,wann_natom
    read(50,*,err=20) wann_iprj(1,i), wann_iprj(2,i)
  enddo
case ('wannier_min_prjao')
  read(50,*,err=20) wannier_min_prjao
case('wannier_prjao')
  read(50,*,err=20) wannier_prjao
case ('wannier_soft_eint')
  wannier_soft_eint=.true.
  allocate(wannier_soft_eint_e1(wann_ntype))
  allocate(wannier_soft_eint_e2(wann_ntype))
  allocate(wannier_soft_eint_w1(wann_ntype))
  allocate(wannier_soft_eint_w2(wann_ntype))
  do i=1,wann_ntype
    read(50,*,err=20) wannier_soft_eint_e1(i),wannier_soft_eint_e2(i),&
      wannier_soft_eint_w1(i),wannier_soft_eint_w2(i)
  enddo  
case('wannier_lc')
  read(50,*,err=20)wannier_lc
  read(50,*,err=20)nwanlc
  allocate(wanlc_norb(nwanlc))
  allocate(wanlc_iorb(4,wanlcmax,nwanlc))
  allocate(wanlc_iorb_alpha(wanlcmax,nwanlc))
  do n=1,nwanlc
    read(50,*,err=20)wanlc_norb(n)
    if (wanlc_norb(n).gt.wanlcmax) then
      write(*,'("Error(readinput): number of orbitals for linear combination&
      & is too big")')
      stop
    endif
    do i=1,wanlc_norb(n)
      read(50,*,err=20)(wanlc_iorb(l,i,n),l=1,4),wanlc_iorb_alpha(i,n)
    enddo
  enddo
case('wannier_plot')
  read(50,*,err=20)(zero3d(i),i=1,3)
  read(50,*,err=20)(bound3d(i,1),i=1,3)
  read(50,*,err=20)(bound3d(i,2),i=1,3)
  read(50,*,err=20)(bound3d(i,3),i=1,3)
  read(50,*,err=20)(nrxyz(i),i=1,3)
  read(50,*,err=20)nwfplot,firstwf  
case('megqwan_dist')
  read(50,*,err=20) megqwan_mindist,megqwan_maxdist
case('megqwan_cutoff')
  read(50,*,err=20) megqwan_cutoff(:)
case('megqwan_include')
  read(50,*,err=20) nwann_include
  allocate(iwann_include(nwann_include))
  read(50,*,err=20) (iwann_include(i),i=1,nwann_include)
case('megqwan_channels')
  read(50,*,err=20) nwfch
  do i=1,nwfch
    read(50,*,err=20) wfch(1,i), wfch(2,i)
  enddo
  lrespwffilter=.true.
case('mpigrid')
  read(50,'(A256)',err=20) str
  str=adjustl(str)
  j=0
  do
    i=index(trim(str)," ")
    if (i.gt.0) then 
      str1=str(:i)
    else
      str1=str
    endif
    j=j+1
    if (j.gt.mpigrid_maxndim) then
      write(*,'("Error(readinput): number of mpigrid dimensions exceeded")')
      stop
    endif
    read(str1,*,iostat=iostat)mpigrid(j)
    if (iostat.ne.0) then
      write(*,'("Error(readinput): error parsing mpigrid")')
      stop
    endif
    str=adjustl(str(i+1:))
    if (i.le.0) exit
  end do
  mpigrid_ndim=j
  lmpigrid=.true.
case('sic')
  read(50,*,err=20) sic
case('sic_apply')
  allocate(sicw(32,wann_natom))
  do i=1,wann_natom
    read(50,*,err=20) (sicw(j,i),j=1,wann_norb(wann_iprj(2,i)))
  enddo
case('nsclsic')
  read(50,*,err=20) nsclsic
case('sic_sphere')
  read(50,*,err=20) lmaxwan,s_nr
case('sic_smesh')
  read(50,*,err=20) sic_smesh_n,sic_bsht_niter
case('sic_rpole')
  read(50,*,err=20) s_nrpole
  s_nrpole=s_nrpole+1
  allocate(s_rpole(s_nrpole))
  s_rpole(1)=0.d0
  read(50,*,err=20)(s_rpole(i),i=2,s_nrpole)
case('sic_debug_level')
  read(50,*,err=20) sic_debug_level
case('disentangle')
  read(50,*,err=20) ldisentangle  
case('scvl')
  read(50,*,err=20) scvl(:,1)
  read(50,*,err=20) scvl(:,2)
  read(50,*,err=20) scvl(:,3)  
case('papi')
  do i=1,maxpapievents
    read(50,'(A256)') papievent(i)
    if (trim(papievent(i)).eq.'') then
      npapievents=i-1
      goto 10
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many papi events")')
  write(*,*)
  stop
case('sic_wan_cutoff')
  read(50,*,err=20) sic_wan_cutoff
case('sic_me_cutoff')
  read(50,*,err=20) sic_me_cutoff
case('debug_level')
  read(50,*,err=20) debug_level  
case('screenu4')
  read(50,*,err=20) screenu4
case('')
  goto 10
case default
  write(*,*)
  write(*,'("Error(readinput): invalid block name : ",A)') trim(block)
  write(*,*)
  stop
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
! scale the speed of light
solsc=sol*solscf
! scale the lattice vectors (scaling not referenced again in code)
avec(:,1)=sc1*avec(:,1)
avec(:,2)=sc2*avec(:,2)
avec(:,3)=sc3*avec(:,3)
avec(:,:)=sc*avec(:,:)
! case of isolated molecule
if (molecule) then
! convert atomic positions from Cartesian to lattice coordinates
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(ainv,atposl(:,ia,is),v)
      atposl(:,ia,is)=v(:)
    end do
  end do
end if
! read in atomic species data
call readspecies
! delete TEST.OUT
!open(50,file='TEST.OUT')
!close(50,status='DELETE')
return
end subroutine
!EOC
