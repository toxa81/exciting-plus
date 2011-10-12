subroutine wann_plot
use modmain
use mod_nrkp
use mod_sic
use mod_madness
use modxcifc
implicit none
real(8) orig(3)
real(4), allocatable :: func(:)
integer i,nrtot,nrloc,irloc
integer i1,i2,i3,ir,n,j,ispn
real(8), allocatable :: vr(:,:)
real(8), allocatable :: veff(:)
complex(8), allocatable :: zfft(:)
real(8) vr1(3),vr2(3),vr3(3)
character*40 fname
real(8) x(2),alph,x0,dx
logical, parameter :: wfprod=.false.
complex(8) zval(2)
real(8) dens(1),e1(1),e2(1),v1(1),v2(1) 
!
call init0
call init1
if (.not.mpi_grid_in()) return
wproc=mpi_grid_root()
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp

call genwfnr(-1,.false.,lmaxapw)
call elk_m_init

if (task.eq.861) then
  nrtot=nrxyz(1)
  nrxyz(2)=1
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1))/2.d0
endif
if (task.eq.862) then
  nrtot=nrxyz(1)*nrxyz(2)
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2))/2.d0
endif
if (task.eq.863) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
endif

allocate(vr(3,nrtot))
allocate(func(nrtot))
allocate(veff(nrtot))

! make (1,2,3)D-grid of r-points
ir=0
do i1=0,nrxyz(1)-1
  vr1=0.d0
  if (nrxyz(1).ne.1) vr1(:)=dble(i1)*bound3d(:,1)/(nrxyz(1)-1)
  do i2=0,nrxyz(2)-1
    vr2=0.d0
    if (nrxyz(2).ne.1) vr2(:)=dble(i2)*bound3d(:,2)/(nrxyz(2)-1)
    do i3=0,nrxyz(3)-1
      vr3=0.d0
      if (nrxyz(3).ne.1) vr3(:)=dble(i3)*bound3d(:,3)/(nrxyz(3)-1)
      ir=ir+1
      vr(:,ir)=orig(:)+vr1(:)+vr2(:)+vr3(:)
    enddo
  enddo
enddo

! Fourier transform potential to G-space
allocate(zfft(ngrtot))
zfft(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft)

nrloc=mpi_grid_map(nrtot,dim1)

do j=1,nwfplot
  n=j+firstwf-1
  call elk_load_wann_unk(n)
  func=0.0 
  do irloc=1,nrloc
    ir=mpi_grid_map(nrtot,dim1,loc=irloc)
    !if (mod(ir,nrxyz(2)*nrxyz(3)).eq.0.and.mpi_grid_root()) then
    !  write(*,*)'r-point : ',ir,' out of ',nrtot
    !endif
    call s_get_wanval(vr(:,ir),zval)
    dens=0.d0
    do ispn=1,nspinor
      dens(1)=dens(1)+abs(zval(ispn))**2
    enddo
    call xcifc(xctype,n=1,rho=dens,ex=e1,ec=e2,vx=v1,vc=v2)
    func(ir)=sngl(dens(1)) !sngl(abs(zval(1)*(v1(1)+v2(1)))) !sngl(dens(1))
    !call rfval(vr(1,ir),lmaxvr,lmmaxvr,veffmt,zfft,t1)
  enddo !ir
  call mpi_grid_reduce(func(1),nrtot,dims=(/dim1/))
! write to file
  if (mpi_grid_root()) then
    if (task.eq.862) then
      x(1)=sqrt(bound3d(1,1)**2+bound3d(2,1)**2+bound3d(3,1)**2)
      x(2)=sqrt(bound3d(1,2)**2+bound3d(2,2)**2+bound3d(3,2)**2)
      alph=dot_product(bound3d(:,1),bound3d(:,2))/x(1)/x(2)
      alph=acos(alph)
    endif
    if (task.eq.861) then
      x0=-sqrt(sum(bound3d(:,1)**2))/2.d0
      dx=sqrt(sum(bound3d(:,1)**2))/nrxyz(1)
      write(fname,'("wf_",I3.3,".dat")')n
      open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
      do ir=1,nrtot
        write(70,'(2G18.10)')x0+(ir-1)*dx,func(ir)
      enddo
      close(70)
      !open(70,file="veff.dat",status="REPLACE",form="FORMATTED")
      !do ir=1,nrtot
      !  write(70,'(2G18.10)')x0+(ir-1)*dx,veff(ir)
      !enddo
      !close(70)
    endif  
    if (task.eq.862.or.task.eq.863) then
      write(fname,'("wf_",I3.3,".dx")')n
      open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
      if (task.eq.862) then
        write(70,'("object 1 class gridpositions counts",2I4)')nrxyz(1),nrxyz(2)
        write(70,'("origin ",2G18.10)')(/0.d0, 0.d0/)
        write(70,'("delta ",2G18.10)')(/x(1),0.d0/)/(nrxyz(1)-1)
        write(70,'("delta ",2G18.10)')(/x(2)*cos(alph),x(2)*sin(alph)/)/(nrxyz(2)-1)
        write(70,'("object 2 class gridconnections counts",2I4)')nrxyz(1),nrxyz(2)
      endif
      if (task.eq.863) then
        write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
        write(70,'("origin ",3G18.10)')orig(:)
        do i=1,3
          write(70,'("delta ",3G18.10)')bound3d(:,i)/(nrxyz(i)-1)
        enddo
        write(70,'("object 2 class gridconnections counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
      endif
      write(70,'("object 3 class array type float rank 1 shape 1 items ",&
        &I8," lsb ieee data 4")')nrtot
      write(70,'("object ""regular positions regular connections"" class field")')
      write(70,'("component ""positions"" value 1")')
      write(70,'("component ""connections"" value 2")')
      write(70,'("component ""data"" value 3")')
      write(70,'("end")')
      close(70)
      open(70,file=trim(fname),status="OLD",form="UNFORMATTED",position="APPEND")
      write(70)(func(ir),ir=1,nrtot)
      close(70)
    endif
  endif
enddo !j

!    fname="veff.dx"
!    open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
!    if (task.eq.862) then
!      write(70,'("object 1 class gridpositions counts",2I4)')nrxyz(1),nrxyz(2)
!      write(70,'("origin ",2G18.10)')(/0.d0, 0.d0/)
!      write(70,'("delta ",2G18.10)')(/x(1),0.d0/)/nrxyz(1)
!      write(70,'("delta ",2G18.10)')(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
!      write(70,'("object 2 class gridconnections counts",2I4)')nrxyz(1),nrxyz(2)
!    endif
!    if (task.eq.863) then
!      write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
!      write(70,'("origin ",3G18.10)')orig(:)
!      do i=1,3
!        write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
!      enddo
!      write(70,'("object 2 class gridconnections counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
!    endif
!    write(70,'("object 3 class array type float rank 1 shape 1 items ",&
!      &I8," lsb ieee data 4")')nrtot
!    write(70,'("object ""regular positions regular connections"" class field")')
!    write(70,'("component ""positions"" value 1")')
!    write(70,'("component ""connections"" value 2")')
!    write(70,'("component ""data"" value 3")')
!    write(70,'("end")')
!    close(70)
!    open(70,file=trim(fname),status="OLD",form="UNFORMATTED",position="APPEND")
!    write(70)(sngl(veff(ir)),ir=1,nrtot)
!    close(70)
!  endif
!endif


return
end
