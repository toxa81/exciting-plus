subroutine wann_plot
use modmain
use mod_nrkp
implicit none

real(8) orig(3)
complex(4), allocatable :: wf(:,:,:)
complex(4), allocatable :: wfval(:,:)
complex(4), allocatable :: wfp(:)
integer i,nrtot
integer i1,i2,i3,ir,n,m
real(8), allocatable :: vr(:,:)
real(8), allocatable :: veff(:)
complex(8), allocatable :: zfft_vir(:)
character*40 fname
real(8) x(2),alph,t1,x0,dx
logical, parameter :: wfprod=.false.
integer ikloc

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

call genwfnr(-1,.false.)

if (task.eq.361) then
  nrtot=nrxyz(1)
  nrxyz(2)=1
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1))/2.d0
endif
if (task.eq.362) then
  nrtot=nrxyz(1)*nrxyz(2)
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2))/2.d0
endif
if (task.eq.363) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
endif

if (mpi_grid_root()) then
  allocate(wf(nspinor,nwfplot,nrtot))
  allocate(wfp(nrtot))
  allocate(veff(nrtot))
endif
allocate(vr(3,nrtot))
allocate(wfval(nspinor,nwfplot))

! make (1,2,3)D-grid of r-points
ir=0
do i1=0,nrxyz(1)-1
  do i2=0,nrxyz(2)-1
    do i3=0,nrxyz(3)-1
      ir=ir+1
      vr(:,ir)=orig(:)+i1*bound3d(:,1)/nrxyz(1)+&
                       i2*bound3d(:,2)/nrxyz(2)+&
                       i3*bound3d(:,3)/nrxyz(3)
    enddo
  enddo
enddo

! Fourier transform potential to G-space
allocate(zfft_vir(ngrtot))
zfft_vir(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft_vir)

do ir=1,nrtot
  if (mod(ir,nrxyz(2)*nrxyz(3)).eq.0.and.mpi_grid_root()) then
    write(*,*)'r-point : ',ir,' out of ',nrtot
  endif
  call wann_val(vr(1,ir),wfval)
  if (mpi_grid_root()) wf(:,:,ir)=wfval(:,:)
  call f_veff_p(vr(1,ir),veffmt,zfft_vir,t1)
  if (mpi_grid_root()) veff(ir)=t1
  call mpi_grid_barrier()
enddo

if (mpi_grid_root()) then
  if (task.eq.362) then
    x(1)=sqrt(bound3d(1,1)**2+bound3d(2,1)**2+bound3d(3,1)**2)
    x(2)=sqrt(bound3d(1,2)**2+bound3d(2,2)**2+bound3d(3,2)**2)
    alph=dot_product(bound3d(:,1),bound3d(:,2))/x(1)/x(2)
    alph=acos(alph)
  endif

  if (task.eq.361) then
    x0=-sqrt(sum(bound3d(:,1)**2))/2.d0
    dx=sqrt(sum(bound3d(:,1)**2))/nrxyz(1)
    do n=1,nwfplot
      write(fname,'("wf_",I3.3,".dat")')n+firstwf-1
      open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
      do ir=1,nrtot
        write(70,'(2G18.10)')x0+(ir-1)*dx,sum(abs(wf(:,n,ir))**2)
      enddo
      close(70)
    enddo
    open(70,file="veff.dat",status="REPLACE",form="FORMATTED")
    do ir=1,nrtot
      write(70,'(2G18.10)')x0+(ir-1)*dx,veff(ir)
    enddo
    close(70)
  endif  
  if (task.eq.362.or.task.eq.363) then
    do n=1,nwfplot
      write(fname,'("wf_",I3.3,".dx")')n+firstwf-1
      open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
      if (task.eq.362) then
        write(70,'("object 1 class gridpositions counts",2I4)')nrxyz(1),nrxyz(2)
        write(70,'("origin ",2G18.10)')(/0.d0, 0.d0/)
        write(70,'("delta ",2G18.10)')(/x(1),0.d0/)/nrxyz(1)
        write(70,'("delta ",2G18.10)')(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
        write(70,'("object 2 class gridconnections counts",2I4)')nrxyz(1),nrxyz(2)
      endif
      if (task.eq.363) then
        write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
        write(70,'("origin ",3G18.10)')orig(:)
        do i=1,3
          write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
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
      write(70)(sum(abs(wf(:,n,ir))**2),ir=1,nrtot)
      close(70)
    enddo !n
    fname="veff.dx"
    open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
    if (task.eq.362) then
      write(70,'("object 1 class gridpositions counts",2I4)')nrxyz(1),nrxyz(2)
      write(70,'("origin ",2G18.10)')(/0.d0, 0.d0/)
      write(70,'("delta ",2G18.10)')(/x(1),0.d0/)/nrxyz(1)
      write(70,'("delta ",2G18.10)')(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
      write(70,'("object 2 class gridconnections counts",2I4)')nrxyz(1),nrxyz(2)
    endif
    if (task.eq.363) then
      write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
      write(70,'("origin ",3G18.10)')orig(:)
      do i=1,3
        write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
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
    write(70)(sngl(veff(ir)),ir=1,nrtot)
    close(70)
  endif
endif
return
end
