subroutine wann_plot_3d
use modmain
use mod_nrkp
implicit none

real(8) orig(3)
complex(4), allocatable :: wf(:,:,:)
real(4), allocatable :: wf2(:)
integer i,nrtot
integer i1,i2,i3,ir,n
real(8) vrc(3)
real(4), allocatable :: veff(:)
real(8) t1
complex(8), allocatable :: zfft(:)
character*40 fname
logical, parameter :: wfprod=.false.
integer recl

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

call genwfnr(-1,.false.,lmaxvr)

if (task.eq.863) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
endif

allocate(wf(nspinor,nwfplot,nrxyz(2)*nrxyz(3)))
allocate(veff(nrxyz(2)*nrxyz(3)))
allocate(wf2(nrxyz(2)*nrxyz(3)))


! Fourier transform potential to G-space
allocate(zfft(ngrtot))
zfft(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft)

!inquire(iolength=recl)veff
recl=4*nrxyz(2)*nrxyz(3)

if (mpi_grid_root()) then
  do n=1,nwfplot
    write(fname,'("wf_",I3.3,".dx")')n+firstwf-1
    open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
    write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,'("origin ",3G18.10)')orig(:)
    do i=1,3
      write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
    enddo
    write(70,'("object 2 class gridconnections counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,'("object 3 class array type float rank 1 shape 1 items ",&
      &I8," lsb ieee data file wf",I3.3,".bin,0")')nrxyz(1)*nrxyz(2)*nrxyz(3),n
    write(70,'("object ""regular positions regular connections"" class field")')
    write(70,'("component ""positions"" value 1")')
    write(70,'("component ""connections"" value 2")')
    write(70,'("component ""data"" value 3")')
    write(70,'("end")')
    close(70)
    write(fname,'("wf",I3.3,".bin")')n
    open(70+n,file=trim(fname),form="unformatted",status="replace",access="direct",recl=recl)
  enddo !n
  fname="veff.dx"
  open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
  write(70,'("object 1 class gridpositions counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
  write(70,'("origin ",3G18.10)')orig(:)
  do i=1,3
    write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
  enddo
  write(70,'("object 2 class gridconnections counts",3I4)')nrxyz(1),nrxyz(2),nrxyz(3)
  write(70,'("object 3 class array type float rank 1 shape 1 items ",&
    &I8," lsb ieee data file veff.bin,0")')nrxyz(1)*nrxyz(2)*nrxyz(3)
  write(70,'("object ""regular positions regular connections"" class field")')
  write(70,'("component ""positions"" value 1")')
  write(70,'("component ""connections"" value 2")')
  write(70,'("component ""data"" value 3")')
  write(70,'("end")')
  close(70)
  open(70,file="veff.bin",form="unformatted",status="replace",access="direct",recl=recl)
endif
do i1=0,nrxyz(1)-1
  if (mpi_grid_root()) write(*,'("slab ",I4," out of ",I4)')i1+1,nrxyz(1)
  ir=1
  do i2=0,nrxyz(2)-1
    do i3=0,nrxyz(3)-1
      vrc(:)=orig(:)+i1*bound3d(:,1)/nrxyz(1)+&
                     i2*bound3d(:,2)/nrxyz(2)+&
                     i3*bound3d(:,3)/nrxyz(3)
      call wann_val(vrc,wf(1,1,ir))
      call rfval(vrc,lmaxvr,lmmaxvr,veffmt,zfft,t1)
      veff(ir)=sngl(t1)
      ir=ir+1
    enddo
  enddo
  call mpi_grid_reduce(wf(1,1,1),nspinor*nwfplot*nrxyz(2)*nrxyz(3),dims=(/dim_k/))
  call mpi_grid_reduce(veff(1),nrxyz(2)*nrxyz(3),dims=(/dim_k/))
  if (mpi_grid_root()) then
    do n=1,nwfplot
      do ir=1,nrxyz(2)*nrxyz(3)
        wf2(ir)=sum(abs(wf(:,n,ir))**2)
      enddo 
      write(70+n,rec=i1+1)wf2
    enddo
    write(70,rec=i1+1)veff  
  endif
  call mpi_grid_barrier()
enddo
if (mpi_grid_root()) then
  do n=0,nwfplot
    close(70+n)
  enddo
endif  
return
end
