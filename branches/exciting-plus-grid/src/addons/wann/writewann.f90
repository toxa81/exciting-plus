subroutine writewann
use modmain
implicit none
integer i,j,ik,ikloc,idm,lm,ir,ispn,n,ias,io
complex(8), allocatable :: wann_rf(:,:,:,:)
character*20 fname
integer, external :: ikglob

if (.not.wannier) then
  if (iproc.eq.0) then
    write(*,*)
    write(*,'("Error(writewann_h) : WF generation is switched off")')
    write(*,*)
    call pstop
  endif
endif

call init0
call init1
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call geturf
call genurfprod


if (task.eq.601) then
  allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))
  allocate(evecsvloc(nstsv,nstsv,nkptloc))
endif
if (task.eq.602) then
  allocate(wann_rf(lmmaxvr,nrmtmax,nspinor,nwann))
  wann_rf=zzero
endif
evalsv=0.d0
do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptloc
      call getwann(ikloc)
      if (task.eq.600) then
        call getevalsv(vkl(1,ikglob(ikloc)),evalsv(1,ikglob(ikloc)))
      endif
      if (task.eq.601) then
        call getevecfv(vkl(1,ikglob(ikloc)),vgkl(1,1,1,ikloc),evecfvloc(1,1,1,ikloc))
        call getevecsv(vkl(1,ikglob(ikloc)),evecsvloc(1,1,ikloc))
      endif
   end do
  end if
  call barrier(comm_world)
end do
do ikloc=1,nkptloc
  if (task.eq.600) call genwann_h(ikloc)
  if (task.eq.601) call genwann_p(ikloc,evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc))
  if (task.eq.602) then
    do n=1,nwann
      ias=iwann(1,n)
      do ispn=1,nspinor
        do ir=1,nrmt(ias2is(ias))
          do lm=1,lmmaxvr
            do io=1,nrfmax
              wann_rf(lm,ir,ispn,n)=wann_rf(lm,ir,ispn,n)+&
                urf(ir,lm2l(lm),io,ias)*wann_unkmt(lm,io,ias,ispn,n,ikloc)
            enddo
          enddo !lm
        enddo !ir
      enddo !ispn
    enddo !n
  endif
enddo
call zsync(wann_h,nwann*nwann*nkpt,.true.,.false.)
call zsync(wann_p,3*nwann*nwann*nkpt,.true.,.false.)
if (task.eq.602) then
  do n=1,nwann
    call zsync(wann_rf(1,1,1,n),lmmaxvr*nrmtmax*nspinor,.true.,.false.)
  enddo
endif
if (iproc.eq.0.and.task.eq.600) then
  open(200,file='WANN_H0.OUT',form='formatted',status='replace')
  do i=1,nwann
    write(200,'(6X,255G18.10)')(dreal(sum(wann_h(i,j,:))/nkpt),j=1,nwann)
  enddo
  write(200,*)
  do i=1,nwann
    write(200,'(6X,255G18.10)')(dimag(sum(wann_h(i,j,:))/nkpt),j=1,nwann)
  enddo
  close(200)
  do i=1,nwann
    wann_h(i,i,:)=wann_h(i,i,:)-efermi
  enddo
  wann_h=wann_h*ha2ev
  open(200,file='WANN_Hk.OUT',form='formatted',status='replace')
  write(200,*)nkpt,nwann
  do ik=1,nkpt
    write(200,*)1.d0 !wtkp(ikp)
    do i=1,nwann
      do j=1,nwann
        write(200,*)dreal(wann_h(i,j,ik)),dimag(wann_h(i,j,ik))
      enddo
    enddo
  enddo	
  close(200)
endif
if (iproc.eq.0.and.task.eq.601) then
  open(200,file='WANN_P.OUT',form='formatted',status='replace')
  do ik=1,nkpt
    write(200,'("ik : ",I4)')ik
    do idm=1,3
      write(200,'("  x : ",I4)')idm
      do i=1,nwann
        write(200,'(6X,255G18.10)')(dreal(wann_p(idm,i,j,ik)),j=1,nwann)
      enddo
      write(200,*)
      do i=1,nwann
        write(200,'(6X,255G18.10)')(dimag(wann_p(idm,i,j,ik)),j=1,nwann)
      enddo
    enddo
  enddo
  write(200,*)
  do idm=1,3
    write(200,'("  x : ",I4)')idm
    do i=1,nwann
      write(200,'(6X,255G18.10)')(dreal(sum(wann_p(idm,i,j,:))/nkpt),j=1,nwann)
    enddo
    write(200,*)
    do i=1,nwann
      write(200,'(6X,255G18.10)')(dimag(sum(wann_p(idm,i,j,:))/nkpt),j=1,nwann)
    enddo
  enddo
  close(200)
endif
if (iproc.eq.0.and.task.eq.602) then
  do n=1,nwann
    write(fname,'("WANN_",I3.3,"_rfmt.OUT")')n
    open(200,file=trim(fname),form='formatted',status='replace')
    do lm=1,16
      do ir=1,nrmt(ias2is(iwann(1,n)))
        write(200,'(2G18.10)')spr(ir,ias2is(iwann(1,n))),abs(wann_rf(lm,ir,1,n))
      enddo
      write(200,*)
    enddo
    close(200)
  enddo
endif
if (task.eq.601) then
  deallocate(evecfvloc)
  deallocate(evecsvloc)
endif
if (task.eq.602) then
  deallocate(wann_rf)
endif
return
end
