subroutine writewann
use modmain
implicit none
integer i,j,ik,ikloc,idm
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
  allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc(iproc)))
  allocate(evecsvloc(nstsv,nstsv,nkptloc(iproc)))
endif
evalsv=0.d0
do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptloc(iproc)
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
do ikloc=1,nkptloc(iproc)
  if (task.eq.600) call genwann_h(ikloc)
  if (task.eq.601) call genwann_p(ikloc,evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc))
enddo
call zsync(wann_h,nwann*nwann*nkpt,.true.,.false.)
call zsync(wann_p,3*nwann*nwann*nkpt,.true.,.false.)
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
if (task.eq.601) then
  deallocate(evecfvloc)
  deallocate(evecsvloc)
endif
return
end
