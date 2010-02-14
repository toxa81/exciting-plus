subroutine init_band_trans
use modmain
implicit none
integer ikloc,ik,i,j,n1,n2,i1,i2,i3
real(8) t1
logical l1
logical, external :: wann_diel

! setup n,n' stuff
!call timer_start(1,reset=.true.)
!if (spinpol) then
!  if (allocated(spinor_ud)) deallocate(spinor_ud)
!  allocate(spinor_ud(2,nstsv,nkptnr))
!  spinor_ud=0
!  do ikloc=1,nkptnrloc
!    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!    do j=1,nstsv
!      t1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
!      if (t1.gt.1d-10) spinor_ud(1,j,ik)=1
!      t1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
!      if (t1.gt.1d-10) spinor_ud(2,j,ik)=1
!    enddo
!  enddo
!  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
!endif
! get limits
call getmeidx(.true.)
if (allocated(nmegqblh)) deallocate(nmegqblh)
allocate(nmegqblh(nkptnrloc))
nmegqblh=0
if (allocated(bmegqblh)) deallocate(bmegqblh)
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
bmegqblh=0
if (wannier_megq) then
  if (allocated(nmegqblhwan)) deallocate(nmegqblhwan)
  allocate(nmegqblhwan(nkptnrloc))
  nmegqblhwan=0
  if (allocated(imegqblhwan)) deallocate(imegqblhwan)
  allocate(imegqblhwan(nmegqblhmax,nkptnrloc))
  imegqblhwan=0
endif
call getmeidx(.false.)
! split interband transitions between second dimension
if (allocated(nmegqblhloc)) deallocate(nmegqblhloc)
allocate(nmegqblhloc(2,nkptnrloc))
nmegqblhloc=0
do ikloc=1,nkptnrloc
  nmegqblhloc(1,ikloc)=mpi_grid_map(nmegqblh(ikloc),dim_b,offs=i)
  nmegqblhloc(2,ikloc)=i
enddo
nmegqblhlocmax=maxval(nmegqblhloc)
!call timer_stop(1)
if (wannier_megq) then
  if (allocated(bmegqwan)) deallocate(bmegqwan)  
  allocate(bmegqwan(2,nwann*nwann))
  bmegqwan=0
  nmegqwan=0
  do n1=1,nwann
    do n2=1,nwann
      l1=.false.
! for integer occupancy numbers take only transitions between occupied and empty bands
      if (wann_diel().and.(abs(wann_occ(n1)-wann_occ(n2)).gt.1d-8)) l1=.true.
! for fractional occupancies or cRPA calculation take all transitions
      if (.not.wann_diel().or.crpa) l1=.true.
      if (l1) then
        nmegqwan=nmegqwan+1
        bmegqwan(1,nmegqwan)=n1
        bmegqwan(2,nmegqwan)=n2
      endif
    enddo
  enddo
! list of translations
  ntrmegqwan=(2*megqwan_maxtr+1)**3
  if (allocated(itrmegqwan)) deallocate(itrmegqwan)
  allocate(itrmegqwan(3,ntrmegqwan))
  i=0
  do i1=-megqwan_maxtr,megqwan_maxtr
    do i2=-megqwan_maxtr,megqwan_maxtr
      do i3=-megqwan_maxtr,megqwan_maxtr
        i=i+1
        itrmegqwan(:,i)=(/i1,i2,i3/)
      enddo
    enddo
  enddo
endif

return
end