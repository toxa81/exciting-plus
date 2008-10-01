subroutine wann_plot
use modmain
use modwann
implicit none

real(8) r(3),r1(3),t(3),d
real(8) bound(3,3),orig(3)
complex(8), allocatable :: wf(:)

integer ntr(3),i,ivec,nrxyz(3),nrtot
integer i1,i2,i3,ir
integer mtord
real(8), allocatable :: ufr(:,:,:,:)

call init0
call init1

call wann_init

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

call get_a_ort

! Cartesian coordinates of boundary and origin
bound(:,1)=(/1.d0,0.d0,0.d0/)
bound(:,2)=(/0.d0,1.d0,0.d0/)
bound(:,3)=(/0.d0,0.d0,1.d0/)
orig(:)=(/0.d0,0.d0,0.d0/)
nrxyz(:)=(/20,20,20/)
nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
allocate(wf(nrtot))


call getmtord(lmaxapw,mtord)
allocate(ufr(nrmtmax,0:lmaxapw,mtord,natmtot))
call getufr(lmaxapw,mtord,ufr)


ir=0
do i1=0,nrxyz(1)-1
  do i2=0,nrxyz(2)-1
    do i3=0,nrxyz(3)-1
! arbitrary r-point
      r(:)=orig(:)+i1*bound(:,1)/nrxyz(1)+ &
                   i2*bound(:,2)/nrxyz(2)+ &
                   i3*bound(:,3)/nrxyz(3)
      ir=ir+1
      call wann_val(r,wf(ir))
    enddo
  enddo
enddo
      








return
end


subroutine wann_val(r,val)
use modmain
use modwann
implicit none
! arguments
real(8), intent(in) :: r(3)
complex(8), intent(out) :: val

integer ivec,ntr(3)
real t(3),d,r1(3)

! reduce r-point to primitive cell and find corresponding translation vector
do ivec=1,3
  d=dot_product(r,bvec(:,ivec))/twopi
  ntr(ivec)=floor(d)
enddo

t(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)

r1(:)=r(:)-t(:)

return
end

subroutine put_a_ort
use modmain
use modwann
implicit none
integer ik,n,j,ispn

open(70,file='A_ORT.OUT',form='unformatted',status='replace')
do ik=1,nkpt
  write(70)(((a_ort(n,j,ispn,ik),n=1,wf_dim),j=1,nstfv),ispn=1,wann_nspins)
enddo
close(70)

return
end

subroutine get_a_ort
use modmain
use modwann
implicit none
integer ik,n,j,ispn

open(70,file='A_ORT.OUT',form='unformatted',status='old')
do ik=1,nkpt
  read(70)(((a_ort(n,j,ispn,ik),n=1,wf_dim),j=1,nstfv),ispn=1,wann_nspins)
enddo
close(70)

return
end



