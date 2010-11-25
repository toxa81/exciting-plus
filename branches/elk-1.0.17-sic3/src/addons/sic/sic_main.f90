subroutine sic_main
use modmain
use mod_nrkp
use mod_hdf5
use mod_sic
implicit none
integer n,sz,i,j,i1,j1,n1,ispn,vtrl(3)
real(8) t1,t2,t3,vtrc(3)
integer vl(3)
! Wannier functions
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: ene(:,:)
complex(8), allocatable :: vwank(:,:)
complex(8) z1
logical exist
integer n2,ik

sic=.true.

! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) return
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp

wproc=mpi_grid_root()
if (wproc) then
  open(151,file="SIC.OUT",form="FORMATTED",status="REPLACE")
endif
if (wproc) then
  sz=lmmaxvr*nmtloc+ngrloc
  sz=16*sz*ntr*nspinor*(2*nwann+2)/1024/1024
  write(151,*)
  write(151,'("Required memory for real-space arrays (MB) : ",I6)')sz
  write(151,*)
  write(151,'("cutoff radius for WF : ",F12.6)')wann_r_cutoff
  write(151,'("number of translations : ",I4)')ntr
  do i=1,ntr
    write(151,'("  i : ",I4,"    vtl(i) : ",3I4)')i,vtl(:,i)
  enddo
  call flushifc(151)
endif
! generate wave-functions for all k-points in BZ
call genwfnr(151,.false.)  
call sic_wan(151)
allocate(ene(4,nwann))
call sic_pot(151,ene)
!----------------------------------!
! matrix elements of SIC potential !
!----------------------------------!
allocate(vwanme_old(nmegqwan))
vwanme_old=vwanme
! compute matrix elements of SIC potential
!  vwanme = <w_n|v_n|w_{n1,T}>
vwanme=zzero
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  do ispn=1,nspinor    
    vwanme(i)=vwanme(i)+sic_dot_ll(wvmt(1,1,1,ispn,n),wvir(1,1,ispn,n),&
      wanmt(1,1,1,ispn,n1),wanir(1,1,ispn,n1),vl,twanmt(1,1,n),twanmt(1,1,n1))
  enddo
enddo
t1=0.d0
t2=-1.d0
t3=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-vl(1),-vl(2),-vl(3))
  t1=t1+abs(vwanme(i)-dconjg(vwanme(j)))
  if (abs(vwanme(i)-dconjg(vwanme(j))).ge.t2) then
    t2=abs(vwanme(i)-dconjg(vwanme(j)))
    i1=i
    j1=j
  endif
  t3=t3+abs(vwanme(i)-vwanme_old(i))**2
enddo
if (wproc) then
  call timestamp(151,"done with matrix elements")
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("Matrix elements of SIC potential &
    &(n n1  <w_n|v_n|w_n1}>)")')
  do i=1,nmegqwan
    vl(:)=imegqwan(3:5,i)
    if (all(vl.eq.0)) then
      write(151,'(I4,4X,I4,4X,2G18.10)')imegqwan(1:2,i),&
        dreal(vwanme(i)),dimag(vwanme(i))
    endif
  enddo
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,'("Average deviation from ""localization criterion"" : ",F12.6)')t1/nmegqwan
  write(151,*)
  write(151,'("Matrix elements with maximum difference : ",2I6)')i1,j1
  write(151,'(I4,4X,I4,4X,3I4,4X,2G18.10)')imegqwan(:,i1),&
        dreal(vwanme(i1)),dimag(vwanme(i1))
  write(151,'(I4,4X,I4,4X,3I4,4X,2G18.10)')imegqwan(:,j1),&
        dreal(vwanme(j1)),dimag(vwanme(j1))
  write(151,*)
  write(151,'("Diagonal matrix elements")')
  write(151,'(2X,"wann",18X,"V_n")')
  write(151,'(44("-"))')
  do n=1,nwann
    j=idxmegqwan(n,n,0,0,0)
    write(151,'(I4,4X,2G18.10)')n,dreal(vwanme(j)),dimag(vwanme(j))
  enddo  
  t3=sqrt(t3/nmegqwan)
  write(151,*)
  write(151,'("SIC matrix elements RMS difference :",G18.10)')t3  
  call flushifc(151)
endif
deallocate(vwanme_old)
! check hermiticity of V_nn'(k)
allocate(vwank(nwann,nwann))
do ik=1,nkpt
  vwank=zzero
  do i=1,nmegqwan
    n1=imegqwan(1,i)
    n2=imegqwan(2,i)
    vtrl(:)=imegqwan(3:5,i)
    vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
    z1=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
    vwank(n1,n2)=vwank(n1,n2)+z1*vwanme(i)
  enddo
  t1=0.d0
  do n1=1,nwann
    do n2=1,nwann
      t1=max(t1,abs(vwank(n1,n2)-dconjg(vwank(n2,n1))))
    enddo
  enddo
  if (wproc) then
    write(151,*)
    write(151,'("ik : ",I4,"   max.herm.err : ",G18.10 )')ik,t1
!    write(151,*)
!    do n1=1,nwann
!      write(151,'(5X,255F12.7)')(dreal(vwank(n1,n2)),n2=1,nwann)
!    enddo
!    write(151,*)
!    do n1=1,nwann
!      write(151,'(5X,255F12.7)')(dimag(vwank(n1,n2)),n2=1,nwann)
!    enddo
  endif
enddo
deallocate(vwank)
if (wproc) close(151)
! flag that now we have computed sic potential and wannier functions
tsic_wv=.true.
! write to HDF5 file after last iteration
if (isclsic.eq.nsclsic) call sic_writevwan
deallocate(ene)
return
end
