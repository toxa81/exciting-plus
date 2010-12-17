subroutine sic_main
use modmain
use mod_nrkp
use mod_hdf5
use mod_sic
use mod_wannier
use mod_linresp
implicit none
integer n,sz,i,j,i1,j1,j2,n1,n2,ik,ispn,vtrl(3),ikloc
real(8) t1,t2,t3,vtrc(3)
real(8) etot_,ekin_
integer vl(3)
complex(8) z1
real(8), allocatable :: laplsv(:) 
complex(8), allocatable :: vwank(:,:)
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: ene(:,:)

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
  sz=16.d0*sz*sic_orbitals%ntr*nspinor*(2*sic_wantran%nwan)/1024/1024
  write(151,*)
  write(151,'("Required memory for real-space arrays (MB) : ",I6)')sz
  write(151,'("number of translations for real-space functions : ",I4)')sic_orbitals%ntr
!  do i=1,ntr
!    write(151,'("  i : ",I4,"    vtl(i) : ",3I4)')i,vtl(:,i)
!  enddo
  write(151,*)
  write(151,'("number of included Wannier functions : ",I4)')sic_wantran%nwan
  do j=1,sic_wantran%nwan
    write(151,'("  j : ",I4,"    n : ",I4)')j,sic_wantran%iwan(j)
  enddo
  write(151,*)
  write(151,'("cutoff radius for Wannier functions : ",F12.6)')sic_wan_cutoff
  write(151,'("cutoff radius for SIC matrix elements : ",F12.6)')sic_me_cutoff
  write(151,'("number of Wannier transitions : ",I6)')sic_wantran%nwt
  write(151,*)
  write(151,'("LDA energies of Wannier functions")')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(151,'("  n : ",I4,"    sic_wann_e0 : ",F12.6)')n,sic_wann_e0(n)
  enddo
  call flushifc(151)
endif
if (wproc) then
  write(151,*)
  write(151,'(80("="))')
  write(151,'("generating wave-functions for all k-points")')
  write(151,'(80("="))')
endif
! generate wave-functions for all k-points in BZ
call genwfnr(151,.false.)  
! compute kinetic energy
allocate(laplsv(nstsv))
ekin_=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call genlapl(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
    wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc),laplsv)
  do j=1,nstsv
    ekin_=ekin_-0.5d0*laplsv(j)*occsvnr(j,ik)
  enddo
enddo
call mpi_grid_reduce(ekin_,dims=(/dim_k/))
ekin_=ekin_/nkptnr
if (wproc) then
  write(151,*)
  etot_=ekin_+0.5d0*engyvcl+engyx+engyc+engymad-sic_epot
  write(151,'("kinetic energy : ",G18.10)')ekin_
  write(151,'("potential      : ",G18.10)')0.5d0*engyvcl+engyx+engyc
  write(151,'("madelung       : ",G18.10)')engymad
  write(151,'("sic correction : ",G18.10)')sic_epot
  write(151,'("total energy   : ",G18.10)')etot_
  open(152,file="SIC_ETOT.OUT",form="FORMATTED",status="REPLACE")
  write(152,'(G18.10)')etot_
  close(152)
endif

call sic_wan(151)
allocate(ene(4,sic_wantran%nwan))
call sic_pot(151,ene)
deallocate(ene)
! save old matrix elements
allocate(vwanme_old(sic_wantran%nwt))
vwanme_old=vwanme
! compute matrix elements of SIC potential
!  vwanme = <w_n|v_n|w_{n1,T}>
vwanme=zzero
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  do ispn=1,nspinor    
    vwanme(i)=vwanme(i)+sic_dot_ll(sic_orbitals%wvmt(1,1,1,ispn,j),&
      sic_orbitals%wvir(1,1,ispn,j),sic_orbitals%wanmt(1,1,1,ispn,j1),&
      sic_orbitals%wanir(1,1,ispn,j1),vl,sic_orbitals%twanmt(1,1,n),&
      sic_orbitals%twanmt(1,1,n1))
  enddo
enddo
t1=0.d0
t2=-1.d0
t3=0.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
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
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,'("Matrix elements with maximum difference : ",2I6)')i1,j1
  write(151,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,i1),dreal(vwanme(i1)),dimag(vwanme(i1))
  write(151,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,j1),dreal(vwanme(j1)),dimag(vwanme(j1))
  write(151,*)
  write(151,'("diagonal matrix elements (<W_n|V_n|W_n>) :")')
  !write(151,'(2X,"wann",18X,"V_n")')
  !write(151,'(44("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    i=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'("  n : ",I4,8X,2G18.10)')n,dreal(vwanme(i)),dimag(vwanme(i))
  enddo  
  t3=sqrt(t3/sic_wantran%nwt)
  write(151,*)
  write(151,'("SIC matrix elements RMS difference :",G18.10)')t3  
  write(151,*)
  call flushifc(151)
endif
deallocate(vwanme_old)
! check hermiticity of V_nn'(k)
allocate(vwank(sic_wantran%nwan,sic_wantran%nwan))
do ik=1,nkpt
  vwank=zzero
  do i=1,sic_wantran%nwt
    n1=sic_wantran%iwt(1,i)
    j1=sic_wantran%idxiwan(n1)
    n2=sic_wantran%iwt(2,i)
    j2=sic_wantran%idxiwan(n2)
    vtrl(:)=sic_wantran%iwt(3:5,i)
    vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
    z1=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
    vwank(j1,j2)=vwank(j1,j2)+z1*vwanme(i)
  enddo
  t1=0.d0
  do j1=1,sic_wantran%nwan
    do j2=1,sic_wantran%nwan
      t1=max(t1,abs(vwank(j1,j2)-dconjg(vwank(j2,j1))))
    enddo
  enddo
  if (wproc) then
    write(151,'("ik : ",I4,"   max.herm.err : ",G18.10 )')ik,t1
  endif
enddo
deallocate(vwank)
! signal that now we have computed sic potential and wannier functions
tsic_wv=.true.
! write to HDF5 file after last iteration
if (isclsic.eq.nsclsic) call sic_writevwan
if (wproc) then
  call timestamp(151,"Done.")
  close(151)
endif
return
end
