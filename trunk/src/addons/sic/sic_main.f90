subroutine sic_main
use modmain
use mod_nrkp
use mod_hdf5
use mod_sic
use mod_wannier
use mod_linresp
use mod_madness
implicit none
integer n,i,j,i1,j1,j2,n1,n2,ik,ispn,vtrl(3),ikloc,ig,nwtloc,iloc
integer ias,lm
real(8) t1,t2,t3,vtrc(3),pos1(3),pos2(3)
integer vl(3)
complex(8) z1
real(8), allocatable :: laplsv(:) 
complex(8), allocatable :: vwank(:,:)
complex(8), allocatable :: vwanme_old(:)
character*20 c1,c2,c3
character, parameter :: orbc(4)=(/'s','p','d','f'/)
character*2, parameter :: spinc(2)=(/'up','dn'/)
!
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
  write(151,*)
  write(151,'("number of included Wannier functions : ",I4)')sic_wantran%nwan
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    if (wannier_lc) then
      write(151,'("  j : ",I4,"    n : ",I4,"    !  at ",3G18.10)')j,&
        sic_wantran%iwan(j),wanpos(:,n)
    else
      ias=wan_info(1,n)
      lm=wan_info(2,n)
      ispn=wan_info(3,n)
      write(c1,'(I6)')ias2ia(ias)
      write(c2,'(I1)')lm2m(lm)+lm2l(lm)+1
      c3=trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
        orbc(lm2l(lm)+1)//trim(adjustl(c2))
      if (spinpol) c3=trim(adjustl(c3))//"-"//spinc(ispn)
      write(151,'("  j : ",I4,"    n : ",I4,"    ! ",A,"  at ",3G18.10)')j,&
        sic_wantran%iwan(j),trim(c3),wanpos(:,n)
    endif
  enddo
  write(151,*)
  write(151,'("cutoff radius for Wannier functions   : ",F12.6)')sic_wan_cutoff
  write(151,'("cutoff radius for SIC matrix elements : ",F12.6)')sic_me_cutoff
  write(151,'("number of Wannier transitions         : ",I6)')sic_wantran%nwt
  write(151,'("number of translations for Bloch sums : ",I4)')sic_orbitals%ntr
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
!allocate(laplsv(nstsv))
!ekin_=0.d0
!do ikloc=1,nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  call genlapl(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
!    wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc),laplsv)
!  do j=1,nstsv
!    ekin_=ekin_+0.5d0*laplsv(j)*occsvnr(j,ik)
!  enddo
!enddo
!call mpi_grid_reduce(ekin_,dims=(/dim_k/))
!ekin_=ekin_/nkptnr
!if (wproc) then
!  write(151,*)
!  etot_=ekin_+engykncr+0.5d0*engyvcl+engyx+engyc+engymad-sic_energy_pot
!  write(151,'("kinetic energy : ",G18.10)')ekin_+engykncr
!  write(151,'("potential      : ",G18.10)')0.5d0*engyvcl+engyx+engyc
!  write(151,'("madelung       : ",G18.10)')engymad
!  write(151,'("sic correction : ",G18.10)')sic_energy_pot
!  write(151,'("total energy   : ",G18.10)')etot_
!  open(152,file="SIC_ETOT.OUT",form="FORMATTED",status="REPLACE")
!  write(152,'(G18.10)')etot_
!  close(152)
!endif

! get maximum number of G-vectors (this is required for the faster generation of 
!  Wannier functions in the interstitial region)
s_ngvec=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do ig=1,ngknr(ikloc)
    s_ngvec=max(s_ngvec,igkignr(ig,ikloc))
  enddo
enddo
! init Madness related variables 
#ifdef _MAD_
if (allocated(m_ngknr)) deallocate(m_ngknr)
allocate(m_ngknr(nkptnr))
m_ngknr=0
if (allocated(m_igkignr)) deallocate(m_igkignr)
allocate(m_igkignr(ngkmax,nkptnr))
m_igkignr=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  m_ngknr(ik)=ngknr(ikloc)
  m_igkignr(:,ik)=igkignr(:,ikloc)
enddo
call mpi_grid_reduce(m_ngknr(1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(m_igkignr(1,1),ngkmax*nkptnr,dims=(/dim_k/),all=.true.)
m_ngvec=0
do ik=1,nkptnr
  do ig=1,m_ngknr(ik)
    m_ngvec=max(m_ngvec,m_igkignr(ig,ik))
  enddo
enddo
if (allocated(m_wann_unkmt)) deallocate(m_wann_unkmt)
allocate(m_wann_unkmt(lmmaxvr,nufrmax,natmtot,nspinor,nkptnr))
if (allocated(m_wann_unkit)) deallocate(m_wann_unkit)
allocate(m_wann_unkit(ngkmax,nspinor,nkptnr))
call madness_init_box
#endif
! generate Wannier functions and corresponding potential
call sic_wan(151)
call timer_start(t_sic_me,reset=.true.)
if (wproc) then
  write(151,*)
  write(151,'(80("="))')
  write(151,'("matrix elements")')
  write(151,'(80("="))')
endif
! save old matrix elements
allocate(vwanme_old(sic_wantran%nwt))
vwanme_old=vwanme
! compute matrix elements of SIC potential
!  vwanme = <(W*V)_n|W_{n1,T}>
vwanme=zzero
nwtloc=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/))
do iloc=1,nwtloc
  i=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/),loc=iloc)
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  do ispn=1,nspinor
    vwanme(i)=vwanme(i)+s_dot_ll(pos1,pos2,s_wvlm(1,1,ispn,j),s_wanlm(1,1,ispn,j1))
  enddo
enddo
call mpi_grid_reduce(vwanme(1),sic_wantran%nwt,all=.true.)
! check localization criterion 
t2=-1.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  t1=abs(vwanme(i)-dconjg(vwanme(j)))
  if (t1.ge.t2) then
    t2=t1
    i1=i
    j1=j
  endif
enddo
! symmetrize the potential matrix elements
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  z1=0.5d0*(vwanme(i)+dconjg(vwanme(j)))
  vwanme(i)=z1
  vwanme(j)=dconjg(z1)
enddo
! compute RMS difference
t3=0.d0
do i=1,sic_wantran%nwt
  t3=t3+abs(vwanme(i)-vwanme_old(i))**2
enddo
deallocate(vwanme_old)
call timer_stop(t_sic_me)
if (wproc) then
  write(151,*)
  write(151,'("maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,'("matrix elements with maximum difference : ",2I6)')i1,j1
  write(151,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,i1),dreal(vwanme(i1)),dimag(vwanme(i1))
  write(151,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,j1),dreal(vwanme(j1)),dimag(vwanme(j1))
  write(151,*)
  write(151,'("diagonal matrix elements (<(W*V)_n|W_n>) :")')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    i=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'("  n : ",I4,8X,2G18.10)')n,dreal(vwanme(i)),dimag(vwanme(i))
  enddo  
  t3=sqrt(t3/sic_wantran%nwt)
  write(151,*)
  write(151,'("matrix elements RMS difference :",G18.10)')t3
  write(151,*)
  write(151,'("done in : ",F8.3," sec.")')timer_get_value(t_sic_me)
  write(151,*)
  call flushifc(151)
endif
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
! write to HDF5 file
call sic_writevwan
if (wproc) then
  call timestamp(151,"Done.")
  close(151)
endif
return
end
