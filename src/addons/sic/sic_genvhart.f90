subroutine sic_genvhart
use modmain
use mod_addons_q
use mod_lf
implicit none
integer nvqloc,iqloc,iq,n,ngqloc,igloc,itr,itloc,ig
real(8) vtrc(3)
complex(8), allocatable ::megqwan1(:,:,:)
complex(8), allocatable :: pwmt(:,:,:)
complex(8), allocatable :: pwir(:)
complex(8) expikt
logical lgamma

lr_e1=100.1d0
lr_e2=-100.1d0
wannier_megq=.true.
lgamma=.true.
call init_qbz(lgamma,1)
call init_q_gq
! create q-directories
!if (mpi_grid_root()) then
!  call system("mkdir -p q")
!  do iq=1,nvq0
!    call getqdir(iq,ivq0m_list(:,iq),qnm)
!    call system("mkdir -p "//trim(qnm))
!  enddo
!endif
! distribute q-vectors along 3-rd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
allocate(megqwan1(nwann,ngqmax,nvq))
megqwan1=zzero
call timer_start(10,reset=.true.)
! loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.false.,.false.)
! save <n,T=0|e^{-i(G+q)r}|n,T=0>
  do n=1,nwann
    megqwan1(n,1:ngq(iq),iq)=megqwan(idxmegqwan(n,n,0,0,0),1:ngq(iq))
  enddo
enddo
call mpi_grid_reduce(megqwan1(1,1,1),nwann*ngqmax*nvq,dims=(/dim_q/), &
  all=.true.)
call timer_stop(10)
! allocate arrays for plane-wave
allocate(pwmt(lmmaxvr,nrmtmax,natmtot))
allocate(pwir(ngrtot))
! generate Hartree potential
call timer_start(11,reset=.true.)
do iq=1,nvq
  ngqloc=mpi_grid_map(ngq(iq),dim_k)
  do igloc=1,ngqloc
    ig=mpi_grid_map(ngq(iq),dim_k,loc=igloc)
    call genpw((/0,0,0/),vgqc(1,ig,iq),pwmt,pwir)
    do itloc=1,ntrloc
      itr=mpi_grid_map(ntr,dim_t,loc=itloc)
      vtrc(:)=vtl(1,itr)*avec(:,1)+vtl(2,itr)*avec(:,2)+vtl(3,itr)*avec(:,3)
      expikt=exp(zi*dot_product(vtrc(:),vqc(:,iq)))
      do n=1,nwann
        vwanmt(:,:,:,itloc,1,n)=vwanmt(:,:,:,itloc,1,n)+&
          megqwan1(n,ig,iq)*vhgq(ig,iq)*pwmt(:,:,:)*expikt
        vwanir(:,itloc,1,n)=vwanir(:,itloc,1,n)+&
          megqwan1(n,ig,iq)*vhgq(ig,iq)*pwir(:)*expikt
      enddo !n
    enddo !itloc
  enddo !igloc
enddo !iq
do n=1,nwann
  do itloc=1,ntrloc
    call mpi_grid_reduce(vwanmt(1,1,1,itloc,1,n),lmmaxvr*nrmtmax*natmtot,&
      dims=(/dim_k/),all=.true.)
    call mpi_grid_reduce(vwanir(1,itloc,1,n),ngrtot,dims=(/dim_k/),&
      all=.true.)
  enddo
enddo
vwanmt=vwanmt/nkptnr/omega
vwanir=vwanir/nkptnr/omega
call timer_stop(11)
deallocate(pwmt,pwir,megqwan1)
return
end