subroutine sic_genvhart(vhwanmt,vhwanir)
use modmain
use mod_addons_q
use mod_lf
implicit none
complex(8), intent(out) :: vhwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwannloc)
complex(8), intent(out) :: vhwanir(ngrtot,ntrloc,nspinor,nwannloc)
integer nvqloc,iqloc,iq,n,nloc,itr,itloc,ig,l,m1,m2
integer ias,jas,ir,i1,i2,i3
real(8) vtrc(3),v2(3),v3(3)
complex(8), allocatable ::megqwan1(:,:,:)
complex(8), allocatable :: pwmt(:,:,:)
complex(8), allocatable :: pwir(:)
complex(8) expikt
!character*100 qnm

vhwanmt=zzero
vhwanir=zzero
wannier_megq=.true.
call init_qbz(tq0bz,1)
call init_q_gq
! create q-directories
!if (mpi_grid_root()) then
!  call system("mkdir -p q")
!  do iq=1,nvq
!    call getqdir(iq,vqm(:,iq),qnm)
!    call system("mkdir -p "//trim(qnm))
!  enddo
!endif
! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
allocate(megqwan1(nwann,ngqmax,nvq))
megqwan1=zzero
chi0_include_bands(:)=(/100.1d0,-100.1d0/)
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
  do ig=1,ngq(iq)
    call genpw(vgqc(1,ig,iq),pwmt,pwir)
    do itloc=1,ntrloc
      itr=mpi_grid_map(ntr,dim_t,loc=itloc)
      vtrc(:)=vtl(1,itr)*avec(:,1)+vtl(2,itr)*avec(:,2)+vtl(3,itr)*avec(:,3)
      expikt=exp(zi*dot_product(vtrc(:),vqc(:,iq)))/nkptnr/omega
      do nloc=1,nwannloc
        n=mpi_grid_map(nwann,dim_k,loc=nloc)
        do ias=1,natmtot
          if (twanmt(ias,itr,n)) then
            vhwanmt(:,:,ias,itloc,1,nloc)=vhwanmt(:,:,ias,itloc,1,nloc)+&
              megqwan1(n,ig,iq)*vhgq(ig,iq)*pwmt(:,:,ias)*expikt          
          endif
        enddo
        vhwanir(:,itloc,1,nloc)=vhwanir(:,itloc,1,nloc)+&
          megqwan1(n,ig,iq)*vhgq(ig,iq)*pwir(:)*expikt
      enddo !nloc
    enddo !itloc
  enddo !igloc
enddo !iq
call timer_stop(11)
deallocate(pwmt,pwir,megqwan1)
return
end