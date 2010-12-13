subroutine sic_genvhart(vhwanmt,vhwanir)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
use mod_linresp
use mod_sic
implicit none
complex(8), intent(out) :: vhwanmt(lmmaxvr,nmtloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan)
complex(8), intent(out) :: vhwanir(ngrloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan)
integer nvqloc,iqloc,it,iq,n,ig,ias,i,j
complex(8), allocatable :: pwmt(:,:)
complex(8), allocatable :: pwir(:)
complex(8), allocatable ::megqwan1(:,:,:)
complex(8) expikt,zt1
character*100 qnm

vhwanmt=zzero
vhwanir=zzero
call init_qbz(tq0bz,1)
call init_q_gq
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq
    call getqdir(iq,vqm(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif
call mpi_grid_barrier()
wannier_megq=.true.
megq_include_bands(:)=(/100.1d0,-100.1d0/)
call deletewantran(megqwantran)
call genwantran(megqwantran,-0.d0,0.01d0,diagwt=.true.)
allocate(megqwan1(sic_wantran%nwan,ngqmax,nvq))
megqwan1=zzero
! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
call timer_start(10,reset=.true.)
! loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.false.)
! save <n,T=0|e^{-i(G+q)r}|n,T=0>
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    megqwan1(j,1:ngq(iq),iq)=megqwan(megqwantran%iwtidx(n,n,0,0,0),1:ngq(iq))
  enddo
enddo
call mpi_grid_reduce(megqwan1(1,1,1),sic_wantran%nwan*ngqmax*nvq,&
  dims=(/dim_q/),all=.true.)
call timer_stop(10)
! allocate arrays for plane-wave
allocate(pwmt(lmmaxvr,nmtloc))
allocate(pwir(ngrloc))
! generate Hartree potential
call timer_start(11,reset=.true.)
do iq=1,nvq
  do ig=1,ngq(iq)
    call sic_genpw(vgqc(1,ig,iq),pwmt,pwir)
    do it=1,sic_orbitals%ntr
      expikt=exp(zi*dot_product(sic_orbitals%vtc(:,it),vqc(:,iq)))/nkptnr/omega
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        if (sic_apply(n).eq.2) then
          zt1=megqwan1(j,ig,iq)*vhgq(ig,iq)*expikt
          do i=1,nmtloc
            ias=(mtoffs+i-1)/nrmtmax+1
            if (sic_orbitals%twanmt(ias,it,n)) then       
              call zaxpy(lmmaxvr,zt1,pwmt(1,i),1,vhwanmt(1,i,it,1,j),1)
            endif
          enddo
          call zaxpy(ngrloc,zt1,pwir,1,vhwanir(1,it,1,j),1)
        endif !sic_apply(n).eq.2
      enddo !n
    enddo !it
  enddo !ig
enddo !iq
call timer_stop(11)
deallocate(pwmt,pwir,megqwan1)
return
end
