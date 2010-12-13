subroutine sic_genfvprj
use modmain
use mod_sic
implicit none
integer ngk1
integer, allocatable :: igkig1(:)
real(8), allocatable :: gkc1(:)
real(8), allocatable :: tpgkc1(:,:)
real(8), allocatable :: vgkl1(:,:)
real(8), allocatable :: vgkc1(:,:)
complex(8), allocatable :: sfacgk1(:,:)
complex(8), allocatable :: evecfv1(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfir_(:)
complex(8), allocatable :: a(:,:,:)
complex(8), allocatable :: b(:,:,:)
complex(8), allocatable :: expikr(:)
integer ik,h,ikloc,ig,ir,n,ispn,istfv,istsv,ias,j,i
complex(8), allocatable :: wffvmt(:,:,:,:)


call timer_start(t_sic_genfvprj)
! on first SIC iteration Wannier functions are generated from LDA Hamiltonian
!  so we can compute overlap between Wannier states and first-variational
!  states analytically
if (.not.tsic_wv) then
  do ikloc=1,nkptloc
    sic_wb(:,:,:,ikloc)=zzero
    sic_wvb(:,:,:,ikloc)=zzero
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        do istfv=1,nstfv
          istsv=istfv+(ispn-1)*nstfv
          do i=1,nstsv
! TODO: zgemm?
            sic_wb(j,istfv,ispn,ikloc)=sic_wb(j,istfv,ispn,ikloc)+&
              dconjg(wann_c(n,i,ikloc)*evecsvloc(istsv,i,ikloc))
          enddo !j
        enddo !ispn
      enddo !i
    enddo !n
  enddo !ikloc  
  call timer_stop(t_sic_genfvprj)
  return
endif

allocate(igkig1(ngkmax))
allocate(gkc1(ngkmax))
allocate(tpgkc1(2,ngkmax))
allocate(vgkl1(3,ngkmax))
allocate(vgkc1(3,ngkmax))
allocate(sfacgk1(ngkmax,natmtot))
allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfir(ngrtot))
allocate(wfmt_(lmmaxvr,nmtloc))
allocate(wfir_(ngrloc))
allocate(a(sic_wantran%nwan,nstfv,nspinor))
allocate(b(sic_wantran%nwan,nstfv,nspinor))
allocate(expikr(ngrloc))
allocate(wffvmt(nstfv,lmmaxvr,nufrmax,natmtot))

do ik=1,nkpt
  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
  if (mpi_grid_x(dim_k).eq.h) evecfv1(:,:,:)=evecfvloc(:,:,:,ikloc)
  call mpi_grid_bcast(evecfv1(1,1,1),nmatmax*nstfv*nspnfv,dims=(/dim_k,dim2/),&
    root=(/h,0,0/))
  call gengpvec(vkl(1,ik),vkc(1,ik),ngk1,igkig1,vgkl1,vgkc1,gkc1,tpgkc1)
  call gensfacgp(ngk1,vgkc1,ngkmax,sfacgk1)
  call match(ngk1,gkc1,tpgkc1,sfacgk1,apwalm)
  call genwffvmt(lmaxvr,lmmaxvr,ngk1,evecfv1,apwalm,wffvmt)
  do ir=1,ngrloc
    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir+groffs)))
  enddo
  a=zzero
  b=zzero
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
  do istfv=1,nstfv
    wfir=zzero
    call timer_start(t_sic_genfvprj_wfmt)
    call sic_wavefmt(wffvmt,istfv,wfmt_)
    call timer_stop(t_sic_genfvprj_wfmt)
    call timer_start(t_sic_genfvprj_wfir)
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig1(ig)))=evecfv1(ig,istfv,1)
    enddo
    call zfftifc(3,ngrid,1,wfir)
    call timer_stop(t_sic_genfvprj_wfir)
    call timer_start(t_sic_genfvprj_dotp)
    call sic_copy_ir_z(.true.,wfir,wfir_)
    do ir=1,ngrloc
      wfir_(ir)=wfir_(ir)*expikr(ir)/sqrt(omega)
    enddo
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        a(j,istfv,ispn)=sic_dot_lb(.false.,vkc(1,ik),sic_orbitals%wanmt(1,1,1,ispn,j),&
          sic_orbitals%wanir(1,1,ispn,j),sic_orbitals%twanmt(1,1,n),wfmt_,wfir_)
        b(j,istfv,ispn)=sic_dot_lb(.false.,vkc(1,ik),sic_orbitals%wvmt(1,1,1,ispn,j),&
          sic_orbitals%wvir(1,1,ispn,j),sic_orbitals%twanmt(1,1,n),wfmt_,wfir_)
      enddo
    enddo !j
    call timer_stop(t_sic_genfvprj_dotp)
  enddo !istfv
  call timer_start(t_sic_genfvprj_dotp)
  call mpi_grid_reduce(a(1,1,1),sic_wantran%nwan*nstfv*nspinor,root=(/h,0,0/))
  call mpi_grid_reduce(b(1,1,1),sic_wantran%nwan*nstfv*nspinor,root=(/h,0,0/))
  call timer_stop(t_sic_genfvprj_dotp)
  if (mpi_grid_x(dim_k).eq.h) then
    sic_wb(:,:,:,ikloc)=a(:,:,:)
    sic_wvb(:,:,:,ikloc)=b(:,:,:)
  endif
enddo !ik
deallocate(igkig1,gkc1,tpgkc1,vgkl1,vgkc1,sfacgk1,evecfv1)
deallocate(apwalm,wfir,wfmt_,wfir_,a,b,expikr,wffvmt)
call timer_stop(t_sic_genfvprj)
return
end
