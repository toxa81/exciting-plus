subroutine sic_genfvprj
use modmain
use mod_sic
implicit none
complex(8), allocatable :: evecfv1(:,:,:)
integer, allocatable :: igkig1(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfir_(:)
complex(8), allocatable :: a(:,:,:)
complex(8), allocatable :: b(:,:,:)
complex(8), allocatable :: expikr(:)
integer ik,h,ikloc,ig,ir,n,ispn,istfv,istsv,ias,j,i
real(8), allocatable :: gkc1(:)
real(8), allocatable :: tpgkc1(:,:)
complex(8), allocatable :: sfacgk1(:,:)

call timer_start(t_sic_genfvprj)
goto 10
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

10 continue
allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(igkig1(ngkmax))
allocate(gkc1(ngkmax))
allocate(tpgkc1(2,ngkmax))
allocate(sfacgk1(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))
allocate(wfmt_(lmmaxvr,nmtloc))
allocate(wfir_(ngrloc))
allocate(a(sic_wantran%nwan,nstfv,nspinor))
allocate(b(sic_wantran%nwan,nstfv,nspinor))
allocate(expikr(ngrloc))

do ik=1,nkpt
  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
  if (mpi_grid_x(dim_k).eq.h) then
    evecfv1(:,:,:)=evecfvloc(:,:,:,ikloc)
    igkig1(:)=igkig(:,1,ikloc)
    gkc1(:)=gkc(:,1,ikloc)
    tpgkc1(:,:)=tpgkc(:,:,1,ikloc)
    sfacgk1(:,:)=sfacgk(:,:,1,ikloc)
! get apw coeffs 
!    call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
!      apwalm)
  endif  
  call mpi_grid_bcast(evecfv1(1,1,1),nmatmax*nstfv*nspnfv,dims=(/dim_k,dim2/),&
    root=(/h,0/))
  call mpi_grid_bcast(igkig1(1),ngkmax,dims=(/dim_k,dim2/),root=(/h,0/))
  call mpi_grid_barrier
  call mpi_grid_bcast(gkc1(1),ngkmax,dims=(/dim_k,dim2/),root=(/h,0/))
  call mpi_grid_barrier
  call mpi_grid_bcast(tpgkc1(1,1),2*ngkmax,dims=(/dim_k,dim2/),root=(/h,0/))
  call mpi_grid_barrier
  call mpi_grid_bcast(sfacgk1(1,1),ngkmax*natmtot,dims=(/dim_k,dim2/),root=(/h,0/))
  call match(ngk(1,ik),gkc1,tpgkc1,sfacgk1,apwalm)
!  call mpi_grid_bcast(apwalm(1,1,1,1),ngkmax*apwordmax*lmmaxapw*natmtot,&
!    dims=(/dim_k,dim2/),root=(/h,0/))
  do ir=1,ngrloc
    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir+groffs)))
  enddo
  a=zzero
  b=zzero
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
  do istfv=1,nstfv
    wfmt=zzero
    wfir=zzero
    call timer_start(t_sic_genfvprj_wfmt)
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
        evecfv1(1,istfv,1),lmmaxvr,wfmt(1,1,ias))
    enddo
    call timer_stop(t_sic_genfvprj_wfmt)
    call timer_start(t_sic_genfvprj_wfir)
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig1(ig)))=evecfv1(ig,istfv,1)
    enddo
    call zfftifc(3,ngrid,1,wfir)
    call timer_stop(t_sic_genfvprj_wfir)
    call timer_start(t_sic_genfvprj_dotp)
    call sic_copy_mt_z(.true.,lmmaxvr,wfmt,wfmt_)
    call sic_copy_ir_z(.true.,wfir,wfir_)
    do ir=1,ngrloc
      wfir_(ir)=wfir_(ir)*expikr(ir)/sqrt(omega)
    enddo
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        a(j,istfv,ispn)=sic_dot_lb(vkc(1,ik),wanmt(1,1,1,ispn,n),&
          wanir(1,1,ispn,n),twanmt(1,1,n),wfmt_,wfir_,h)
        b(j,istfv,ispn)=sic_dot_lb(vkc(1,ik),wvmt(1,1,1,ispn,j),&
          wvir(1,1,ispn,j),twanmt(1,1,n),wfmt_,wfir_,h)
      enddo
    enddo !j
    call timer_stop(t_sic_genfvprj_dotp)
  enddo !istfv
  if (mpi_grid_x(dim_k).eq.h) then
    sic_wb(:,:,:,ikloc)=a(:,:,:)
    sic_wvb(:,:,:,ikloc)=b(:,:,:)
  endif
enddo !ik
deallocate(evecfv1,igkig1,apwalm,wfmt,wfir,wfmt_,wfir_,a,b)
deallocate(expikr,gkc1,tpgkc1,sfacgk1)
call timer_stop(t_sic_genfvprj)
return
end
