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
integer ik,h,ikloc,ig,ir,n,ispn,ist,ias,i1,i2,i3
real(8) v2(3),v3(3)
logical exist

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(igkig1(ngkmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))
allocate(wfmt_(lmmaxvr,nmtloc))
allocate(wfir_(ngrloc))
allocate(a(nwann,nstfv,nspinor))
allocate(b(nwann,nstfv,nspinor))
allocate(expikr(ngrloc))

do ik=1,nkpt
  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
  if (mpi_grid_x(dim_k).eq.h) then
    evecfv1(:,:,:)=evecfvloc(:,:,:,ikloc)
    igkig1(:)=igkig(:,1,ikloc)
! get apw coeffs 
    call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
      apwalm)
  endif  
  call mpi_grid_bcast(evecfv1(1,1,1),nmatmax*nstfv*nspnfv,dims=(/dim_k,dim2/),&
    root=(/h,0/))
  call mpi_grid_bcast(igkig1(1),ngkmax,dims=(/dim_k,dim2/),root=(/h,0/))
  call mpi_grid_bcast(apwalm(1,1,1,1),ngkmax*apwordmax*lmmaxapw*natmtot,&
    dims=(/dim_k,dim2/),root=(/h,0/))
  do ir=1,ngrloc
    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir+groffs)))
  enddo
  a=zzero
  b=zzero
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
  do ist=1,nstfv
    wfmt=zzero
    wfir=zzero
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
        evecfv1(1,ist,1),lmmaxvr,wfmt(1,1,ias))
    enddo
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig1(ig)))=evecfv1(ig,ist,1)
    enddo
    call zfftifc(3,ngrid,1,wfir)
    call sic_copy_mt_z(.true.,lmmaxvr,wfmt,wfmt_)
    call sic_copy_ir_z(.true.,wfir,wfir_)
    do ir=1,ngrloc
      wfir_(ir)=wfir_(ir)*expikr(ir)/sqrt(omega)
    enddo
    do n=1,nwann
      do ispn=1,nspinor
        a(n,ist,ispn)=sic_dot_lb(vkc(1,ik),wanmt(1,1,1,ispn,n),&
          wanir(1,1,ispn,n),twanmt(1,1,n),wfmt_,wfir_)
        b(n,ist,ispn)=sic_dot_lb(vkc(1,ik),wvmt(1,1,1,ispn,n),&
          wvir(1,1,ispn,n),twanmt(1,1,n),wfmt_,wfir_)
      enddo
    enddo !n
  enddo !ist
  if (mpi_grid_x(dim_k).eq.h) then
    sic_wb(:,:,:,ikloc)=a(:,:,:)
    sic_wvb(:,:,:,ikloc)=b(:,:,:)
  endif
enddo !ik
deallocate(evecfv1,igkig1,apwalm,wfmt,wfir,wfmt_,wfir_,a,b)
deallocate(expikr)
return
end