subroutine sic_genfvprj
use modmain
use mod_lf
implicit none
complex(8), allocatable :: evecfv1(:,:,:)
integer, allocatable :: igkig1(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: a(:,:,:)
complex(8), allocatable :: b(:,:,:)
integer ik,h,ikloc,ig,ir,n,nloc,ispn,ist,ias,i1,i2,i3
real(8) v2(3),v3(3)

allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(igkig1(ngkmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))
allocate(a(nwann,nstfv,nspinor))
allocate(b(nwann,nstfv,nspinor))

do ik=1,nkpt
  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
  if (mpi_grid_x(dim_k).eq.h) then
    evecfv1(:,:,:)=evecfvloc(:,:,:,ikloc)
    igkig1(:)=igkig(:,1,ikloc)
! get apw coeffs 
    call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
      apwalm)
  endif  
  call mpi_grid_bcast(evecfv1(1,1,1),nmatmax*nstfv*nspnfv,dims=(/dim_k,dim_t/),&
    root=(/h,0/))
  call mpi_grid_bcast(igkig1(1),ngkmax,dims=(/dim_k,dim_t/),root=(/h,0/))
  call mpi_grid_bcast(apwalm(1,1,1,1),ngkmax*apwordmax*lmmaxapw*natmtot,&
    dims=(/dim_k,dim_t/),root=(/h,0/))
  a=zzero
  b=zzero
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
  do ist=1,nstfv
    wfmt=zzero
    wfir=zzero
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
        evecfv1(:,ist,1),lmmaxvr,wfmt(1,1,ias))
    enddo
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig1(ig)))=evecfv1(ig,ist,1)
    enddo
    call zfftifc(3,ngrid,1,wfir)
    wfir(:)=wfir(:)/sqrt(omega)
    ir=0
    do i3=0,ngrid(3)-1
      v2(3)=dble(i3)/dble(ngrid(3))
      do i2=0,ngrid(2)-1
        v2(2)=dble(i2)/dble(ngrid(2))
        do i1=0,ngrid(1)-1
          v2(1)=dble(i1)/dble(ngrid(1))
          ir=ir+1
          call r3mv(avec,v2,v3)
          wfir(ir)=exp(zi*dot_product(vkc(:,ik),v3(:)))*wfir(ir)
        enddo
      enddo
    enddo
    do nloc=1,nwannloc
      n=mpi_grid_map(nwann,dim_k,loc=nloc)
      do ispn=1,nspinor
        a(n,ist,ispn)=lf_dot_blh(.true.,vkc(1,ik),wanmt(1,1,1,1,ispn,nloc),&
          wanir(1,1,ispn,nloc),wfmt,wfir)
        b(n,ist,ispn)=lf_dot_blh(.true.,vkc(1,ik),wvmt(1,1,1,1,ispn,nloc),&
          wvir(1,1,ispn,nloc),wfmt,wfir)
      enddo
    enddo !nloc
  enddo !ist
  call mpi_grid_reduce(a(1,1,1),nwann*nstfv*nspinor,dims=(/dim_k/),root=(/h/))
  call mpi_grid_reduce(b(1,1,1),nwann*nstfv*nspinor,dims=(/dim_k/),root=(/h/))
  if (mpi_grid_x(dim_k).eq.h) then
    sic_wb(:,:,:,ikloc)=a(:,:,:)
    sic_wvb(:,:,:,ikloc)=b(:,:,:)
  endif
enddo !ik
deallocate(evecfv1,igkig1,apwalm,wfmt,wfir,a,b)
return
end