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
complex(8), allocatable :: wffvmt(:,:,:,:)
integer ik,h,ikloc,ig,ir,n,ispn,ist,ias,i,j,m1,l,lm1,lm,is,io1,ic
logical exist
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)

if (.not.tsic_wv) return
!if (.not.tsic_wv) then
!  do ikloc=1,nkptloc
!    sic_wb(:,:,:,ikloc)=0.8d0*sic_wb(:,:,:,ikloc) !zzero
!    sic_wvb(:,:,:,ikloc)=zzero
!    do n=1,nwantot
!      do ispn=1,nspinor
!        do ist=1,nstfv
!          i=ist+(ispn-1)*nstfv
!          do j=1,nstsv
!! TODO: zgemm?
!            sic_wb(n,ist,ispn,ikloc)=sic_wb(n,ist,ispn,ikloc)+&
!              0.2d0*dconjg(wann_c(n,j,ikloc)*evecsvloc(i,j,ikloc))
!          enddo !j
!        enddo !ispn
!      enddo !i
!    enddo !n
!  enddo !ikloc  
!  return
!endif

allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(igkig1(ngkmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))
allocate(wfmt_(lmmaxvr,nmtloc))
allocate(wfir_(ngrloc))
allocate(a(nwantot,nstfv,nspinor))
allocate(b(nwantot,nstfv,nspinor))
allocate(expikr(ngrloc))
!allocate(wffvmt(nstfv,lmmaxvr,nufrmax,natmtot))

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
!  if (.not.tsic_wv) then
!    call genwffvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv1,apwalm,wffvmt)
!  endif
  a=zzero
  b=zzero
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
  do ist=1,nstfv
    wfmt=zzero
    wfir=zzero
    if (tsic_wv) then
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
      do n=1,nwantot
        do ispn=1,nspinor
          a(n,ist,ispn)=sic_dot_lb(vkc(1,ik),wanmt(1,1,1,ispn,n),&
            wanir(1,1,ispn,n),twanmt(1,1,n),wfmt_,wfir_)
          b(n,ist,ispn)=sic_dot_lb(vkc(1,ik),wvmt(1,1,1,ispn,n),&
            wvir(1,1,ispn,n),twanmt(1,1,n),wfmt_,wfir_)
        enddo
      enddo !n
    else
!      do n=1,nwantot
!        ias=wan_info(1,n) 
!        lm=wan_info(2,n)
!        ispn=wan_info(3,n)
!        l=lm2l(lm)
!        is=ias2is(ias)
!        ic=ias2ic(ias)
!        do io1=1,nufr(l,is)
!          do ir=1,nrmt(is)
!            fr(ir)=ufr(ir,l,io1,ic)*(1+cos(pi*spr(ir,is)/rmt(is)))*(spr(ir,is)**2)                                                        
!          enddo
!          call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
!          do m1=-l,l
!            lm1=idxlm(l,m1)
!            a(n,ist,ispn)=a(n,ist,ispn)+dconjg(wffvmt(ist,lm1,io1,ias))*&
!              gr(nrmt(is))*rylm_lps(lm,lm1,ias)
!          enddo
!        enddo
!      enddo
    endif
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
