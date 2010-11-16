subroutine sic_genwann(vtrl,ngknr,igkignr,wanmt_,wanir_)
use modmain
use mod_lf
implicit none
integer, intent(in) :: vtrl(3)
integer, intent(in) :: ngknr(nkptnrloc)
integer, intent(in) :: igkignr(ngkmax,nkptnrloc)
complex(8), intent(out) :: wanmt_(lmmaxvr,nrmtmax,natmtot,nspinor,nwannloc)
complex(8), intent(out) :: wanir_(ngrtot,nspinor,nwannloc)
integer is,ias,ir,ig,ikloc,ik,l
integer io,lm,n,ispn,h,nloc
real(8) v1(3),vtrc(3)
complex(8), allocatable :: zfir(:,:)
complex(8), allocatable :: wmt(:,:,:,:)
complex(8), allocatable :: wir(:,:)
complex(8) expikr,zt1
logical, external :: vrinmt
wanmt_=zzero
wanir_=zzero
vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
allocate(wmt(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(wir(ngrtot,nspinor))
allocate(zfir(ngrtot,nspinor))
do n=1,nwann
  nloc=mpi_grid_map(nwann,dim_k,x=h,glob=n)
  wmt=zzero
  wir=zzero
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! muffin-tin part
    call timer_start(1)
    expikr=exp(zi*dot_product(vkcnr(:,ik),vtrc(:)))/nkptnr
    do ias=1,natmtot
      is=ias2is(ias)
      if (twanmt(ias,ivtit(vtrl(1),vtrl(2),vtrl(3)),n)) then
        do ispn=1,nspinor
          do ir=1,nrmt(is)
            do l=0,lmaxvr
              do io=1,nufr(l,is)
                zt1=expikr*ufr(ir,l,io,ias2ic(ias))
                do lm=l**2+1,(l+1)**2
                  wmt(lm,ir,ias,ispn)=wmt(lm,ir,ias,ispn)+&
                    zt1*wann_unkmt(lm,io,ias,ispn,n,ikloc)
                enddo !lm
              enddo !io
            enddo !l
          enddo !ir
        enddo !ispn
      endif
    enddo !ias 
    call timer_stop(1)
! interstitial part
    call timer_start(2)
    zfir=zzero
    do ispn=1,nspinor
      do ig=1,ngknr(ikloc)
        zfir(igfft(igkignr(ig,ikloc)),ispn)=wann_unkit(ig,ispn,n,ikloc)
      enddo
      call zfftifc(3,ngrid,1,zfir(:,ispn))
    enddo !ispn
    do ir=1,ngrtot
      v1(:)=vgrc(:,ir)+vtrc(:)
      zt1=exp(zi*dot_product(vkcnr(:,ik),v1(:)))
      do ispn=1,nspinor
        wir(ir,ispn)=wir(ir,ispn)+zt1*zfir(ir,ispn)
      enddo      
    enddo
    call timer_stop(2)
  enddo !ikloc
  call mpi_grid_reduce(wmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor,&
    dims=(/dim_k/),root=(/h/))
  call mpi_grid_reduce(wir(1,1),ngrtot*nspinor,dims=(/dim_k/),root=(/h/)) 
  wir=wir/sqrt(omega)/nkptnr
  if (mpi_grid_x(dim_k).eq.h) then
    wanmt_(:,:,:,:,nloc)=wmt(:,:,:,:)
    wanir_(:,:,nloc)=wir(:,:)
  endif
enddo !n
deallocate(wmt,wir,zfir)
! cutoff in the interstitial
do nloc=1,nwannloc
  n=mpi_grid_map(nwann,dim_k,loc=nloc)
  do ir=1,ngrtot
    v1(:)=vgrc(:,ir)+vtrc(:)-atposc(:,ias2ia(iwann(1,n)),ias2is(iwann(1,n)))
    if (sqrt(sum(v1(:)**2)).gt.wann_r_cutoff) then
      wanir_(ir,:,nloc)=zzero
    endif
  enddo
enddo
return
end

