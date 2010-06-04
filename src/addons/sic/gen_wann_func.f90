subroutine gen_wann_func(vtrl,ngknr,vgkcnr,igkignr,wanmt_,wanir_)
use modmain
implicit none
integer, intent(in) :: vtrl(3)
integer, intent(in) :: ngknr(nkptnrloc)
real(8), intent(in) :: vgkcnr(3,ngkmax,nkptnrloc)
integer, intent(in) :: igkignr(ngkmax,nkptnrloc)
complex(8), intent(out) :: wanmt_(lmmaxvr,nrmtmax,natmtot,nspinor,nwann)
complex(8), intent(out) :: wanir_(ngrtot,nspinor,nwann)
integer ia,is,ias,ir,ir0,i1,i2,i3,ig,ikloc,ik
integer io,lm,n,ispn,itmp(3)
real(8) v2(3),v3(3),r0,vr0(3),vtrc(3)
complex(8), allocatable :: zfir(:,:,:)
complex(8) expikr
logical, external :: vrinmt
wanmt_=zzero
wanir_=zzero
vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
! muffin-tin part
call timer_start(1)
do ias=1,natmtot
  is=ias2is(ias)
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    expikr=exp(zi*dot_product(vkcnr(:,ik),vtrc(:)))/nkptnr
    do ir=1,nrmt(is)
      do lm=1,lmmaxvr
        do io=1,nufr(lm2l(lm),is)
          do ispn=1,nspinor
            do n=1,nwann
              wanmt_(lm,ir,ias,ispn,n)=wanmt_(lm,ir,ias,ispn,n)+&
                expikr*wann_unkmt(lm,io,ias,ispn,n,ikloc)*&
                ufr(ir,lm2l(lm),io,ias)
            enddo !n
          enddo !ispn
        enddo !io
      enddo !lm
    enddo !ir
  enddo !ikloc
enddo !ias 
call timer_stop(1)
! interstitial part
call timer_start(2)
allocate(zfir(ngrtot,nspinor,nwann))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  zfir=zzero
  do n=1,nwann
    do ispn=1,nspinor
      do ig=1,ngknr(ikloc)
        zfir(igfft(igkignr(ig,ikloc)),ispn,n)=wann_unkit(ig,ispn,n,ikloc)
      enddo
      call zfftifc(3,ngrid,1,zfir(:,ispn,n))
    enddo !ispn
  enddo !n
  ir=0
  do i3=0,ngrid(3)-1
    v2(3)=dble(i3)/dble(ngrid(3))
    do i2=0,ngrid(2)-1
      v2(2)=dble(i2)/dble(ngrid(2))
      do i1=0,ngrid(1)-1
        v2(1)=dble(i1)/dble(ngrid(1))
        ir=ir+1
        call r3mv(avec,v2,v3)
        v3(:)=v3(:)+vtrc(:)
        expikr=exp(zi*dot_product(vkcnr(:,ik),v3(:)))
        do ispn=1,nspinor
          do n=1,nwann
            wanir_(ir,ispn,n)=wanir_(ir,ispn,n)+expikr*zfir(ir,ispn,n)
          enddo !in
        enddo
      enddo
    enddo
  enddo
enddo
wanir_(:,:,:)=wanir_(:,:,:)/sqrt(omega)/nkptnr
call timer_stop(2)
call mpi_grid_reduce(wanmt_(1,1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor*nwann,&
  dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(wanir_(1,1,1),ngrtot*nspinor*nwann,dims=(/dim_k/),all=.true.)  
deallocate(zfir)
return
end

