subroutine sic_genwan
use modmain
use mod_nrkp
use mod_sic
implicit none
integer is,ic,ias,ir,ig,ikloc,ik,l
integer io,lm,n,ispn,h,nloc,it
real(8) v1(3),vtrc(3)
complex(8), allocatable :: zfir(:,:,:)
complex(8), allocatable :: wmt(:,:,:,:)
complex(8), allocatable :: wir(:,:)
complex(8), allocatable :: expikr(:,:)
complex(8) expikt,zt1
allocate(wmt(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(wir(ngrtot,nspinor))
allocate(zfir(ngrtot,nspinor,nkptnrloc))
allocate(expikr(ngrtot,nkptnrloc))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do ir=1,ngrtot
    expikr(ir,ikloc)=exp(zi*dot_product(vkcnr(:,ik),vgrc(:,ir)))
  enddo
enddo
wanmt=zzero
wanir=zzero
do n=1,nwann
! precompute interstitial wave-functions
  zfir=zzero
  do ikloc=1,nkptnrloc
    do ispn=1,nspinor
      do ig=1,ngknr(ikloc)
        zfir(igfft(igkignr(ig,ikloc)),ispn,ikloc)=&
          wann_unkit(ig,ispn,n,ikloc)
      end do
      call zfftifc(3,ngrid,1,zfir(1,ispn,ikloc))
    end do !ispn
  end do  
  do it=1,ntr
    wmt=zzero
    wir=zzero    
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      expikt=exp(zi*dot_product(vkcnr(:,ik),vtc(:,it)))/nkptnr
! muffin-tin part    
      call timer_start(1)
      do ias=1,natmtot
        if (twanmt(ias,it,n)) then
          is=ias2is(ias)
          ic=ias2ic(ias)
          do ispn=1,nspinor
            do ir=1,nrmt(is)
              do l=0,lmaxvr
                do io=1,nufr(l,is)
                  zt1=expikt*ufr(ir,l,io,ic)
                  do lm=l**2+1,(l+1)**2
                    wmt(lm,ir,ias,ispn)=wmt(lm,ir,ias,ispn)+&
                      zt1*wann_unkmt(lm,io,ias,ispn,n,ikloc)
                  end do !lm
                end do !io
              end do !l
            end do !ir
          end do !ispn
        end if !twanmt
      end do !ias
      call timer_stop(1)
      call timer_start(2)
! interstitial part      
      do ir=1,ngrtot
        zt1=expikt*expikr(ir,ikloc)/sqrt(omega)
        do ispn=1,nspinor
          wir(ir,ispn)=wir(ir,ispn)+zt1*zfir(ir,ispn,ikloc)
        end do !ispn      
      end do !ir
      call timer_stop(2)
    end do !ikloc
    if (mpi_grid_side(dims=(/dim_k/))) then
      call mpi_grid_reduce(wmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor,&
        dims=(/dim_k/),all=.true.)
      call mpi_grid_reduce(wir(1,1),ngrtot*nspinor,dims=(/dim_k/),&
        all=.true.)
    end if
    call mpi_grid_bcast(wmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor,&
      dims=ortdims((/dim_k/)))
    call mpi_grid_bcast(wir(1,1),ngrtot*nspinor,dims=ortdims((/dim_k/)))
    do ispn=1,nspinor
      call sic_copy_mt_z(.true.,lmmaxvr,wmt(1,1,1,ispn),&
        wanmt(1,1,it,ispn,n))
      call sic_copy_ir_z(.true.,wir(1,ispn),wanir(1,it,ispn,n))
    end do
  end do !it
end do !n
deallocate(wmt,wir,zfir,expikr)
! cutoff in the interstitial
do it=1,ntr
  do n=1,nwann
    do ir=1,ngrloc
      v1(:)=vgrc(:,ir+groffs)+vtc(:,it)-&
        atposc(:,ias2ia(iwann(1,n)),ias2is(iwann(1,n)))
      if (sqrt(sum(v1(:)**2)).gt.wann_r_cutoff) then
        wanir(ir,it,:,n)=zzero
      endif
    enddo
  enddo
end do !it
return
end

