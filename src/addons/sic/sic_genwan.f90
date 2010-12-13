subroutine sic_genwan
use modmain
use mod_nrkp
use mod_sic
implicit none
integer is,ic,ias,ir,ig,ikloc,ik,l
integer io,lm,n,ispn,it,j
real(8) v1(3)
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
sic_orbitals%wanmt=zzero
sic_orbitals%wanir=zzero
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
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
  do it=1,sic_orbitals%ntr
    wmt=zzero
    wir=zzero    
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      expikt=exp(zi*dot_product(vkcnr(:,ik),sic_orbitals%vtc(:,it)))/nkptnr
! muffin-tin part    
      call timer_start(1)
      do ias=1,natmtot
        if (sic_orbitals%twanmt(ias,it,n)) then
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
        end if !sic_orbitals%twanmt
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
        dims=(/dim_k/))
      call mpi_grid_reduce(wir(1,1),ngrtot*nspinor,dims=(/dim_k/))
    end if
    call mpi_grid_bcast(wmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor)
    call mpi_grid_bcast(wir(1,1),ngrtot*nspinor)
    do ispn=1,nspinor
      call sic_copy_mt_z(.true.,lmmaxvr,wmt(1,1,1,ispn),&
        sic_orbitals%wanmt(1,1,it,ispn,j))
      call sic_copy_ir_z(.true.,wir(1,ispn),sic_orbitals%wanir(1,it,ispn,j))
    end do
  end do !it
end do !j
deallocate(wmt,wir,zfir,expikr)
! cutoff in the interstitial
do it=1,sic_orbitals%ntr
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    do ir=1,ngrloc
      v1(:)=vgrc(:,ir+groffs)+sic_orbitals%vtc(:,it)-&
        atposc(:,ias2ia(wan_info(1,n)),ias2is(wan_info(1,n)))
      if (sqrt(sum(v1(:)**2)).gt.sic_wan_cutoff) then
        sic_orbitals%wanir(ir,it,:,j)=zzero
      endif
    enddo
  enddo
end do !it
return
end

