subroutine sic_convert_unk(ngp,igpig,unkmt_in,unkit_in,unkmt_out,unkit_out)
use modmain
use mod_sic
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: unkmt_in(lmmaxapw,nufrmax,natmtot,nspinor)
complex(8), intent(in) :: unkit_in(ngkmax,nspinor)
complex(8), intent(inout) :: unkmt_out(nrmtmax,lmmaxapw,natmtot,nspinor)
complex(8), intent(inout) :: unkit_out(ngkmax,nspinor)
!
integer ispn,ias,is,ic,lm,l,io,ig
complex(8), allocatable :: zv(:)
complex(8), allocatable :: unkmt_old(:,:,:,:)
complex(8), allocatable :: unkit_old(:,:)
!

allocate(unkmt_old(nrmtmax,lmmaxapw,natmtot,nspinor))
allocate(unkit_old(ngkmax,nspinor))
unkmt_old=unkmt_out
unkit_old=unkit_out
unkmt_out=zzero
unkit_out=zzero

allocate(zv(ngrtot))
do ispn=1,nspinor
! muffin-tin part
  do ias=1,natmtot
    is=ias2is(ias)
    ic=ias2ic(ias)
    do lm=1,lmmaxapw
      l=lm2l(lm)
      do io=1,nufr(l,is)
        unkmt_out(:,lm,ias,ispn)=unkmt_out(:,lm,ias,ispn)+&
              &unkmt_in(lm,io,ias,ispn)*ufr(:,l,io,ic)
      enddo !io
    enddo !lm
  enddo !ias
! interstitial
  zv=zzero
  do ig=1,ngp
    zv(igfft(igpig(ig)))=unkit_in(ig,ispn)
  enddo
  call zfftifc(3,ngrid,1,zv)
  zv(:)=zv(:)*cfunir(:)
  call zfftifc(3,ngrid,-1,zv)
  do ig=1,ngp
    unkit_out(ig,ispn)=zv(igfft(igpig(ig)))
  enddo !ig
enddo !ispn
!write(*,*)sum(abs(unkmt_out-unkmt_old))
!write(*,*)sum(abs(unkit_out-unkit_old))
!unkmt_out=0.5d0*unkmt_out+0.5d0*unkmt_old
!unkit_out=0.5d0*unkit_out+0.5d0*unkit_old
deallocate(unkmt_old,unkit_old)
deallocate(zv)
return
end subroutine
