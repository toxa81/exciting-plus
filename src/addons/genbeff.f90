subroutine genbeff
use modmain
implicit none
integer is,ia,ias,nr,ir,i,ig
integer l1,l2,l3,m2,m3,lm3
integer ilo,ilo1,ilo2,io,io1,io2
real(8) t1,cb
real(8), allocatable :: bmt(:,:,:,:)
complex(8), allocatable :: zfft(:) 
! automatic arrays
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
!
if (.not.spinpol.or.tsveqn) return
allocate(bmt(lmmaxvr,nrmtmax,natmtot,ndmag))
bmt=0.d0
cb=gfacte/(4.d0*solsc)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do i=1,ndmag
    bmt(:,:,ias,i)=bxcmt(:,:,ias,i)
  enddo
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bmt(1,:,ias,ndmag)=bmt(1,:,ias,ndmag)+t1/y00
enddo
baa=0.d0
bloa=0.d0
blolo=0.d0
! begin loops over atoms and species
do is=1,nspecies
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
    do l1=0,lmaxapw
      do io1=1,apword(l1,is)
        do l2=0,lmaxapw
          do io2=1,apword(l2,is)
            do lm3=1,lmmaxvr
              do i=1,ndmag
                do ir=1,nr
                  t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir) 
                  fr(ir)=t1*bmt(lm3,ir,ias,i)
                enddo
                call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                baa(lm3,io1,l1,io2,l2,ias,i)=gr(nr)
              end do
            end do
          end do
        end do
      end do
    end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
    do ilo=1,nlorb(is)
      l1=lorbl(ilo,is)
      do l2=0,lmaxapw
        do io2=1,apword(l2,is)
          do lm3=1,lmmaxvr
            do i=1,ndmag
              do ir=1,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir)
                fr(ir)=t1*bmt(lm3,ir,ias,i)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              bloa(lm3,ilo,io2,l2,ias,i)=gr(nr)
            end do
          end do
        end do
      end do
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do ilo1=1,nlorb(is)
      l1=lorbl(ilo1,is)
      do ilo2=1,nlorb(is)
        l2=lorbl(ilo2,is)
        do lm3=1,lmmaxvr
          do i=1,ndmag
            do ir=1,nr
              t1=lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
              fr(ir)=t1*bmt(lm3,ir,ias,i)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            blolo(lm3,ilo1,ilo2,ias,i)=gr(nr)
          end do
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do
deallocate(bmt)
allocate(zfft(ngrtot))
beffig=zzero
do i=1,ndmag
  do ir=1,ngrtot
    zfft(ir)=zone*(bxcir(ir,i)+cb*bfieldc(3))*cfunir(ir)
  enddo
  call zfftifc(3,ngrid,-1,zfft)
  do ig=1,ngvec
    beffig(ig,i)=zfft(igfft(ig))
  enddo
enddo
deallocate(zfft)
return
end subroutine
