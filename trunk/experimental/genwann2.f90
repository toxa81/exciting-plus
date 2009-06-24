subroutine genwann2(ik,apwalm,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
!complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: wannmt_new(:,:,:,:,:)
complex(8), allocatable :: wannit_new(:,:,:)
integer ispn,n,m,i,is,ia,ias,ir,irc,io,l,lm,ig,ifg
integer, external :: ikglob
complex(8) z1
complex(8), external :: zfmtinp,zfinp
complex(8), external :: zfint

!allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)))

allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))


  


allocate(wannmt_new(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nmax))
allocate(wannit_new(ngrtot,nspinor,wann_nmax))

wannmt_new=zzero
wannit_new=zzero

do ispn=1,nspinor
  do n=1,nwann(ispn)
  
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        irc=0
        do ir=1,nrmt(is),lradstp
          irc=irc+1
          do io=1,nrfmax
            do l=0,lmaxvr
              do m=-l,l
                lm=idxlm(l,m)
                wannmt_new(lm,irc,ias,ispn,n)=wannmt_new(lm,irc,ias,ispn,n)+&
                  wann_unkmt(lm,io,ias,n,ispn,ik)*urf(ir,l,io,ias)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
   
    do ig=1,ngk(1,ikglob(ik))
      ifg=igfft(igkig(ig,1,ik))
      wannit_new(ifg,ispn,n)=wann_unkit(ig,n,ispn,ik)/sqrt(omega)
    end do
    call zfftifc(3,ngrid,1,wannit_new(:,ispn,n))
  enddo
enddo

wannmt(:,:,:,:,:,ik)=0.85*wannmt(:,:,:,:,:,ik)+0.15*wannmt_new(:,:,:,:,:)
wannit(:,:,:,ik)=0.85*wannit(:,:,:,ik)+0.15*wannit_new(:,:,:)

!do ispn=1,nspinor
!  do n=1,nwann(ispn)
!    do m=1,nwann(ispn)
!      z1=zzero
!      z1=zfinp(.true.,wannmt(:,:,:,ispn,n,ik),wannmt(:,:,:,ispn,m,ik),wannit(:,ispn,n,ik), &
!        wannit(:,ispn,m,ik))
!      write(*,*)'ispn=',ispn,'n,m=',n,m,'norm=',z1
!    enddo
!  enddo
!enddo
!
!wannmt(:,:,:,:,:,ik)=zzero
!wannit(:,:,:,ik)=zzero
!write(*,*)'WF2'
!call genwfsv(.false.,ngk(1,ikglob(ik)),igkig(:,1,ik),evalsv(1,ikglob(ik)), &
!  apwalm,evecfv,evecsv,wfmt,wfir)
!do ispn=1,nspinor
!  do n=1,nwann(ispn)
!    do i=1,nstfv
!      wannmt(:,:,:,ispn,n,ik)=wannmt(:,:,:,ispn,n,ik)+&
!        wann_c(n,i,ispn,ik)*wfmt(:,:,:,ispn,i+(ispn-1)*nstfv)
!      wannit(:,ispn,n,ik)=wannit(:,ispn,n,ik)+&
!        wann_c(n,i,ispn,ik)*wfir(:,ispn,i+(ispn-1)*nstfv)
!    enddo
!  enddo
!enddo
!do ispn=1,nspinor
!  do n=1,nwann(ispn)
!    do m=1,nwann(ispn)
!      z1=zzero
!      z1=zfinp(.false.,wannmt(:,:,:,ispn,n,ik),wannmt(:,:,:,ispn,m,ik),wannit(:,ispn,n,ik), &
!        wannit(:,ispn,m,ik))
!      write(*,*)'ispn=',ispn,'n,m=',n,m,'norm=',z1
!    enddo
!  enddo
!enddo
  
  



deallocate(wfmt,wfir,zrhomt,zrhoir,wannmt_new,wannit_new)
return
end