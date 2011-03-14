subroutine sic_genwanlm(fout,n,wantp,wanlm)
use modmain
use mod_sic
implicit none
integer, intent(in) :: fout
integer, intent(in) :: n
complex(8), intent(inout) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(inout) :: wanlm(lmmaxwan,s_nr,nspinor)
! local variables
integer ispn,lm,lm1,flm(lmmaxwan,nspinor),j,ir,l,m
real(8) a,b,t1,t2
complex(8) zt1,zlm(lmmaxwan)
complex(8), external :: zdotc
complex(8), allocatable :: wantp1(:,:,:)
! convert to spherical harmonics
do ispn=1,nspinor
  call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,&
    wantp(1,1,ispn),s_ntp,zzero,wanlm,lmmaxwan)
  !call nfsft_bt(lmaxwan,s_ntp,s_nr,wantp,wanlm) 
enddo
! convert back to spherical coordinates
allocate(wantp1(s_ntp,s_nr,nspinor))
do ispn=1,nspinor
  call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,wanlm(1,1,ispn),&
    lmmaxwan,zzero,wantp1(1,1,ispn),s_ntp)
enddo
if (wproc) then
  write(fout,*)
  write(fout,'(" n : ",I4)')n
  write(fout,'(80("-"))')
  write(fout,'("imaginary part of original wantp (total, average)         : ",2G18.10)')&
    sum(abs(dimag(wantp))),sum(abs(dimag(wantp)))/s_nr/s_ntp
  write(fout,'("imaginary part of back-transformed wantp (total, average) : ",2G18.10)')&
    sum(abs(dimag(wantp1))),sum(abs(dimag(wantp1)))/s_nr/s_ntp
  write(fout,'("average difference of functions                           : ",G18.10)')&
    sum(abs(wantp-wantp1))/s_nr/s_ntp
  write(fout,'("total difference of functions                           : ",G18.10)')&
    sum(abs(wantp-wantp1))
endif
!if (wproc) then
!  do ir=1,s_nr
!    write(400,*)s_r(ir),sum(abs(wantp(:,ir,1)-wantp1(:,ir,1)))
!  enddo
!endif
!if (wproc) then
!  do ir=1,s_nr
!    write(400,*)s_r(ir),sum(abs(dimag(wantp(:,ir,1))))
!  enddo
!endif
deallocate(wantp1)

if (wproc) then
  t2=0.d0
  do l=0,lmaxwan
    do ispn=1,nspinor
      t1=0.d0
      do m=-l,l
        lm=idxlm(l,m)
        do ir=1,s_nr
          t1=t1+(abs(wanlm(lm,ir,ispn))**2)*s_rw(ir)
        enddo
      enddo
      write(fout,'("   l,ispn : ",2I4,"    norm : ",G18.10)')l,ispn,t1
      t2=t2+t1
    enddo
  enddo
  call flushifc(fout)
endif
call mpi_grid_bcast(t2)
wanlm=wanlm/sqrt(t2)
wantp=wantp/sqrt(t2)
return

! filter out small lm- contributions
flm=1
t2=0.d0
do lm=1,lmmaxwan
  do ispn=1,nspinor
    t1=0.d0
    do ir=1,s_nr
      t1=t1+(abs(wanlm(lm,ir,ispn))**2)*s_rw(ir)
    enddo
    if (t1.le.1d-8) then
      !wanlm(lm,:,ispn)=zzero
      flm(lm,ispn)=0
    !else
    !  t2=t2+t1
    endif
  enddo !ispn
enddo !lm
!write(*,*)"remaining weight : ",t2
! renormalize
!wanlm=wanlm/sqrt(t2)
!! fix tails
zt1=zzero
do ir=1,s_nr
  do ispn=1,nspinor
    zt1=zt1+zdotc(lmmaxwan,wanlm(1,ir,ispn),1,wanlm(1,ir,ispn),1)*s_rw(ir)
  enddo
  if (abs(zt1-1.d0).lt.0.01d0) exit
enddo
write(*,*)"tail radius : ",s_r(ir)
do ispn=1,nspinor
  do lm=1,lmmaxwan
    if (flm(lm,ispn).eq.1) then
      b=abs(wanlm(lm,ir,ispn))/abs(wanlm(lm,ir-1,ispn))
      b=log(b)/(s_r(ir-1)-s_r(ir))
      a=exp(b*s_r(ir))*abs(wanlm(lm,ir,ispn))
      write(*,*)"lm=",lm,"a,b=",a,b
      do j=ir,s_nr
        zt1=wanlm(lm,j,ispn)
        if (abs(zt1).gt.1d-10) then
          wanlm(lm,j,ispn)=a*exp(-b*s_r(j))*(zt1/abs(zt1))
        else
          wanlm(lm,j,ispn)=zzero
        endif
      enddo
    endif
  enddo
enddo
! renormalize
!zt1=zzero
!do ispn=1,nspinor
!  do ir=1,s_nr
!    zt1=zt1+zdotc(lmmaxwan,wanlm(1,ir,ispn),1,wanlm(1,ir,ispn),1)*s_rw(ir)
!  enddo
!enddo
!wanlm=wanlm/sqrt(abs(zt1))
!! convert to spherical coordinates
!call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,wanlm,&
!  lmmaxwan,zzero,wantp,s_ntp)
return
end
