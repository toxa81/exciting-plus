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
    wantp(1,1,ispn),s_ntp,zzero,wanlm(1,1,ispn),lmmaxwan)
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
  do ispn=1,nspinor
    do l=0,lmaxwan
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
!call mpi_grid_bcast(t2)
!wanlm=wanlm/sqrt(t2)
!wantp=wantp/sqrt(t2)
return
end
