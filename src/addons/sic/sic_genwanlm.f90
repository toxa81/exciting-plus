subroutine sic_genwanlm(fout,n,wantp,wanlm)
use modmain
use mod_sic
implicit none
integer, intent(in) :: fout
integer, intent(in) :: n
complex(8), intent(in) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(inout) :: wanlm(lmmaxwan,s_nr,nspinor)
! local variables
integer ispn,lm,ir,l,m
real(8) t1
complex(8), allocatable :: wantp1(:,:,:)
real(8), allocatable :: f(:)
!
do ispn=1,nspinor
  call sic_zbsht(s_nr,wantp(1,1,ispn),wanlm(1,1,ispn))
enddo
allocate(wantp1(s_ntp,s_nr,nspinor))
allocate(f(s_nr))
! convert to spherical coordinates
do ispn=1,nspinor
  call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,wanlm(1,1,ispn),&
    &lmmaxwan,zzero,wantp1(1,1,ispn),s_ntp)
enddo
if (wproc) then
  write(fout,*)
  write(fout,'(" n : ",I4)')n
  write(fout,'(80("-"))')
  write(fout,*)
  write(fout,'("imaginary part of the original wf (total, average)      : ",2G18.10)')&
    &sum(abs(dimag(wantp))),sum(abs(dimag(wantp)))/s_nr/s_ntp/nspinor
  write(fout,'("imaginary part of the reconstructed wf (total, average) : ",2G18.10)')&
    &sum(abs(dimag(wantp1))),sum(abs(dimag(wantp1)))/s_nr/s_ntp/nspinor
  write(fout,'("difference of functions (total,average)                 : ",2G18.10)')&
    &sum(abs(wantp-wantp1)),sum(abs(wantp-wantp1))/s_nr/s_ntp/nspinor
  write(fout,*)
  do ispn=1,nspinor
    do l=0,lmaxwan
      t1=0.d0
      do m=-l,l
        lm=idxlm(l,m)
        do ir=1,s_nr
          f(ir)=(abs(wanlm(lm,ir,ispn))**2)
        enddo
        t1=t1+rintegrate(s_nr,s_r,f)
      enddo !m
      write(fout,'("   l,ispn : ",2I4,"    norm : ",G18.10)')l,ispn,t1
    enddo !l
  enddo !ispn
  call flushifc(fout)
endif
deallocate(wantp1,f)
return
end

