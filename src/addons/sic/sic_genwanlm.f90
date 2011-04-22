subroutine sic_genwanlm(fout,n,wantp,wanlm)
use modmain
use mod_sic
implicit none
integer, intent(in) :: fout
integer, intent(in) :: n
complex(8), intent(inout) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(inout) :: wanlm(lmmaxwan,s_nr,nspinor)
! local variables
integer ispn,lm,lm1,flm(lmmaxwan,nspinor),j,ir,l,m,iter
integer nrloc,roffs,niter
real(8) t1,t2
real(8), allocatable :: tdiff(:)
complex(8), allocatable :: wantp1(:,:,:)
!
niter=5
nrloc=mpi_grid_map(s_nr,dim_k,offs=roffs)
wanlm=zzero
! convert to spherical harmonics
do ispn=1,nspinor
  call zgemm('T','N',lmmaxwan,nrloc,s_ntp,zone,s_ylmb,s_ntp,&
    wantp(1,roffs+1,ispn),s_ntp,zzero,wanlm(1,roffs+1,ispn),lmmaxwan)
enddo
allocate(wantp1(s_ntp,nrloc,nspinor))
allocate(tdiff(niter))
tdiff=0.d0
do iter=1,niter
! convert back to spherical coordinates
  do ispn=1,nspinor
    call zgemm('T','N',s_ntp,nrloc,lmmaxwan,zone,s_ylmf,lmmaxwan,wanlm(1,1+roffs,ispn),&
      lmmaxwan,zzero,wantp1(1,1,ispn),s_ntp)
  enddo
  wantp1(:,1:nrloc,:)=wantp(:,roffs+1:roffs+nrloc,:)-wantp1(:,1:nrloc,:)
  tdiff(iter)=sum(abs(wantp1))
! convert to spherical harmonics
  do ispn=1,nspinor
    call zgemm('T','N',lmmaxwan,nrloc,s_ntp,zone,s_ylmb,s_ntp,&
      wantp1(1,1,ispn),s_ntp,zone,wanlm(1,roffs+1,ispn),lmmaxwan)
  enddo
enddo
deallocate(wantp1)
call mpi_grid_reduce(tdiff(1),niter,dims=(/dim_k/))
call mpi_grid_reduce(wanlm(1,1,1),lmmaxwan*s_nr*nspinor,dims=(/dim_k/),all=.true.)
allocate(wantp1(s_ntp,s_nr,nspinor))
! convert back to spherical coordinates
do ispn=1,nspinor
  call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,wanlm(1,1,ispn),&
    lmmaxwan,zzero,wantp1(1,1,ispn),s_ntp)
enddo
if (wproc) then
  write(fout,*)
  write(fout,'(" n : ",I4)')n
  write(fout,'(80("-"))')
  do iter=1,niter
    write(fout,'(" iteration : ",I4,"  difference : ",G18.10)')iter,tdiff(iter)
  enddo
  write(fout,*)
  write(fout,'("imaginary part of original wantp (total, average)         : ",2G18.10)')&
    sum(abs(dimag(wantp))),sum(abs(dimag(wantp)))/s_nr/s_ntp
  write(fout,'("imaginary part of back-transformed wantp (total, average) : ",2G18.10)')&
    sum(abs(dimag(wantp1))),sum(abs(dimag(wantp1)))/s_nr/s_ntp
  write(fout,'("average difference of functions                           : ",G18.10)')&
    sum(abs(wantp-wantp1))/s_nr/s_ntp
  write(fout,'("total difference of functions                           : ",G18.10)')&
    sum(abs(wantp-wantp1))
  write(fout,*)
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
deallocate(wantp1)
deallocate(tdiff)
return
end
