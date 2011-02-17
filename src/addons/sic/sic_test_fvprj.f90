subroutine sic_test_fvprj(fout)
use modmain
use mod_sic
use mod_nrkp
implicit none
integer, intent(in) :: fout
integer ikloc,ik,ist,ir,itp,i,j,n,ispn,istsv
real(8) vrc(3),t1
complex(8) zt1,zt2
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: evecfvnr(:,:,:)
complex(8), allocatable :: evecsvnr(:,:)
complex(8), allocatable :: wffvtp(:,:)
complex(8), allocatable :: wffvlm(:,:)
complex(8), external :: zdotc
!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wffvmt(lmmaxvr,nufrmax,natmtot,nstfv))
allocate(evecfvnr(nmatmax,nstfv,nspnfv))
allocate(evecsvnr(nstsv,nstsv))
allocate(wffvtp(s_ntp,s_nr))
allocate(wffvlm(lmmaxwan,s_nr))
open(220,file="sic_test_fvprj.out",form="formatted",status="replace")
t1=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(220,'("ik : ",I4,"  vklnr : ",3G18.10)')ik,vklnr(:,ik)
  call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnr)
  call getevecsv(vklnr(1,ik),evecsvnr)
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),&
    sfacgknr(1,1,ikloc),apwalm)
  call genwffvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvnr,apwalm,wffvmt)
  do ist=1,nstfv
    write(220,'("  ist : ",I4)')ist
    do ir=1,s_nr
      do itp=1,s_ntp
        vrc(:)=s_spx(:,itp)*s_r(ir)
        call s_get_wffvval(vrc,ngknr(ikloc),vkcnr(1,ik),vgkcnr(1,1,ikloc), &
          wffvmt(1,1,1,ist),evecfvnr(1,ist,1),wffvtp(itp,ir))
      enddo
    enddo
    call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,wffvtp,&
      s_ntp,zzero,wffvlm,lmmaxwan)
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        zt1=zzero
        do ir=1,s_nr
          zt1=zt1+zdotc(lmmaxwan,s_wanlm(1,ir,ispn,j),1,wffvlm(1,ir),1)*s_rw(ir)
        enddo
        zt2=zzero
        istsv=ist+(ispn-1)*nstfv
        do i=1,nstsv
          zt2=zt2+dconjg(wanncnrloc(n,i,ikloc)*evecsvnr(istsv,i))
        enddo !i
        t1=max(t1,abs(zt1-zt2))
        write(220,'("    n : ",I4,"  ispn : ",I4,&
          &"  numerical : ",2G18.10,"  analytical : ",2G18.10,&
          &"  diff : ",G18.10)')n,ispn,dreal(zt1),dimag(zt1),&
          dreal(zt2),dimag(zt2),abs(zt1-zt2)
      enddo !ispn
    enddo !j
  enddo !ist
  write(220,*)
  call flushifc(220)
enddo !ikloc
call mpi_grid_reduce(t1,dims=(/dim_k/),op=op_max)
if (wproc) then
  write(fout,*)
  write(fout,'("Maximum deviation for <W_n|\phi_{jk}> : ",G18.10)')t1
endif
close(220)
deallocate(apwalm)
deallocate(wffvmt)
deallocate(evecfvnr)
deallocate(evecsvnr)
deallocate(wffvtp)
deallocate(wffvlm)
return
end

