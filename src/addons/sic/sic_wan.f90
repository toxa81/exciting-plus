subroutine sic_wan(fout)
use modmain
use mod_sic
use mod_nrkp
use mod_wannier
use mod_linresp
implicit none
! arguments
integer, intent(in) :: fout
! local variables
complex(8), allocatable :: ovlp(:)
integer n,ispn,vl(3),n1,i,j,j1
real(8) t1,t2,t3
complex(8) z1

if (wproc) then
  write(fout,*)
  write(fout,'("generating Wannier functions on a grid")')
  write(fout,'(80("-"))')
endif
call timer_reset(1)
call timer_reset(2)
call timer_start(3,reset=.true.)
call sic_genwan
call timer_stop(3)
deallocate(wann_unkmt)
deallocate(wann_unkit)
if (wproc) then
  write(fout,*)
  write(fout,'("local time for muffin-tin part (sec.)   : ",F8.3)')timer_get_value(1)
  write(fout,'("local time for interstitial part (sec.) : ",F8.3)')timer_get_value(2)
  write(fout,'("total global time : ",F8.3)')timer_get_value(3)
  call flushifc(fout)
endif
allocate(ovlp(sic_wantran%nwt))
ovlp=zzero
! compute overlap integrals 
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  do ispn=1,nspinor
    ovlp(i)=ovlp(i)+sic_dot_ll(sic_orbitals%wanmt(1,1,1,ispn,j),&
      sic_orbitals%wanir(1,1,ispn,j),sic_orbitals%wanmt(1,1,1,ispn,j1),&
      sic_orbitals%wanir(1,1,ispn,j1),vl,sic_orbitals%twanmt(1,1,n),&
      sic_orbitals%twanmt(1,1,n1))
  enddo
enddo
! check orthonormality
t1=0.d0
t2=0.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  z1=ovlp(i)
  if (n.eq.n1.and.all(vl.eq.0)) then
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
  t1=t1+abs(z1)
enddo
if (wproc) then
  write(fout,*)
!  write(fout,'("Wannier overlap integrals (n n1 <w_n|w_n1>)")')
!  do i=1,sic_wantran%nwt
!    vl(:)=sic_wantran%iwt(3:5,i)
!    if (all(vl.eq.0)) then
!      write(151,'(I4,4X,I4,4X,2G18.10)')sic_wantran%iwt(1:2,i),&
!        dreal(ovlp(i)),dimag(ovlp(i))
!    endif
!  enddo
  write(fout,'("Wannier overlap integrals: (<W_n|W_n>)")')
  do i=1,sic_wantran%nwan
    n=sic_wantran%iwan(i)
    j=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'("  n : ",I4,8X,2G18.10)')n,dreal(ovlp(j)),dimag(ovlp(j))
  enddo
  write(fout,*)
  write(fout,'("Maximum deviation from norm                 : ",F12.6)')t2
  call timestamp(fout,"done with Wannier functions")
  call flushifc(151)
endif
deallocate(ovlp)
return
end