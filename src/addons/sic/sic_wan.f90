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
integer n,ispn,vl(3),n1,i,j
real(8) t1,t2,t3
complex(8) z1

if (wproc) then
  write(fout,*)
  write(fout,'("sic_wan.f90")')
  write(fout,'("generate Wannier functions on a grid")')
  write(fout,'(80("-"))')
endif
call timer_reset(1)
call timer_reset(2)
call sic_genwan
deallocate(wann_unkmt)
deallocate(wann_unkit)
if (wproc) then
  write(fout,*)
  write(fout,'("time for muffin-tin part (sec.)   : ",F8.3)')timer_get_value(1)
  write(fout,'("time for interstitial part (sec.) : ",F8.3)')timer_get_value(2)
  call flushifc(fout)
endif
allocate(ovlp(sic_wantran%nwt))
ovlp=zzero
! compute overlap integrals
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  do ispn=1,nspinor
    ovlp(i)=ovlp(i)+sic_dot_ll(wanmt(1,1,1,ispn,n),wanir(1,1,ispn,n),&
      wanmt(1,1,1,ispn,n1),wanir(1,1,ispn,n1),vl,twanmt(1,1,n),&
      twanmt(1,1,n1))
  enddo
enddo
! check orthonormality
t1=0.d0
t2=0.d0
t3=0.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  z1=ovlp(i)
  if (n.eq.n1.and.vl(1).eq.0.and.vl(2).eq.0.and.vl(3).eq.0) then
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
  t1=t1+abs(z1)
  t3=max(t3,abs(ovlp(i)-dconjg(ovlp(j))))
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
  write(fout,'("Wannier overlap integrals (<w_n|w_n>)")')
  do n=1,nwantot
    j=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'(I4,4X,2G18.10)')n,dreal(ovlp(j)),dimag(ovlp(j))
  enddo
  write(fout,*)
  write(fout,'("Maximum deviation from norm                 : ",F12.6)')t2
  write(fout,'("Maximum of <w_n|w_{n1,T}> - <w_n1|w_{n,-T}> : ",G18.10)')t3
  call timestamp(fout,"done with Wannier functions")
  call flushifc(151)
endif
deallocate(ovlp)
return
end