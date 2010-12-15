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
real(8) t1,t2,vrc(3)
complex(8) z1,z2(4)
real(8), allocatable :: xmt(:,:,:,:)
real(8), allocatable :: xir(:,:,:)
real(8), allocatable :: spread(:)

if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("generating Wannier functions on a grid")')
  write(fout,'(80("="))')
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
  write(fout,'("overlap integrals",13X,"<W_n|W_n>")')
  write(fout,'(60("-"))')
  do i=1,sic_wantran%nwan
    n=sic_wantran%iwan(i)
    j=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'("  n : ",I4,8X,2G18.10)')n,dreal(ovlp(j)),dimag(ovlp(j))
  enddo
  write(fout,'(60("-"))')
  write(fout,'("maximum deviation from norm : ",F12.6)')t2
endif
deallocate(ovlp)
allocate(xmt(lmmaxvr,nmtloc,sic_orbitals%ntr,4))
allocate(xir(ngrloc,sic_orbitals%ntr,4))
call sic_gen_r(xmt,xir)
allocate(spread(sic_wantran%nwan))
spread=0.d0
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  z2=zzero
  do i=1,4
    do ispn=1,nspinor
      z2(i)=z2(i)+sic_int_zdz(sic_orbitals%wanmt(1,1,1,ispn,j),&
        sic_orbitals%wanir(1,1,ispn,j),xmt(1,1,1,i),xir(1,1,i),&
        sic_orbitals%wanmt(1,1,1,ispn,j),sic_orbitals%wanir(1,1,ispn,j),&
        sic_orbitals%twanmtuc(1,n))
    enddo
  enddo
  do i=1,3
    vrc(i)=dreal(z2(i))
  enddo
  spread(j)=dreal(z2(4))-dot_product(vrc(:),vrc(:))
enddo
if (wproc) then
  write(fout,*)
  write(fout,'("quadratic spreads",3X,"<r^2> - <r>^2  [a.u.]^2")')
  write(fout,'(60("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(151,'("  n : ",I4,8X,G18.10)')n,spread(j)
  enddo
  write(fout,'(60("-"))')
  write(fout,'("total spread : ",F12.6)')sum(spread)
  call timestamp(fout,"done with Wannier functions")
  call flushifc(151)  
endif
deallocate(xmt,xir,spread)
return
end