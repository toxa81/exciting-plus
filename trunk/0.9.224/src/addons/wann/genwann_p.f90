subroutine genwann_p(ikloc,evecfv,evecsv)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)

complex(8), allocatable :: zt2(:,:,:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
integer m1,m2,j1,j2
integer, external :: ikglob

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(pmat(3,nstsv,nstsv))
call match(ngk(1,ikglob(ikloc)),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
  sfacgk(1,1,1,ikloc),apwalm)
call genpmat(ngk(1,ikglob(ikloc)),igkig(1,1,ikloc),gkc(1,1,ikloc),&
  apwalm,evecfv,evecsv,pmat)
!do m1=1,nstsv
!  do m2=m1,nstsv
!    write(*,*)'ik=',ikloc,'m1=',m1,'m2=',m2,'grad=',pmat(:,m1,m2)
!  enddo
!enddo
!call pstop
! compute p_nn'(k)=<W_n|\grad|W_n'> 
allocate(zt2(3,nwann,nwann))
zt2=zzero
do m1=1,nwann
  do m2=1,nwann
    do j1=1,nstsv
      do j2=1,nstsv
        zt2(:,m1,m2)=zt2(:,m1,m2)+dconjg(wann_c(m1,j1,ikloc))*wann_c(m2,j2,ikloc)*&
          pmat(:,j1,j2)
      enddo
    enddo
  enddo
enddo
wann_p(:,:,:,ikglob(ikloc))=zt2(:,:,:)
deallocate(zt2,apwalm,pmat)
return
end
