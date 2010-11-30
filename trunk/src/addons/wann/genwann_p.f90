subroutine genwann_p(ikloc,evecfv,evecsv)
use modmain
use mod_wannier
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)

complex(8), allocatable :: zt2(:,:,:)
complex(8), allocatable :: pmat(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
integer m1,m2,j1,j2

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(pmat(3,nstsv,nstsv))
call match(ngk(1,ikloc),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
  sfacgk(1,1,1,ikloc),apwalm)
call genpmat(ngk(1,ikloc),igkig(1,1,ikloc),vgkc(1,1,1,ikloc),&
  apwalm,evecfv,evecsv,pmat)
! compute p_nn'(k)=<W_n|\grad|W_n'> 
allocate(zt2(3,nwantot,nwantot))
zt2=zzero
do m1=1,nwantot
  do m2=1,nwantot
    do j1=1,nstsv
      do j2=1,nstsv
        zt2(:,m1,m2)=zt2(:,m1,m2)+dconjg(wann_c(m1,j1,ikloc))*wann_c(m2,j2,ikloc)*&
          pmat(:,j1,j2)
      enddo
    enddo
  enddo
enddo
wann_p(:,:,:,ikloc)=zt2(:,:,:)
deallocate(zt2,apwalm,pmat)
return
end
