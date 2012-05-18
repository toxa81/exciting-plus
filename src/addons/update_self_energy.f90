subroutine update_self_energy(iq)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer, intent(in) :: iq
!
integer ikloc,iw,ig,ig1,ig2,n,jk,i,j,ik
real(8) ppa_e0
complex(8) zw
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: ppa_w(:,:)
complex(8), allocatable :: ppa_r(:,:)
!
allocate(chi0(ngvecme,ngvecme))
allocate(epsinv(ngvecme,ngvecme,2))
allocate(megqblh2(nstsv*nstsv,ngvecme))
allocate(ppa_w(ngvecme,ngvecme))
allocate(ppa_r(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
!
ppa_e0=1.d0
!
krnl=zzero
do ig=1,ngvecme
  krnl(ig,ig)=vhgq(ig,iq)
enddo
epsinv=zzero
do iw=1,2
  if (iw.eq.1) then
    zw=dcmplx(0.d0,lr_eta)
  else
    zw=dcmplx(0.d0,ppa_e0)
  endif
  chi0=zzero
  do ikloc=1,nkptnrloc
    if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call genchi0blh(ikloc,zw,chi0)
    endif
  enddo
  call mpi_grid_reduce(chi0(1,1),ngvecme*ngvecme,dims=(/dim_k/),all=.true.)
  chi0(:,:)=chi0(:,:)/nkptnr/omega
! compute epsilon=1-chi0*v
  do ig=1,ngvecme
    epsinv(ig,ig,iw)=zone
  enddo
  call zgemm('N','N',ngvecme,ngvecme,ngvecme,-zone,chi0,ngvecme,&
    &krnl,ngvecme,zone,epsinv(1,1,iw),ngvecme)
! inverse epsilon: eps^{-1}=(1-chi0*v)^{-1} = 1+chi*v
  call invzge(epsinv(1,1,iw),ngvecme)
! find frequency dependent part of inverse epsilon
  do ig=1,ngvecme
    epsinv(ig,ig,iw)=epsinv(ig,ig,iw)-zone
  enddo
enddo
! compute coefficients of the plasmone-pole approximation
do ig1=1,ngvecme
  do ig2=1,ngvecme
    ppa_w(ig1,ig2)=ppa_e0*sqrt(epsinv(ig1,ig2,2)/(epsinv(ig1,ig2,1)-epsinv(ig1,ig2,2)))
    ppa_r(ig1,ig2)=-epsinv(ig1,ig2,1)*ppa_w(ig1,ig2)/2.d0
  enddo
enddo
if (all(vqm(:,iq).eq.0)) then
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      if (igqig(ig1,iq).eq.1.and.ig1.eq.ig2) then
        ppa_r(ig1,ig1)=ppa_r(ig1,ig1)*aq0(iq)
      else
        ppa_r(ig1,ig2)=ppa_r(ig1,ig2)*0.125d0
      endif
    enddo !ig2      
  enddo !ig1
endif
!
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(3,ik)
  do i=1,namegqblh(ikloc)
    n=bamegqblh(2,i,ikloc)
    j=bamegqblh(1,i,ikloc)
    do ig1=1,ngvecme
      do ig2=1,ngvecme
        do iw=1,lr_nw
          self_energy_c(iw,n,ikloc)=self_energy_c(iw,n,ikloc)-dconjg(amegqblh(i,ig1,ikloc))*vhgq(ig1,iq)*&
            &ppa_r(ig1,ig2)*(occsvnr(j,jk)/(evalsvnr(j,jk)-dconjg(lr_w(iw))-ppa_w(ig1,ig2)) - &
            &(occmax-occsvnr(j,jk))*(evalsvnr(j,jk)-lr_w(iw)+ppa_w(ig1,ig2)))*amegqblh(i,ig2,ikloc)
        enddo
      enddo
    enddo
  enddo
enddo

deallocate(chi0)
deallocate(megqblh2)
deallocate(ppa_w)
deallocate(ppa_r)
deallocate(krnl)
deallocate(epsinv)
return
end subroutine
