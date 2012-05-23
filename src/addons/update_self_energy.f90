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
integer ist1,ist2
real(8) ppa_e0
complex(8) zw,zt1
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: ppa_w(:,:)
complex(8), allocatable :: ppa_r(:,:)
complex(8), allocatable :: zm(:,:)
complex(8), allocatable :: ame(:,:)
!
allocate(chi0(ngq(iq),ngq(iq)))
allocate(epsinv(ngq(iq),ngq(iq),2))
allocate(megqblh2(nstsv*nstsv,ngq(iq)))
allocate(ppa_w(ngq(iq),ngq(iq)))
allocate(ppa_r(ngq(iq),ngq(iq)))
allocate(krnl(ngq(iq),ngq(iq)))
allocate(zm(ngq(iq),ngq(iq)))
allocate(ame(ngq(iq),nstsv*nstsv))
!
ppa_e0=1.d0
!
krnl=zzero
do ig=1,ngq(iq)
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
      call genchi0blh(ikloc,ngq(iq),zw,chi0)
    endif
  enddo
  call mpi_grid_reduce(chi0(1,1),ngq(iq)*ngq(iq),dims=(/dim_k/),all=.true.)
  chi0(:,:)=chi0(:,:)/nkptnr/omega
! compute epsilon=1-chi0*v
  do ig=1,ngq(iq)
    epsinv(ig,ig,iw)=zone
  enddo
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0,ngq(iq),&
    &krnl,ngq(iq),zone,epsinv(1,1,iw),ngq(iq))
! inverse epsilon: eps^{-1}=(1-chi0*v)^{-1} = 1+chi*v
  call invzge(epsinv(1,1,iw),ngq(iq))
! find frequency dependent part of inverse epsilon
  do ig=1,ngq(iq)
    epsinv(ig,ig,iw)=epsinv(ig,ig,iw)-zone
  enddo
enddo
! compute coefficients of the plasmone-pole approximation
do ig1=1,ngq(iq)
  do ig2=1,ngq(iq)
    ppa_w(ig1,ig2)=ppa_e0*sqrt(epsinv(ig1,ig2,2)/(epsinv(ig1,ig2,1)-epsinv(ig1,ig2,2)))
    ppa_r(ig1,ig2)=-epsinv(ig1,ig2,1)*ppa_w(ig1,ig2)/2.d0
  enddo
enddo

! restore q->0 matrix elements
if (vq_gamma(iq)) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do ig=1,ngq(iq)
      if (igqig(ig,iq).eq.1) then
        amegqblh(:,ig,ikloc)=zzero
        do i=1,nmegqblh(ikloc)
          ist1=bmegqblh(1,i,ikloc)
          ist2=bmegqblh(2,i,ikloc)
          if (ist1.eq.ist2) then
            amegqblh(i,ig,ikloc)=zone
          endif
        enddo
      endif
    enddo
  enddo
endif

do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(3,ik)
  ! change order of indices
  do ig=1,ngq(iq)
    ame(ig,:)=amegqblh(:,ig,ikloc)
  enddo
  do iw=1,lr_nw
    j=-1
    do i=1,namegqblh(ikloc)
      n=bamegqblh(2,i,ikloc)
      if (bamegqblh(1,i,ikloc).ne.j) then
        j=bamegqblh(1,i,ikloc)
        do ig1=1,ngq(iq)
          do ig2=1,ngq(iq)
            zm(ig1,ig2)=-wtvhgq(ig1,iq)*ppa_r(ig1,ig2)*(occsvnr(j,jk)/(evalsvnr(j,jk)-dconjg(lr_w(iw))-ppa_w(ig1,ig2))-&
            &(occmax-occsvnr(j,jk))/(evalsvnr(j,jk)-lr_w(iw)+ppa_w(ig1,ig2)))
          enddo
        enddo
      endif
      do ig2=1,ngq(iq)
        zt1=zzero
        do ig1=1,ngq(iq)
          zt1=zt1+dconjg(ame(ig1,i))*zm(ig1,ig2)
        enddo
        self_energy_c(iw,n,ikloc)=self_energy_c(iw,n,ikloc)+&
          &zt1*ame(ig2,i)
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
deallocate(zm)
deallocate(ame)
return
end subroutine
