subroutine sic_hunif_sv(ikloc,evecfd)
use modmain
use mod_sic
use mod_seceqn
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: evecfd(nspinor*nmatmax,nstsv)
!
integer ik,ig,ist,ispn,lm,ir,is,ic,ias,j,io,l,j1,j2
integer vtrc(3),n1,n2,i1,i2,i
real(8) vtrl(3),en
complex(8) zt1,zt2,expikt
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: h(:,:),vk(:,:)
complex(8), allocatable :: evecfd_new(:,:)
!
if (.not.tsic_wv) return
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)   
allocate(wfsvmt(lmmaxapw,nufrmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
wfsvmt=zzero
wfsvit=zzero
call genapwalm(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
  sfacgk(:,:,1,ikloc),apwalm)
call genwfsvc(lmaxapw,lmmaxapw,ngk(1,ik),nstsv,apwalm,evecfd,wfsvmt,wfsvit) 
deallocate(apwalm)
!sic_wb(:,:,:,ikloc)=zzero
!sic_wvb(:,:,:,ikloc)=zzero
!do j=1,sic_wantran%nwan
!  do ispn=1,nspinor
!    do ias=1,natmtot
!      is=ias2is(ias)
!      ic=ias2ic(ias)
!      do lm=1,lmmaxapw
!        l=lm2l(lm)
!        do io=1,nufr(l,is)
!          zt1=zzero
!          zt2=zzero
!          do ir=1,nrmt(is)
!            zt1=zt1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*ufr(ir,l,io,ic)*mt_rw(ir,is)
!            zt2=zt2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*ufr(ir,l,io,ic)*mt_rw(ir,is)
!          enddo
!          do ist=1,nstsv
!            sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+zt1*wfsvmt(lm,io,ias,ispn,ist)
!            sic_wvb(j,ist,1,ikloc)=sic_wvb(j,ist,1,ikloc)+zt2*wfsvmt(lm,io,ias,ispn,ist)
!          enddo
!        enddo
!      enddo
!    enddo
!    do ist=1,nstsv
!      do ig=1,ngk(1,ik)
!        sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+&
!          dconjg(s_wkit(ig,ispn,j,ikloc))*wfsvit(ig,ispn,ist)
!        sic_wvb(j,ist,1,ikloc)=sic_wvb(j,ist,1,ikloc)+&
!          dconjg(s_wvkit(ig,ispn,j,ikloc))*wfsvit(ig,ispn,ist)
!      enddo
!    enddo
!  enddo
!enddo
call sic_genbprj(ikloc,wfsvmt=wfsvmt,wfsvit=wfsvit)

write(*,*)"in sic_hunif_sv"
do j1=1,sic_wantran%nwan
  do j2=1,sic_wantran%nwan
    zt1=zzero
    do j=1,nstsv
      zt1=zt1+dconjg(sic_wb(j1,j,1,ikloc))*sic_wb(j2,j,1,ikloc)
    enddo
    write(*,*)j1,j2,zt1
  enddo
enddo

deallocate(wfsvmt,wfsvit)

sic_wan_h0k(:,:,ikloc)=zzero
do j1=1,sic_wantran%nwan 
  do j2=1,sic_wantran%nwan 
    do ist=1,nstsv
      sic_wan_h0k(j1,j2,ikloc)=sic_wan_h0k(j1,j2,ikloc)+&
        sic_wb(j1,ist,1,ikloc)*dconjg(sic_wb(j2,ist,1,ikloc))*evalsv(ist,ik)
    enddo
  enddo
enddo

allocate(h(nstsv,nstsv))
h=zzero
do j=1,nstsv
  h(j,j)=evalsv(j,ik)
enddo

! compute V_{nn'}(k)
allocate(vk(sic_wantran%nwan,sic_wantran%nwan))
vk=zzero
do i=1,sic_wantran%nwt
  n1=sic_wantran%iwt(1,i)
  j1=sic_wantran%idxiwan(n1)
  n2=sic_wantran%iwt(2,i)
  j2=sic_wantran%idxiwan(n2)
  vtrl(:)=sic_wantran%iwt(3:5,i)
  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  vk(j1,j2)=vk(j1,j2)+expikt*sic_vme(i)
enddo







do i1=1,nstsv
  do i2=1,nstsv
    do j1=1,sic_wantran%nwan
      n1=sic_wantran%iwan(j1)
      en=dreal(sic_vme(sic_wantran%iwtidx(n1,n1,0,0,0)))!+sic_wan_e0(n1)
      do j2=1,sic_wantran%nwan
! 2-nd term : -\sum_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'}     
          !hunif(ist1,ist2)=hunif(ist1,ist2)-sic_wan_h0k(j1,j2,ikloc)*&
          !  dconjg(sic_wb(j1,istfv1,ispn1,ikloc))*sic_wb(j2,istfv2,ispn2,ikloc)
! 5-th term, first part: 
!  -\sum_{\alpha,\alpha'} ( P_{\alpha}V_{\alpha}P_{\alpha'} +
!                           P_{\alpha'}V_{\alpha}P_{\alpha} )
        h(i1,i2)=h(i1,i2)-&
          vk(j1,j2)*dconjg(sic_wb(j1,i1,1,ikloc))*&
          sic_wb(j2,i2,1,ikloc)-dconjg(vk(j1,j2))*&
          dconjg(sic_wb(j2,i1,1,ikloc))*sic_wb(j1,i2,1,ikloc)
      enddo !j2
! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
      h(i1,i2)=h(i1,i2)+&
        dconjg(sic_wb(j1,i1,1,ikloc))*sic_wb(j1,i2,1,ikloc)*en
! 5-th term, second part: \sum_{\alpha} P_{\alpha}V_{\alpha}+V_{\alpha}P_{\alpha}
      h(i1,i2)=h(i1,i2)+&
        dconjg(sic_wb(j1,i1,1,ikloc))*sic_wvb(j1,i2,1,ikloc)+&
        dconjg(sic_wvb(j1,i1,1,ikloc))*sic_wb(j1,i2,1,ikloc)
    enddo !j1
  enddo !i2 
enddo !i1

call diagzhe(nstsv,h,evalsv(1,ik))
allocate(evecfd_new(nspinor*nmatmax,nstsv))
evecfd_new=zzero
do j=1,nstsv
  do i=1,nstsv
    evecfd_new(:,j)=evecfd_new(:,j)+h(i,j)*evecfd(:,i)
  enddo
enddo
evecfd=evecfd_new
deallocate(h,vk,evecfd_new)










return
end subroutine
