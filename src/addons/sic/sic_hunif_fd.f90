subroutine sic_hunif_fd(ikloc,nmatp,hunif,om)
use modmain
use mod_sic
use mod_wannier
! arguments
implicit none
integer, intent(in) :: ikloc
integer, intent(in) :: nmatp
complex(8), intent(inout) :: hunif(nmatp,nmatp,nspinor)
complex(8), intent(inout) :: om(nmatp,nmatp)
! local variables
integer i,ik,vtrl(3),n1,n2,i1,i2,j1,j2,ispn
real(8) vtrc(3),vn
complex(8) expikt,zt1
complex(8), allocatable :: vk(:,:)
character*500 msg,fname
logical, parameter :: tcheckherm=.false.
complex(8), allocatable :: hm1(:,:,:),om1(:,:),zm1(:,:),zm0(:,:)
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_sic_hunif)
sic_wan_h0k(:,:,ikloc)=zzero
! restore full hermitian matrix
do ispn=1,nspinor
  do j1=2,nmatp
    do j2=1,j1-1
      hunif(j1,j2,ispn)=dconjg(hunif(j2,j1,ispn))
    enddo
  enddo
  do j1=1,nmatp
    hunif(j1,j1,ispn)=zone*dreal(hunif(j1,j1,ispn))
  enddo
enddo
do j1=2,nmatp
  do j2=1,j1-1
    om(j1,j2)=dconjg(om(j2,j1))
  enddo
enddo
do j1=1,nmatp
  om(j1,j1)=zone*dreal(om(j1,j1))
enddo
allocate(hm1(nmatp,nmatp,nspinor),om1(nmatp,nmatp))
om1=om
call invzge(om1,nmatp)

!write(*,*)"ik=",ik
!do j1=1,sic_wantran%nwan
!  do j2=1,sic_wantran%nwan
!    zt1=zzero
!    do i1=1,nmatp
!      do i2=1,nmatp
!        zt1=zt1+sic_wb(j1,i1,1,ikloc)*om1(i1,i2)*dconjg(sic_wb(j2,i2,1,ikloc))
!      enddo
!    enddo
!    write(*,*)"j1,j2=",j1,j2,"  prod=",zt1
!  enddo
!enddo

allocate(zm1(nmatp,nmatp))
do ispn=1,nspinor
  call zgemm('N','N',nmatp,nmatp,nmatp,zone,om1,nmatp,hunif(1,1,ispn),nmatp,zzero,zm1,nmatp)
  call zgemm('N','N',nmatp,nmatp,nmatp,zone,zm1,nmatp,om1,nmatp,zzero,hm1(1,1,ispn),nmatp)
enddo
deallocate(zm1)
deallocate(om1)
allocate(zm1(sic_wantran%nwan,nmatp))
do ispn=1,nspinor
  call zgemm('N','N',sic_wantran%nwan,nmatp,nmatp,zone,sic_wb(1,1,ispn,ikloc),&
    &sic_wantran%nwan,hm1(1,1,ispn),nmatp,zzero,zm1,sic_wantran%nwan)
  call zgemm('N','C',sic_wantran%nwan,sic_wantran%nwan,nmatp,zone,zm1,&
    &sic_wantran%nwan,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,zone,&
    &sic_wan_h0k(1,1,ikloc),sic_wantran%nwan)
enddo
deallocate(hm1,zm1)




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

allocate(zm1(sic_wantran%nwan,nmatp))
allocate(zm0(sic_wantran%nwan,sic_wantran%nwan))

zm0(:,:)=sic_wan_h0k(:,:,ikloc)+vk(:,:)

do ispn=1,nspinor

  call zgemm('N','N',sic_wantran%nwan,nmatp,sic_wantran%nwan,-zone,zm0,sic_wantran%nwan,&
    &sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,zzero,zm1,sic_wantran%nwan)
  
  !call zgemm('C','N',nmatp,nmatp,sic_wantran%nwan,zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,&
  !  &zm1,sic_wantran%nwan,zone,hunif(1,1,ispn),nmatp)


  call zgemm('C','N',sic_wantran%nwan,nmatp,sic_wantran%nwan,-zone,vk,sic_wantran%nwan,&
    &sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,zone,zm1,sic_wantran%nwan)
  
  do i1=1,nmatp
    do j1=1,sic_wantran%nwan
      zm1(j1,i1)=zm1(j1,i1)+sic_wb(j1,i1,ispn,ikloc)*(sic_wan_h0k(j1,j1,ikloc)+vk(j1,j1)) + &
        sic_wvb(j1,i1,ispn,ikloc)
    enddo
  enddo
  
  call zgemm('C','N',nmatp,nmatp,sic_wantran%nwan,zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,&
    &zm1,sic_wantran%nwan,zone,hunif(1,1,ispn),nmatp)


  
  !call zgemm('C','N',nmatp,nmatp,sic_wantran%nwan,zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,&
  !  &sic_wvb(1,1,ispn,ikloc),sic_wantran%nwan,zone,hunif(1,1,ispn),nmatp)
  
  call zgemm('C','N',nmatp,nmatp,sic_wantran%nwan,zone,sic_wvb(1,1,ispn,ikloc),sic_wantran%nwan,&
    &sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,zone,hunif(1,1,ispn),nmatp)

  !do i1=1,nmatp
  !  do j1=1,sic_wantran%nwan
  !    zm1(j1,i1)=sic_wb(j1,i1,ispn,ikloc)*(sic_wan_h0k(j1,j1,ikloc)+vk(j1,j1))
  !  enddo
  !enddo
  !
  !call zgemm('C','N',nmatp,nmatp,sic_wantran%nwan,zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,&
  !  &zm1,sic_wantran%nwan,zone,hunif(1,1,ispn),nmatp)




  !do i1=1,nmatp
  !  do i2=1,nmatp


      !do j1=1,sic_wantran%nwan
        !do j2=1,sic_wantran%nwan
         
          !hunif(i1,i2,ispn)=hunif(i1,i2,ispn)-sic_wan_h0k(j1,j2,ikloc)*&
          !  &dconjg(sic_wb(j1,i1,ispn,ikloc))*sic_wb(j2,i2,ispn,ikloc)
        
          
          !hunif(i1,i2,ispn)=hunif(i1,i2,ispn)-&
          !  &dconjg(sic_wb(j1,i1,ispn,ikloc))*vk(j1,j2)*sic_wb(j2,i2,ispn,ikloc)-&
          !  &dconjg(sic_wb(j2,i1,ispn,ikloc))*dconjg(vk(j1,j2))*sic_wb(j1,i2,ispn,ikloc)
       
        !enddo !j2
       
        !hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
        !  &dconjg(sic_wb(j1,i1,ispn,ikloc))*sic_wvb(j1,i2,ispn,ikloc)+&
        !  &dconjg(sic_wvb(j1,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
       
        !hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
        !  &dconjg(sic_wb(j1,i1,ispn,ikloc))*&
        !  &sic_wb(j1,i2,ispn,ikloc)*(sic_wan_h0k(j1,j1,ikloc)+vk(j1,j1))
     
      !enddo !j1
    
    !enddo
  !enddo
enddo


deallocate(zm1,zm0)








! compute V_{nn'}(k)
!allocate(vk(sic_wantran%nwan,sic_wantran%nwan))
!vk=zzero
!do i=1,sic_wantran%nwt
!  n1=sic_wantran%iwt(1,i)
!  j1=sic_wantran%idxiwan(n1)
!  n2=sic_wantran%iwt(2,i)
!  j2=sic_wantran%idxiwan(n2)
!  vtrl(:)=sic_wantran%iwt(3:5,i)
!  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
!  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
!  vk(j1,j2)=vk(j1,j2)+expikt*sic_vme(i)
!enddo
!! setup unified Hamiltonian
!do ispn=1,nspinor
!  do i1=1,nmatp
!    do i2=1,nmatp
!      do j1=1,sic_wantran%nwan
!        n1=sic_wantran%iwan(j1)
!        vn=dreal(sic_vme(sic_wantran%iwtidx(n1,n1,0,0,0)))
!        do j2=1,sic_wantran%nwan
!! -\sum_{\alpha,\alpha'} ( P_{\alpha}V_{\alpha}P_{\alpha'} +
!!                           P_{\alpha'}V_{\alpha}P_{\alpha} )
!          hunif(i1,i2,ispn)=hunif(i1,i2,ispn)-&
!            vk(j1,j2)*dconjg(sic_wb(j1,i1,ispn,ikloc))*&
!            sic_wb(j2,i2,ispn,ikloc)-dconjg(vk(j1,j2))*&
!            dconjg(sic_wb(j2,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
!        enddo !j2
!! \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
!        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
!          dconjg(sic_wb(j1,i1,ispn,ikloc))*&
!          sic_wb(j1,i2,ispn,ikloc)*vn
!! \sum_{\alpha} P_{\alpha}V_{\alpha}+V_{\alpha}P_{\alpha}
!        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
!          dconjg(sic_wb(j1,i1,ispn,ikloc))*sic_wvb(j1,i2,ispn,ikloc)+&
!          dconjg(sic_wvb(j1,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
!      enddo !j1
!    enddo !i2 
!  enddo !i1
!enddo !ispn
deallocate(vk)
call timer_stop(t_sic_hunif)
return
end subroutine

