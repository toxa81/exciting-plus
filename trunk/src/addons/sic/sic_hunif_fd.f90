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
integer i,j,ik,vtrl(3),n1,n2,i1,i2,j1,j2,ispn1,ispn2,istfv1,istfv2,ist1,ist2,ispn
real(8) vtrc(3),en
complex(8) expikt
complex(8), allocatable :: vk(:,:),zm1(:,:),zm2(:,:)
character*500 msg,fname
logical, parameter :: tcheckherm=.false.
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_sic_hunif)
! restore full hermitian matrix
!do ispn=1,nspinor
!  do j1=2,nmatp
!    do j2=1,j1-1
!      hunif(j1,j2,ispn)=dconjg(hunif(j2,j1,ispn))
!    enddo
!  enddo
!  do j1=1,nmatp
!    hunif(j1,j1,ispn)=zone*dreal(hunif(j1,j1,ispn))
!  enddo
!enddo
!do j1=2,nmatp
!  do j2=1,j1-1
!    om(j1,j2)=dconjg(om(j2,j1))
!  enddo
!enddo
!do j1=1,nmatp
!  om(j1,j1)=zone*dreal(om(j1,j1))
!enddo
! compute H_{nn'}^{0}(k); remember that on input hunif=H0
!allocate(zm1(sic_wantran%nwan,nmatp))
!allocate(zm2(sic_wantran%nwan,sic_wantran%nwan))
!sic_wan_h0k(:,:,ikloc)=zzero
!do ispn=1,nspinor
!  call zgemm('N','N',sic_wantran%nwan,nmatp,nmatp,zone,sic_wb(1,1,ispn,ikloc),&
!    sic_wantran%nwan,hunif(1,1,ispn),nmatp,zzero,zm1,sic_wantran%nwan)
!  call zgemm('N','C',sic_wantran%nwan,sic_wantran%nwan,nmatp,zone,zm1,&
!    sic_wantran%nwan,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan,zone,&
!    sic_wan_h0k(1,1,ikloc),sic_wantran%nwan)
!enddo

!call invzge(om,nmatp)
!
!  call zgemm('N','N',sic_wantran%nwan,nmatp,nmatp,zone,sic_wb(1,1,1,ikloc),&
!    sic_wantran%nwan,om,nmatp,zzero,zm1,sic_wantran%nwan)
!  call zgemm('N','C',sic_wantran%nwan,sic_wantran%nwan,nmatp,zone,zm1,&
!    sic_wantran%nwan,sic_wb(1,1,1,ikloc),sic_wantran%nwan,zzero,&
!    zm2,sic_wantran%nwan)
!
!do j1=1,sic_wantran%nwan
!  do j2=1,sic_wantran%nwan
!    write(*,*)j1,j2,zm2(j1,j2)
!  enddo
!enddo
!deallocate(zm1)
!deallocate(zm2)
!call pstop

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
!
!!if (mpi_grid_root((/dim2/))) then
!!  write(fname,'("hlda_n",I2.2,"_k",I4.4".txt")')nproc,ik
!!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!!  write(fname,'("sic_wann_h0k_n",I2.2,"_k",I4.4".txt")')nproc,ik
!!  call wrmtrx(fname,nwantot,nwantot,sic_wann_h0k,nwantot)
!!  write(fname,'("sic_vwank_np",I2.2,"_kp",I4.4".txt")')nproc,ik
!!  call wrmtrx(fname,nwantot,nwantot,vwank,nwantot)
!!endif
!
! setup unified Hamiltonian
! 1-st term: LDA Hamiltonian itself
do ispn=1,nspinor
  do i1=1,nmatp
    do i2=1,nmatp
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
          hunif(i1,i2,ispn)=hunif(i1,i2,ispn)-&
              vk(j1,j2)*dconjg(sic_wb(j1,i1,ispn,ikloc))*&
              sic_wb(j2,i2,ispn,ikloc)-dconjg(vk(j1,j2))*&
              dconjg(sic_wb(j2,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
        enddo !j2
! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
          dconjg(sic_wb(j1,i1,ispn,ikloc))*&
          sic_wb(j1,i2,ispn,ikloc)*en
! 5-th term, second part: \sum_{\alpha} P_{\alpha}V_{\alpha}+V_{\alpha}P_{\alpha}
        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
          dconjg(sic_wb(j1,i1,ispn,ikloc))*sic_wvb(j1,i2,ispn,ikloc)+&
          dconjg(sic_wvb(j1,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
      enddo !j1
    enddo !i2 
  enddo !i1
enddo !ispn
!
!!if (mpi_grid_root((/dim2/))) then
!!  write(fname,'("hunif_n",I2.2,"_k",I4.4".txt")')nproc,ik
!!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!!endif
!
!if (tcheckherm) then
!  call checkherm(nstsv,hunif,i,j)
!  if (i.gt.0) then
!    write(msg,'("matrix hunif is not Hermitian at k-point ",I4,&
!      &" : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
!      ik,i,j,dreal(hunif(i,j)),dimag(hunif(i,j)),dreal(hunif(j,i)),dimag(hunif(j,i))
!    call mpi_grid_msg("sic_hunif",msg)
!  endif
!endif
deallocate(vk)
call timer_stop(t_sic_hunif)
return
end subroutine

