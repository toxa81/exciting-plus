subroutine sic_hunif(ikloc,hunif)
use modmain
use mod_sic
use mod_hdf5
use mod_wannier
use mod_linresp
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: hunif(nstsv,nstsv)
! local variables
integer i,j,ik,vtrl(3),n1,n2,j1,j2,ispn1,ispn2,istfv1,istfv2,ist1,ist2
real(8) vtrc(3),en
complex(8) expikt
complex(8), allocatable :: vk(:,:),zm1(:,:)
!complex(8), allocatable :: wfmt1(:,:),wfmt2(:,:)
character*500 msg,fname
logical, parameter :: tcheckherm=.false.
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_sic_hunif)
! restore full hermitian matrix
do j1=2,nstsv
  do j2=1,j1-1
    hunif(j1,j2)=dconjg(hunif(j2,j1))
  enddo
enddo
do j1=1,nstsv
  hunif(j1,j1)=zone*dreal(hunif(j1,j1))
enddo
! compute H_{nn'}^{0}(k); remember that on input hunif=H0
allocate(zm1(sic_wantran%nwan,nstsv))
call zgemm('N','N',sic_wantran%nwan,nstsv,nstsv,zone,sic_wb(1,1,1,ikloc),&
  &sic_wantran%nwan,hunif,nstsv,zzero,zm1,sic_wantran%nwan)
call zgemm('N','C',sic_wantran%nwan,sic_wantran%nwan,nstsv,zone,zm1,&
  &sic_wantran%nwan,sic_wb(1,1,1,ikloc),sic_wantran%nwan,zzero,&
  &sic_wan_h0k(1,1,ikloc),sic_wantran%nwan)
deallocate(zm1)
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

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hlda_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!  write(fname,'("sic_wann_h0k_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwantot,nwantot,sic_wann_h0k,nwantot)
!  write(fname,'("sic_vwank_np",I2.2,"_kp",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwantot,nwantot,vwank,nwantot)
!endif

! setup unified Hamiltonian
! 1-st term: LDA Hamiltonian itself
do istfv1=1,nstfv
  do ispn1=1,nspinor
    ist1=istfv1+(ispn1-1)*nstfv
    do istfv2=1,nstfv
      do ispn2=1,nspinor
        ist2=istfv2+(ispn2-1)*nstfv
        do j1=1,sic_wantran%nwan
          !n1=sic_wantran%iwan(j1)
          !en=dreal(sic_vme(sic_wantran%iwtidx(n1,n1,0,0,0)))!+sic_wan_e0(n1)
          do j2=1,sic_wantran%nwan
! 2-nd term : -\sum_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'}     
            hunif(ist1,ist2)=hunif(ist1,ist2)-sic_wan_h0k(j1,j2,ikloc)*&
              dconjg(sic_wb(j1,istfv1,ispn1,ikloc))*sic_wb(j2,istfv2,ispn2,ikloc)
! 5-th term, first part: 
!  -\sum_{\alpha,\alpha'} ( P_{\alpha}V_{\alpha}P_{\alpha'} +
!                           P_{\alpha'}V_{\alpha}P_{\alpha} )
            !hunif(ist1,ist2)=hunif(ist1,ist2)-&
            !  vk(j1,j2)*dconjg(sic_wb(j1,istfv1,ispn1,ikloc))*&
            !  sic_wb(j2,istfv2,ispn2,ikloc)-dconjg(vk(j1,j2))*&
            !  dconjg(sic_wb(j2,istfv1,ispn1,ikloc))*sic_wb(j1,istfv2,ispn2,ikloc)
          enddo !j2
! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
          hunif(ist1,ist2)=hunif(ist1,ist2)+&
            dconjg(sic_wb(j1,istfv1,ispn1,ikloc))*&
            sic_wb(j1,istfv2,ispn2,ikloc)*(sic_wan_h0k(j1,j1,ikloc)+vk(j1,j1))
! 5-th term, second part: \sum_{\alpha} P_{\alpha}V_{\alpha}+V_{\alpha}P_{\alpha}
          !hunif(ist1,ist2)=hunif(ist1,ist2)+&
          !  dconjg(sic_wb(j1,istfv1,ispn1,ikloc))*sic_wvb(j1,istfv2,ispn2,ikloc)+&
          !  dconjg(sic_wvb(j1,istfv1,ispn1,ikloc))*sic_wb(j1,istfv2,ispn2,ikloc)
        enddo !j1
      enddo !ispn2
    enddo ! istfv2
  enddo !ispn1
enddo !istfv1

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hunif_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!endif

if (tcheckherm) then
  call checkherm(nstsv,hunif,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif is not Hermitian at k-point ",I4,&
      &" : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(hunif(i,j)),dimag(hunif(i,j)),dreal(hunif(j,i)),dimag(hunif(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif
deallocate(vk)
call timer_stop(t_sic_hunif)
return
end
