subroutine sic_hunif(ikloc,hunif)
use modmain
use mod_sic
use mod_hdf5
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: hunif(nstsv,nstsv)
! local variables
logical exist
integer i,j,n,ik,vtrl(3),n1,j1,j2,ispn1,ispn2,jst1,jst2,n2
real(8) vtrc(3),t1
complex(8) expikt
complex(8), allocatable :: vwank(:,:)
complex(8), allocatable :: zm1(:,:)
character*500 fname,msg
logical, parameter :: tcheckherm=.true.

! simplified potential correction if we don't have SIC potential yet
!if (.not.tsic_wv) then
!  do j1=1,nstfv
!    do ispn1=1,nspinor
!      jst1=j1+(ispn1-1)*nstfv
!      do j2=1,nstfv
!        do ispn2=1,nspinor
!          jst2=j2+(ispn2-1)*nstfv
!          do n=1,nwann
!            hunif(jst1,jst2)=hunif(jst1,jst2)+dconjg(sic_wb(n,j1,ispn1,ikloc))*&
!              sic_wb(n,j2,ispn2,ikloc)*dreal(vwanme(idxmegqwan(n,n,0,0,0)))
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!  return
!endif
! restore full hermitian matrix
do j1=2,nstsv
  do j2=1,j1-1
    hunif(j1,j2)=dconjg(hunif(j2,j1))
  enddo
enddo
do j1=1,nstsv
  hunif(j1,j1)=zone*dreal(hunif(j1,j1))
enddo

! compute V_{nn'}(k)
allocate(vwank(nwann,nwann))
vwank=zzero
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
do i=1,nmegqwan
  n1=imegqwan(1,i)
  n2=imegqwan(2,i)
  vtrl(:)=imegqwan(3:5,i)
  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  vwank(n1,n2)=vwank(n1,n2)+expikt*vwanme(i)
enddo
! compute H_{nn'}^{0}(k); remember that on input hunif=H0
allocate(zm1(nwann,nstsv))
call zgemm('N','N',nwann,nstsv,nstsv,zone,sic_wb(1,1,1,ikloc),nwann,hunif,&
  nstsv,zzero,zm1,nwann)
call zgemm('N','C',nwann,nwann,nstsv,zone,zm1,nwann,sic_wb(1,1,1,ikloc),nwann,&
  zzero,sic_wann_h0k(1,1,ikloc),nwann)
deallocate(zm1)
!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hlda_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!  write(fname,'("sic_wann_h0k_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nwann,sic_wann_h0k,nwann)
!  write(fname,'("sic_vwank_np",I2.2,"_kp",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nwann,vwank,nwann)
!endif

!allocate(zm1(nstsv,nstsv))
!zm1=zzero
!do j1=1,nstfv
!  do ispn1=1,nspinor
!    jst1=j1+(ispn1-1)*nstfv
!    do j2=1,nstfv
!      do ispn2=1,nspinor
!        jst2=j2+(ispn2-1)*nstfv
!        do n=1,nwann
!          zm1(jst1,jst2)=zm1(jst1,jst2)+dconjg(a(n,j1,ispn1))*&
!            a(n,j2,ispn2)*(vn(n)+sic_wann_e0(n))
!        enddo
!      enddo
!    enddo
!  enddo
!enddo
!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hunif2_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,zm1,nstsv)
!endif

if (tcheckherm) then
  call checkherm(nwann,sic_wann_h0k(1,1,ikloc),i,j)
  if (i.gt.0) then
    write(msg,'("matrix sic_wann_h0k is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(sic_wann_h0k(i,j,ikloc)),dimag(sic_wann_h0k(i,j,ikloc)),&
      dreal(sic_wann_h0k(j,i,ikloc)),dimag(sic_wann_h0k(j,i,ikloc))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif

allocate(zm1(nstsv,nstsv))
! setup unified Hamiltonian
! 1-st term: LDA Hamiltonian itself
! 2-nd term : -\sum_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'}
zm1=zzero
do j1=1,nstfv
  do ispn1=1,nspinor
    jst1=j1+(ispn1-1)*nstfv
    do j2=1,nstfv
      do ispn2=1,nspinor
        jst2=j2+(ispn2-1)*nstfv
        do n1=1,nwann
          do n2=1,nwann
            zm1(jst1,jst2)=zm1(jst1,jst2)-sic_wann_h0k(n1,n2,ikloc)*&
              dconjg(sic_wb(n1,j1,ispn1,ikloc))*sic_wb(n2,j2,ispn2,ikloc)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
if (tcheckherm) then
  call checkherm(nstsv,zm1,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif2 is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(zm1(i,j)),dimag(zm1(i,j)),dreal(zm1(j,i)),dimag(zm1(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif
hunif=hunif+zm1
! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
zm1=zzero
do j1=1,nstfv
  do ispn1=1,nspinor
    jst1=j1+(ispn1-1)*nstfv
    do j2=1,nstfv
      do ispn2=1,nspinor
        jst2=j2+(ispn2-1)*nstfv
        do n=1,nwann
          zm1(jst1,jst2)=zm1(jst1,jst2)+dconjg(sic_wb(n,j1,ispn1,ikloc))*&
            sic_wb(n,j2,ispn2,ikloc)*(dreal(vwanme(idxmegqwan(n,n,0,0,0)))+&
            sic_wann_e0(n))
        enddo
      enddo
    enddo
  enddo
enddo
if (tcheckherm) then
  call checkherm(nstsv,zm1,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif34 is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(zm1(i,j)),dimag(zm1(i,j)),dreal(zm1(j,i)),dimag(zm1(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif
hunif=hunif+zm1
! 5-th term : \sum_{\alpha} P_{\alpha} V_{\alpha} Q + 
!             \sum_{\alpha} Q V_{\alpha} P_{\alpha} 
!  where Q=1-\sum_{\alpha'}P_{\alpha'}
zm1=zzero
do j1=1,nstfv
  do ispn1=1,nspinor
    jst1=j1+(ispn1-1)*nstfv
    do j2=1,nstfv
      do ispn2=1,nspinor
        jst2=j2+(ispn2-1)*nstfv
        do n=1,nwann
          zm1(jst1,jst2)=zm1(jst1,jst2)+&
            dconjg(sic_wb(n,j1,ispn1,ikloc))*sic_wvb(n,j2,ispn2,ikloc)+&
            dconjg(sic_wvb(n,j1,ispn1,ikloc))*sic_wb(n,j2,ispn2,ikloc)
        enddo
      enddo
    enddo
  enddo
enddo
if (tcheckherm) then
  call checkherm(nstsv,zm1,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif51 is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(zm1(i,j)),dimag(zm1(i,j)),dreal(zm1(j,i)),dimag(zm1(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif
hunif=hunif+zm1
zm1=zzero
do j1=1,nstfv
  do ispn1=1,nspinor
    jst1=j1+(ispn1-1)*nstfv
    do j2=1,nstfv
      do ispn2=1,nspinor
        jst2=j2+(ispn2-1)*nstfv
        do n1=1,nwann
          do n2=1,nwann
            hunif(jst1,jst2)=hunif(jst1,jst2)-&
              vwank(n1,n2)*dconjg(sic_wb(n1,j1,ispn1,ikloc))*&
              sic_wb(n2,j2,ispn2,ikloc)-dconjg(vwank(n1,n2))*&
              dconjg(sic_wb(n2,j1,ispn1,ikloc))*sic_wb(n1,j2,ispn2,ikloc)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
if (tcheckherm) then
  call checkherm(nstsv,zm1,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif52 is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(zm1(i,j)),dimag(zm1(i,j)),dreal(zm1(j,i)),dimag(zm1(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif
hunif=hunif+zm1

!
!
!do j1=1,nstfv
!  do ispn1=1,nspinor
!    jst1=j1+(ispn1-1)*nstfv
!    do j2=1,nstfv
!      do ispn2=1,nspinor
!        jst2=j2+(ispn2-1)*nstfv
!        do n1=1,nwann
!          do n2=1,nwann
!            hunif(jst1,jst2)=hunif(jst1,jst2)-sic_wann_h0k(n1,n2,ikloc)*&
!              dconjg(sic_wb(n1,j1,ispn1,ikloc))*sic_wb(n2,j2,ispn2,ikloc)
!          enddo
!        enddo
!! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
!! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
!        do n=1,nwann
!          hunif(jst1,jst2)=hunif(jst1,jst2)+dconjg(sic_wb(n,j1,ispn1,ikloc))*&
!            sic_wb(n,j2,ispn2,ikloc)*(dreal(vwanme(idxmegqwan(n,n,0,0,0)))+&
!            sic_wann_e0(n))
!        enddo
!! 5-th term : \sum_{\alpha} P_{\alpha} V_{\alpha} Q + 
!!             \sum_{\alpha} Q V_{\alpha} P_{\alpha} 
!!  where Q=1-\sum_{\alpha'}P_{\alpha'}
!        do n=1,nwann
!          hunif(jst1,jst2)=hunif(jst1,jst2)+&
!            dconjg(sic_wb(n,j1,ispn1,ikloc))*sic_wvb(n,j2,ispn2,ikloc)+&
!            dconjg(sic_wvb(n,j1,ispn1,ikloc))*sic_wb(n,j2,ispn2,ikloc)
!        enddo
!        do n1=1,nwann
!          do n2=1,nwann
!            hunif(jst1,jst2)=hunif(jst1,jst2)-&
!              vwank(n1,n2)*dconjg(sic_wb(n1,j1,ispn1,ikloc))*&
!              sic_wb(n2,j2,ispn2,ikloc)-dconjg(vwank(n1,n2))*&
!              dconjg(sic_wb(n2,j1,ispn1,ikloc))*sic_wb(n1,j2,ispn2,ikloc)
!          enddo
!        enddo
!      enddo !ispn2
!    enddo !j2
!  enddo !ispn1
!enddo !j1

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hunif_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!endif

if (tcheckherm) then
  call checkherm(nstsv,hunif,i,j)
  if (i.gt.0) then
    write(msg,'("matrix hunif is not Hermitian at k-point ",I4,&
      " : i, j, m(i,j), m(j,i)=",I5,",",I5,", (",2G18.10,"), (",2G18.10,")")')&
      ik,i,j,dreal(hunif(i,j)),dimag(hunif(i,j)),dreal(hunif(j,i)),dimag(hunif(j,i))
    call mpi_grid_msg("sic_hunif",msg)
  endif
endif

! check hermiticity
!if (mpi_grid_root((/dim2/))) then  
!  do j1=1,nstsv
!    do j2=1,nstsv
!      t1=abs(hunif(j1,j2)-dconjg(hunif(j2,j1)))
!      if (t1.gt.epsherm) then
!        write(*,*)
!        write(*,'("Warning(sic_hunif) : unified Hamiltonian is not hermitian")')
!        write(*,'("  k-point : ",I4)')ik
!        write(*,'("  j1, j2, diff : ",2I4,G18.10)')j1,j2,t1
!      endif
!    enddo
!  enddo
!endif
deallocate(vwank)
return
end