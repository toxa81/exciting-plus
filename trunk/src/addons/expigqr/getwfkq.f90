subroutine getwfkq(ikstep,ngknr_jk,igkignr_jk,wfsvmt_jk,wfsvit_jk)
use modmain
use mod_nrkp
use mod_wannier
use mod_expigqr
implicit none
integer, intent(in) :: ikstep
integer, intent(out) :: ngknr_jk
integer, intent(out) :: igkignr_jk(ngkmax)
complex(8), intent(out) :: wfsvmt_jk(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfsvit_jk(ngkmax,nspinor,nstsv)

integer i,ik,jk,nkptnrloc1,jkloc,j,tag

! each proc knows that it needs wave-functions at jk=idxkq(1,ik) (at k'=k+q)
!
! the distribution of k-points could look like this
!                p0          p1          p2
!          +-----------+-----------+-----------+
! ikstep=1 | ik=1 jk=3 | ik=4 jk=2 | ik=7 jk=5 |
! ikstep=2 | ik=2 jk=4 | ik=5 jk=7 | ik=8 jk=6 |
! ikstep=3 | ik=3 jk=1 | ik=6 jk=8 |  -        |
!          +-----------+-----------+-----------+

! two actions are required:
! 1) each processor scans trough other processors and determines, which
!    processors require its part of k-points; during this phase it
!    executes non-blocking 'send'
! 2) each processor must know the index of other processor, from which 
!    it gets jk-point; during this phase it executes blocking 'recieve'

do i=0,mpi_grid_dim_size(dim_k)-1
! number of k-points on the processor i
  nkptnrloc1=mpi_grid_map(nkptnr,dim_k,x=i)
  if (ikstep.le.nkptnrloc1) then
! for the step ikstep processor i computes matrix elements between k-point ik 
    ik=mpi_grid_map(nkptnr,dim_k,x=i,loc=ikstep)
! and k-point jk
    jk=idxkq(1,ik)
! find the processor j and local index of k-point jkloc for the k-point jk
    jkloc=mpi_grid_map(nkptnr,dim_k,glob=jk,x=j)
    if (mpi_grid_dim_pos(dim_k).eq.j.and.mpi_grid_dim_pos(dim_k).ne.i) then
! send to i
      tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
      call mpi_grid_send(wfsvmtnrloc(1,1,1,1,1,jkloc),&
        lmmaxvr*nufrmax*natmtot*nspinor*nstsv,(/dim_k/),(/i/),tag)
      call mpi_grid_send(wfsvitnrloc(1,1,1,jkloc),ngkmax*nspinor*nstsv,&
        (/dim_k/),(/i/),tag+1)
      call mpi_grid_send(ngknr(jkloc),1,(/dim_k/),(/i/),tag+2)
      call mpi_grid_send(igkignr(1,jkloc),ngkmax,(/dim_k/),(/i/),tag+3)
      if (wannier_megq) then
        call mpi_grid_send(wanncnrloc(1,1,jkloc),nwantot*nstsv,(/dim_k/),(/i/),tag+4)
      endif
    endif
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      if (j.ne.i) then
! recieve from j
        tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
        call mpi_grid_recieve(wfsvmt_jk(1,1,1,1,1),&
          lmmaxvr*nufrmax*natmtot*nspinor*nstsv,(/dim_k/),(/j/),tag)
        call mpi_grid_recieve(wfsvit_jk(1,1,1),ngkmax*nspinor*nstsv,&
          (/dim_k/),(/j/),tag+1)
        call mpi_grid_recieve(ngknr_jk,1,(/dim_k/),(/j/),tag+2)
        call mpi_grid_recieve(igkignr_jk(1),ngkmax,(/dim_k/),(/j/),tag+3)
        if (wannier_megq) then
          call mpi_grid_recieve(wann_c_jk(1,1,ikstep),nwantot*nstsv,&
            (/dim_k/),(/j/),tag+4)
        endif
      else
! local copy
        wfsvmt_jk(:,:,:,:,:)=wfsvmtnrloc(:,:,:,:,:,jkloc)
        wfsvit_jk(:,:,:)=wfsvitnrloc(:,:,:,jkloc)
        ngknr_jk=ngknr(jkloc)
        igkignr_jk(:)=igkignr(:,jkloc)
        if (wannier_megq) wann_c_jk(:,:,ikstep)=wanncnrloc(:,:,jkloc)
      endif
    endif
  endif   
enddo
call mpi_grid_barrier((/dim_k/))
return
end
