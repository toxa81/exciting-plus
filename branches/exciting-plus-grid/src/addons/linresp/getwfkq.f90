subroutine getwfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
  wfsvit2,ngknr2,igkignr2,wann_c2) !,evecfv2,evecsv2)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: ikstep
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,*)
complex(8), intent(in) :: wfsvitloc(ngkmax,nspinor,nstsv,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
complex(8), intent(out) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfsvit2(ngkmax,nspinor,nstsv)
integer, intent(out) :: ngknr2
integer, intent(out) :: igkignr2(ngkmax)
complex(8), intent(out) :: wann_c2(nwann,nstsv)
!complex(8), intent(out) :: evecfv2(nmatmax,nstfv,nspnfv)
!complex(8), intent(out) :: evecsv2(nstsv,nstsv)

integer idx0,bs,ip,jx
logical lsend,lrecv
integer i,ik,jk,ikloc,nkptnrloc1,jkloc,j
integer rank,tag,ierr,req
integer, allocatable :: stat(:)

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

do i=0,mpi_grid_size(dim_k)-1
! number of k-points on the processor i
  nkptnrloc1=mpi_grid_map(nkptnr,dim_k,x=i)
  if (ikstep.le.nkptnrloc1) then
! for the step ikstep processor i computes matrix elements between k-point ik 
    ik=mpi_grid_map(nkptnr,dim_k,x=i,loc=ikstep)
! and k-point jk
    jk=idxkq(1,ik)
! find the processor j and local index of k-point jkloc for the k-point jk
    jkloc=mpi_grid_map(nkptnr,dim_k,glob=jk,x=j)
    if (mpi_grid_x(dim_k).eq.j.and.mpi_grid_x(dim_k).ne.i) then
! send to i
      tag=(ikstep*mpi_grid_size(dim_k)+i)*10
      call mpi_grid_send(wfsvmtloc(1,1,1,1,1,jkloc),&
        lmmaxvr*nrfmax*natmtot*nspinor*nstsv,(/dim_k/),(/i/),tag)
!      write(*,*)mpi_grid_x(dim_k), ' --> ',i,sum(wfsvmtloc(:,:,:,:,:,jkloc))
      call mpi_grid_send(wfsvitloc(1,1,1,jkloc),ngkmax*nspinor*nstsv,&
        (/dim_k/),(/i/),tag+1)
      call mpi_grid_send(ngknr(jkloc),1,(/dim_k/),(/i/),tag+2)
      call mpi_grid_send(igkignr(1,jkloc),ngkmax,(/dim_k/),(/i/),tag+3)
      if (wannier) then
        call mpi_grid_send(wann_c(1,1,jkloc),nwann*nstsv,(/dim_k/),(/i/),tag+4)
      endif
    endif
    if (mpi_grid_x(dim_k).eq.i) then
      if (j.ne.i) then
! recieve from j
        tag=(ikstep*mpi_grid_size(dim_k)+i)*10
        call mpi_grid_recieve(wfsvmt2(1,1,1,1,1),&
          lmmaxvr*nrfmax*natmtot*nspinor*nstsv,(/dim_k/),(/j/),tag)
!        write(*,*)mpi_grid_x(dim_k), ' <-- ',j,sum(wfsvmt2(:,:,:,:,:))
        call mpi_grid_recieve(wfsvit2(1,1,1),ngkmax*nspinor*nstsv,&
          (/dim_k/),(/j/),tag+1)
        call mpi_grid_recieve(ngknr2,1,(/dim_k/),(/j/),tag+2)
        call mpi_grid_recieve(igkignr2(1),ngkmax,(/dim_k/),(/j/),tag+3)
        if (wannier) then
          call mpi_grid_recieve(wann_c2(1,1),nwann*nstsv,(/dim_k/),(/j/),tag+4)
        endif
      else
! local copy
        wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,jkloc)
        wfsvit2(:,:,:)=wfsvitloc(:,:,:,jkloc)
        ngknr2=ngknr(jkloc)
        igkignr2(:)=igkignr(:,jkloc)
        if (wannier) wann_c2(:,:)=wann_c(:,:,jkloc)
      endif
    endif
  endif   
enddo
call mpi_grid_barrier((/dim_k,dim_g/))
return
end
