subroutine get_adjoint_megqblh
use modmain
use mod_expigqr
implicit none
integer ik,jk,ikstep,nkstep,jkloc,i,j,tag
logical need_to_recieve 
integer, allocatable :: jkmap(:,:)
!
if (allocated(amegqblh)) deallocate(amegqblh)
allocate(amegqblh(nstsv*nstsv,ngvecme,nkptnrloc))
if (allocated(namegqblh)) deallocate(namegqblh)
allocate(namegqblh(nkptnrloc))
if (allocated(bamegqblh)) deallocate(bamegqblh)
allocate(bamegqblh(2,nstsv*nstsv,nkptnrloc))
!
nkstep=nkptnrloc
call mpi_grid_bcast(nkstep,dims=(/dim_k/))
allocate(jkmap(2,0:mpi_grid_dim_size(dim_k)-1))
do ikstep=1,nkstep
  jkmap=-1
  need_to_recieve=.false.
  ! if this processor has a k-point for this step
  if (ikstep.le.nkptnrloc) then
    ! k-point global index
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikstep)
    ! k-q point global index
    jk=idxkq(3,ik)
    ! find index of processor and a local jk index
    jkloc=mpi_grid_map(nkptnr,dim_k,x=j,glob=jk)
    ! save index of processor from which k-q point is recieved and local index of k-q point
    jkmap(1,mpi_grid_dim_pos(dim_k))=j
    jkmap(2,mpi_grid_dim_pos(dim_k))=jkloc
    ! make a local copy if jk is on the same processor
    if (j.eq.mpi_grid_dim_pos(dim_k)) then
      amegqblh(:,:,ikstep)=megqblh(:,:,jkloc)
      namegqblh(ikstep)=nmegqblh(jkloc)
      bamegqblh(:,:,ikstep)=bmegqblh(:,:,jkloc)
    else
      need_to_recieve=.true.
    endif
  endif
  call mpi_grid_reduce(jkmap(1,0),2*mpi_grid_dim_size(dim_k),dims=(/dim_k/),all=.true.,op=op_max)
  ! check who needs k-point which is stored on this processor
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (jkmap(1,i).eq.mpi_grid_dim_pos(dim_k).and.mpi_grid_dim_pos(dim_k).ne.i) then
      jkloc=jkmap(2,i)
      ! send to proc i
      tag=(ikstep*mpi_grid_dim_size(dim_k)+i)*10
      call mpi_grid_send(megqblh(1,1,jkloc),nstsv*nstsv*ngvecme,(/dim_k/),(/i/),tag)
      call mpi_grid_send(nmegqblh(jkloc),1,(/dim_k/),(/i/),tag+1)
      call mpi_grid_send(bmegqblh(1,1,jkloc),2*nstsv*nstsv,(/dim_k/),(/i/),tag+2)
    endif
  enddo
  if (need_to_recieve) then
    j=jkmap(1,mpi_grid_dim_pos(dim_k))
    tag=(ikstep*mpi_grid_dim_size(dim_k)+mpi_grid_dim_pos(dim_k))*10
    call mpi_grid_recieve(amegqblh(1,1,ikstep),nstsv*nstsv*ngvecme,(/dim_k/),(/j/),tag)
    call mpi_grid_recieve(namegqblh(ikstep),1,(/dim_k/),(/j/),tag+1)
    call mpi_grid_recieve(bamegqblh(1,1,ikstep),2*nstsv*nstsv,(/dim_k/),(/j/),tag+2)
  endif
enddo
deallocate(jkmap)
call mpi_grid_barrier((/dim_k/))
return
end subroutine
