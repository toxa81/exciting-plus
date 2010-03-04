subroutine calc_uuj(uuj,lmaxexp,gq0)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: lmaxexp
real(8), intent(in) :: gq0(ngvecme)
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
! local variables
integer ia,is,ias,ig,l1,l2,l3,io1,io2,ir
real(8), allocatable :: jlgq0r(:,:,:,:)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) t1,jl(0:lmaxexp)
! for parallel
integer cart_size,cart_rank,cart_group,ierr
integer tmp_comm,tmp_group,tmp_size
integer, allocatable :: ranks(:)
integer i,idx0,bs

allocate(jlgq0r(nrmtmax,0:lmaxexp,nspecies,ngvecme))

idx0=0
bs=ngvecme
#ifdef _MPI_
call mpi_comm_size(comm_cart_110,cart_size,ierr)
call mpi_comm_rank(comm_cart_110,cart_rank,ierr)
call mpi_comm_group(comm_cart_110,cart_group,ierr)
tmp_size=min(cart_size,ngvecme)
tmp_group=MPI_GROUP_EMPTY
allocate(ranks(tmp_size))
do i=1,tmp_size
  ranks(i)=i-1
enddo
call mpi_group_incl(cart_group,tmp_size,ranks,tmp_group,ierr) 
deallocate(ranks)
call mpi_comm_create(comm_cart_110,tmp_group,tmp_comm,ierr)
idx0=0
bs=0
if (cart_rank.lt.tmp_size) then
  call idxbos(ngvecme,tmp_size,cart_rank+1,idx0,bs)
endif
#endif

! generate Bessel functions j_l(|G+q'|x)
jlgq0r=0.d0
do ig=idx0+1,idx0+bs
  do is=1,nspecies
    do ir=1,nrmt(is)
      t1=gq0(ig)*spr(ir,is)
      call sbessel(lmaxexp,t1,jl)
      jlgq0r(ir,:,is,ig)=jl(:)
    enddo
  enddo
enddo

uuj=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=idx0+1,idx0+bs
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do l3=0,lmaxexp
            do io1=1,nrfmax
              do io2=1,nrfmax
                do ir=1,nrmt(is)
                  fr(ir)=urf(ir,l1,io1,ias)*urf(ir,l2,io2,ias)*jlgq0r(ir,l3,is,ig)*(spr(ir,is)**2)
                enddo
                call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
                uuj(l1,l2,l3,io1,io2,ias,ig)=gr(nrmt(is))
              enddo !io2
            enddo !io1
          enddo !l3
        enddo !l2
      enddo !l1   
    enddo !ig
  enddo !ia
enddo !is

deallocate(jlgq0r)

#ifdef _MPI_
if (cart_rank.lt.tmp_size) then
  do ig=1,ngvecme
    call dsync2(tmp_comm,uuj(0,0,0,1,1,1,ig), &
      (lmaxvr+1)*(lmaxvr+1)*(lmaxexp+1)*nrfmax*nrfmax*natmtot,.true.,.false.)
    call barrier(tmp_comm)
  enddo
endif
do ig=1,ngvecme
  call dsync2(comm_cart_110,uuj(0,0,0,1,1,1,ig), &
    (lmaxvr+1)*(lmaxvr+1)*(lmaxexp+1)*nrfmax*nrfmax*natmtot,.false.,.true.)
  call barrier(comm_cart_110)
enddo
if (cart_rank.lt.tmp_size) then
  call mpi_comm_free(tmp_comm,ierr)
  call mpi_group_free(tmp_group,ierr)
endif
call mpi_group_free(cart_group,ierr)
#endif

return
end   
