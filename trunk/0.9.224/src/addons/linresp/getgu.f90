subroutine getgu(req,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: lmaxexp
real(8), intent(in) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
complex(8), intent(in) :: ylmgq0((lmaxexp+1)**2,ngvecme)
complex(8), intent(in) :: sfacgq0(ngvecme,natmtot)
integer, intent(inout) :: ngumax
integer, intent(out) :: ngu(natmtot,ngvecme)
complex(4), intent(out) :: gu(ngumax,natmtot,ngvecme)
integer, intent(out) :: igu(4,ngumax,natmtot,ngvecme)

integer ig,ias,i,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3
real(8) t1
! for parallel
integer cart_size,cart_rank,cart_group,ierr
integer tmp_comm,tmp_group,tmp_size
integer, allocatable :: ranks(:)
integer idx0,bs

real(8), external :: gaunt

idx0=0
bs=ngvecme
#ifdef _MPI_
call mpi_comm_size(comm_cart,cart_size,ierr)
call mpi_comm_rank(comm_cart,cart_rank,ierr)
call mpi_comm_group(comm_cart,cart_group,ierr)
tmp_size=min(cart_size,ngvecme)
tmp_group=MPI_GROUP_EMPTY
allocate(ranks(tmp_size))
do i=1,tmp_size
  ranks(i)=i-1
enddo
call mpi_group_incl(cart_group,tmp_size,ranks,tmp_group,ierr) 
deallocate(ranks)
call mpi_comm_create(comm_cart,tmp_group,tmp_comm,ierr)
idx0=0
bs=0
if (cart_rank.lt.tmp_size) then
  call idxbos(ngvecme,tmp_size,cart_rank+1,idx0,bs)
endif
#endif

if (req) ngumax=0 
if (.not.req) then
  igu=0
  ngu=0
  gu=cmplx(0.0,0.0)
endif
do ig=idx0+1,idx0+bs
  do ias=1,natmtot
    i=0
    do io1=1,nrfmax
      do io2=1,nrfmax
        do l1=0,lmaxvr
        do m1=-l1,l1 
          lm1=idxlm(l1,m1)
          do l2=0,lmaxvr
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            do l3=0,lmaxexp
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              t1=gaunt(l2,l1,l3,m2,m1,m3)*uuj(l1,l2,l3,io1,io2,ias,ig)
              if (abs(t1).gt.1d-10) then
                i=i+1
                if (.not.req) then
                  gu(i,ias,ig)=t1*ylmgq0(lm3,ig)*dconjg(zi**l3)*fourpi*dconjg(sfacgq0(ig,ias))
                  igu(1,i,ias,ig)=lm1
                  igu(2,i,ias,ig)=lm2
                  igu(3,i,ias,ig)=io1
                  igu(4,i,ias,ig)=io2
                endif
              endif
            enddo
            enddo
          enddo
          enddo
        enddo
        enddo
      enddo
    enddo
    if (.not.req) ngu(ias,ig)=i
    if (req) ngumax=max(ngumax,i)
  enddo
enddo

#ifdef _MPI_
if (req) then
  if (cart_rank.lt.tmp_size) then 
    call mpi_reduce(ngumax,i,1,MPI_INTEGER,MPI_MAX,0,tmp_comm,ierr)
    ngumax=i
  endif
  call mpi_bcast(ngumax,1,MPI_INTEGER,0,comm_cart,ierr)
else
  if (cart_rank.lt.tmp_size) then
    do ig=1,ngvecme
      call csync2(tmp_comm,gu(1,1,ig),ngumax*natmtot,.true.,.false.)
      call isync2(tmp_comm,igu(1,1,1,ig),4*ngumax*natmtot,.true.,.false.)
    enddo
    call isync2(tmp_comm,ngu,natmtot*ngvecme,.true.,.false.)
  endif
  do ig=1,ngvecme
    call csync2(comm_cart,gu(1,1,ig),ngumax*natmtot,.false.,.true.)
    call isync2(comm_cart,igu(1,1,1,ig),4*ngumax*natmtot,.false.,.true.)
  enddo
  call isync2(comm_cart,ngu,natmtot*ngvecme,.false.,.true.)
endif
if (cart_rank.lt.tmp_size) then
  call mpi_comm_free(tmp_comm,ierr)
  call mpi_group_free(tmp_group,ierr)
endif
call mpi_group_free(cart_group,ierr)
#endif

return
end
