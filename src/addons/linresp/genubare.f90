subroutine genubare(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)
real(8), allocatable :: vcgq(:)
integer ig,ntloc,itloc,it,n1,n2
real(8) vtc(3),vgq0c(3),gq0
complex(8) expiqt

ntloc=mpi_grid_map(ntr_uscrn,dim_b)

! setup sqrt(4Pi)/|G+q| array
allocate(vcgq(ngvecme))
do ig=1,ngvecme
! generate G+q vectors  
  vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
  gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
  if (ig.eq.1.and.ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
    vcgq(ig)=0.d0
  else
    vcgq(ig)=sqrt(fourpi)/gq0
  endif
enddo !ig

if (mpi_grid_x(dim_k).eq.0) then
  do itloc=1,ntloc
    it=mpi_grid_map(ntr_uscrn,dim_b,loc=itloc)
    vtc(:)=vtl_uscrn(1,it)*avec(:,1)+vtl_uscrn(2,it)*avec(:,2)+&
           vtl_uscrn(3,it)*avec(:,3)
    expiqt=exp(-zi*dot_product(vq0c,vtc))
    do n1=1,nwann
      do n2=1,nwann
        do ig=1,ngvecme
          ubarewan(n1,n2,it)=ubarewan(n1,n2,it)+&
            expiqt*dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig))*(vcgq(ig)**2)*&
            megqwan(idxmegqwan(n2,n2,0,0,0),ig)
        enddo
      enddo
    enddo
  enddo
endif

deallocate(vcgq)

return
end