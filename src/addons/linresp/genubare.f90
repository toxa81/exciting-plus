subroutine genubare(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)
real(8), allocatable :: vcgq(:)
integer ig,ntloc,itloc,it,n1,n2
real(8) vtc(3),vgq0c(3),gq0
complex(8) expiqt

ntloc=mpi_grid_map(ntr_uscrn,dim_b)

! setup 4Pi/|G+q|^2 array
allocate(vcgq(ngvecme))
do ig=1,ngvecme
! generate G+q vectors  
  vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
  gq0=dot_product(vgq0c,vgq0c)
  if (ig.eq.1.and.ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
    vcgq(ig)=0.d0
    do n1=1,8
      vgq0c(:)=q0gamma(:,n1)
      gq0=dot_product(vgq0c,vgq0c)
      vcgq(ig)=vcgq(ig)+0.125*(2*Pi)**3*a0gamma(n1)*fourpi/gq0
    enddo
  else
    vcgq(ig)=fourpi/gq0
  endif
enddo !ig
vcgq(:)=vcgq(:)/omega/nkptnr

if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
  write(*,*)'gamma contribution to integral',vcgq(1)
endif

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
            expiqt*dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig))*vcgq(ig)* &
            megqwan(idxmegqwan(n2,n2,0,0,0),ig)
        enddo
      enddo
    enddo
  enddo
endif

deallocate(vcgq)

return
end