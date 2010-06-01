subroutine genubare(iq)
use modmain
implicit none
integer, intent(in) :: iq
integer ivq0m(3)
real(8), allocatable :: vcgq(:)
integer ig,ntloc,itloc,it,n1,n2
real(8) vtc(3),vgq0c(3),gq0,a0
complex(8) expiqt
logical l1

ivq0m(:)=ivq0m_list(:,iq)
ntloc=mpi_grid_map(ntr_uscrn,dim_b)

l1=.true.
! setup 4Pi/|G+q|^2 array
allocate(vcgq(ngvecme))
do ig=1,ngvecme
  if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
    if (ig.eq.1) then
      vgq0c(:)=vgc(:,ig+gvecme1-1)+q0gamma(:,iq)
      a0=a0gamma(iq)
    else
      vgq0c(:)=vgc(:,ig+gvecme1-1)
      a0=0.125d0
    endif      
  else
    vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
    a0=1.d0
  endif
  gq0=dot_product(vgq0c,vgq0c)
  if (gq0.gt.vhgqmax) then
    if (ig.eq.ngvecme.and.l1) then
      write(*,'("Warning(genubare) : not enough G-vectors")')
      write(*,'(" ig : ",I4)')ig
      write(*,'(" iq : ",I4)')iq        
    endif
    vcgq(iq)=0.d0
    l1=.false.
  else
    vcgq(ig)=a0*fourpi/gq0
  endif
enddo

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