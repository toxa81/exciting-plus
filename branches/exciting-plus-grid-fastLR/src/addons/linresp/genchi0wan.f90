subroutine genchi0wan(igq0,mexp,chi0wan_k,chi0wan)
use modmain
implicit none
integer, intent(in) :: igq0
complex(8), intent(in) :: chi0wan_k(nmegqwan,nmegqwan,nkptnrloc)
complex(8), intent(in) :: mexp(nmegqwan,nmegqwan,nkptnrloc) !nkptnrloc,ntrchi0wan)
complex(8), intent(out) :: chi0wan(nmegqwan,nmegqwan)
integer ikloc,ik,i1,i2,n1,n2,n3,n4
integer it1(3),it2(3),it(3)
real(8) vtrc(3)
complex(8) zt1
!complex(8), allocatable :: mexp(:,:)


chi0wan(:,:)=zzero
do ikloc=1,nkptnrloc
  chi0wan(:,:)=chi0wan(:,:)+mexp(:,:,ikloc)*chi0wan_k(:,:,ikloc)
enddo
!do i1=1,nmegqwan
!  do i2=1,nmegqwan
!    !n1=imegqwan(1,i1)
!    !n2=imegqwan(2,i1)
!    it1(:)=imegqwan(3:5,i1)
!    !n3=imegqwan(1,i2)
!    !n4=imegqwan(2,i2)
!    it2(:)=imegqwan(3:5,i2)
!    it(:)=it1(:)-it2(:)
!    vtrc(:)=avec(:,1)*it(1)+avec(:,2)*it(2)+avec(:,3)*it(3)
!    do ikloc=1,nkptnrloc
!      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  ! phase e^{i(k+q)T}
!      zt1=exp(dcmplx(0.d0,dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
!  ! chi0wan=chi0wan+e^{i(k+q)T}*chi0wan_k(k)
!      chi0wan(i1,i2)=chi0wan(i1,i2)+zt1*chi0wan_k(i1,i2,ikloc)
!    enddo !ikloc
!  enddo
!enddo
    
    
    
    


!call zgemm('N','N',nmegqwan*nmegqwan,ntrchi0wan,nkptnrloc,zone,chi0wan_k,&
!  nmegqwan*nmegqwan,mexp,nkptnrloc,zzero,chi0wan,nmegqwan*nmegqwan)


!chi0wan(:,:,:)=zzero
! todo: utilize second dimension
! loop over translations
!do it2=1,ntrchi0wan
!! translation vector
!  vtrc(:)=avec(:,1)*itrchi0wan(1,it2)+&
!          avec(:,2)*itrchi0wan(2,it2)+&
!          avec(:,3)*itrchi0wan(3,it2)
!  do ikloc=1,nkptnrloc
!    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!! phase e^{i(k+q)T}
!    zt1=exp(dcmplx(0.d0,dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
!! chi0wan=chi0wan+e^{i(k+q)T}*chi0wan_k(k)
!    call zaxpy(nmegqwan*nmegqwan,zt1,chi0wan_k(1,1,ikloc),1,chi0wan(1,1,it2),1)
!  enddo !ikloc
!enddo !it2
! sum chi0wan over all k-points
call mpi_grid_reduce(chi0wan(1,1),nmegqwan*nmegqwan,dims=(/dim_k/))
chi0wan(:,:)=chi0wan(:,:)/nkptnr/omega
! compute chi0_GqGq using the Wannier functions expansion
!  chi0=\sum_{T,T'}\sum_{n,m,n',m'} A^{*}_{nmT}(q,Gq)*chi0wan(n,m,n',m',T-T')*A_{n'm'T'}(q,Gq)
!  where A_{nmT}(q,G)=<n,0|e^{-i(G+q)x}|m,T>
!deallocate(mexp)
return
end