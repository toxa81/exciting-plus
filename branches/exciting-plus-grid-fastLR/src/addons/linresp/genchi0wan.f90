subroutine genchi0wan(igq0,mexp,chi0wan_k,chi0wan,chi0_GqGq_wan_full)
use modmain
implicit none
integer, intent(in) :: igq0
complex(8), intent(in) :: chi0wan_k(nmegqwan,nmegqwan,nkptnrloc)
complex(8), intent(in) :: mexp(nkptnrloc,ntrchi0wan)
complex(8), intent(out) :: chi0wan(nmegqwan,nmegqwan,ntrchi0wan)
complex(8), intent(out) :: chi0_GqGq_wan_full
integer ikloc,ik,it2,i1,i2,n1,n2
real(8) vtrc(3)
complex(8) zt1
!complex(8), allocatable :: mexp(:,:)

call zgemm('N','N',nmegqwan*nmegqwan,ntrchi0wan,nkptnrloc,zone,chi0wan_k,&
  nmegqwan*nmegqwan,mexp,nkptnrloc,zzero,chi0wan,nmegqwan*nmegqwan)


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
call mpi_grid_reduce(chi0wan(1,1,1),nmegqwan*nmegqwan*ntrchi0wan,dims=(/dim_k/))
chi0wan(:,:,:)=chi0wan(:,:,:)/nkptnr/omega
! compute chi0_GqGq using the Wannier functions expansion
!  chi0=\sum_{T,T'}\sum_{n,m,n',m'} A^{*}_{nmT}(q,Gq)*chi0wan(n,m,n',m',T-T')*A_{n'm'T'}(q,Gq)
!  where A_{nmT}(q,G)=<n,0|e^{-i(G+q)x}|m,T>
if (mpi_grid_root(dims=(/dim_k/))) then
  chi0_GqGq_wan_full=zzero
!  do i1=1,ntrmegqwan
!    do n1=1,nmegqwan    
!      do i2=1,ntrmegqwan
!        do n2=1,nmegqwan
!          chi0_GqGq_wan_full=chi0_GqGq_wan_full+&
!            megqwan(n1,i1,igq0)*chi0wan(n1,n2,itridxwan(i1,i2))*&
!            dconjg(megqwan(n2,i2,igq0))
!        enddo !n2
!      enddo !i2
!    enddo !n1
!  enddo !i1
endif
!deallocate(mexp)
return
end