subroutine sic_genmti
use modmain
use mod_sic
implicit none
integer ikloc,ik,j,n,ispn,ias,ir,itp
real(8) x(3)
complex(8) ufrval(lmmaxvr,nufrmax)
complex(8), allocatable :: wtp(:,:)
complex(8), allocatable :: wvtp(:,:)
complex(8), external :: zdotu
!
allocate(wtp(sic_wantran%nwan,nspinor))
allocate(wvtp(sic_wantran%nwan,nspinor))
sic_wuy=zzero
sic_wvuy=zzero
do ir=1,s_nr
  do itp=1,s_ntp
    x(:)=s_spx(:,itp)*s_r(ir)
    wtp=zzero
    wvtp=zzero
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
        wtp(j,ispn)=zdotu(lmmaxwan,s_wanlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)*&
          s_tpw(itp)*s_rw(ir)
        wvtp(j,ispn)=zdotu(lmmaxwan,s_wvlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)*&
          s_tpw(itp)*s_rw(ir)
      enddo
    enddo
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        call s_get_ufrval(x,vkc(:,ik),ias,ufrval)
        if (ias.ne.-1) then
          do ispn=1,nspinor
            sic_wuy(:,:,ias,j,ispn,ikloc)=sic_wuy(:,:,ias,j,ispn,ikloc)+&
              dconjg(wtp(j,ispn))*ufrval(:,:)
            sic_wvuy(:,:,ias,j,ispn,ikloc)=sic_wvuy(:,:,ias,j,ispn,ikloc)+&
              dconjg(wvtp(j,ispn))*ufrval(:,:)
          enddo !ispn
        endif
      enddo !j
    enddo !ikloc
  enddo !itp
enddo !ir
deallocate(wtp,wvtp)

!allocate(wantp(s_ntp,s_nr,nspinor))
!allocate(wvtp(s_ntp,s_nr,nspinor))
!do ikloc=1,nkptloc
!  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!  do j=1,sic_wantran%nwan
!    n=sic_wantran%iwan(j)
!    do ispn=1,nspinor
!! convert functions to spherical coordinates
!      call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,&
!        s_wanlm(1,1,ispn,j),lmmaxwan,zzero,wantp(1,1,ispn),s_ntp)
!      call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,&
!        s_wvlm(1,1,ispn,j),lmmaxwan,zzero,wvtp(1,1,ispn),s_ntp)
!    enddo
!    do ias=1,natmtot
!      do ir=1,s_nr
!        do itp=1,s_ntp
!          x(:)=s_spx(:,itp)*s_r(ir) !+pos of W_n
!          call s_get_ufrval(ias,x,vkc(:,ik),ufrval)
!          do ispn=1,nspinor
!            sic_wufr(:,:,ias,ispn,j,ikloc)=sic_wufr(:,:,ias,ispn,j,ikloc)+&
!              dconjg(wantp(itp,ir,ispn))*ufrval(:,:)*s_tpw(itp)*s_rw(ir)
!            sic_wvufr(:,:,ias,ispn,j,ikloc)=sic_wvufr(:,:,ias,ispn,j,ikloc)+&
!              dconjg(wvtp(itp,ir,ispn))*ufrval(:,:)*s_tpw(itp)*s_rw(ir)
!          enddo
!        enddo
!      enddo !ir
!    enddo !ias
!  enddo !j    
!enddo !ikloc
!deallocate(wantp,wvtp)
return
end
