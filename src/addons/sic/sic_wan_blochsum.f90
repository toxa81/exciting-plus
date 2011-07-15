subroutine sic_wan_blochsum
use modmain
use mod_sic
use mod_nrkp
use mod_madness
implicit none
!
integer j,n,ikloc,ik,jk,ispn,ias,is,ic,lm,l,io,ig
complex(8), allocatable :: zv(:)
!
allocate(zv(ngrtot))
s_wkmt=zzero
s_wkit=zzero
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  call elk_load_wann_unk(n)
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! quick-fix: find index of non-reduced point
    do jk=1,nkptnr
      if (sum(abs(vklnr(:,jk)-vkl(:,ik))).lt.epslat) exit
    enddo
    do ispn=1,nspinor
! muffin-tin part
      do ias=1,natmtot
        is=ias2is(ias)
        ic=ias2ic(ias)
        do lm=1,lmmaxapw
          l=lm2l(lm)
          do io=1,nufr(l,is)
            s_wkmt(:,lm,ias,ispn,j,ikloc)=s_wkmt(:,lm,ias,ispn,j,ikloc)+&
              m_wann_unkmt(lm,io,ias,ispn,jk)*ufr(:,l,io,ic)
          enddo !io
        enddo !lm
      enddo !ias
! interstitial
      zv=zzero
      do ig=1,m_ngknr(jk)
        zv(igfft(m_igkignr(ig,jk)))=m_wann_unkit(ig,ispn,jk)
      enddo
      call zfftifc(3,ngrid,1,zv)
      zv(:)=zv(:)*cfunir(:)
      call zfftifc(3,ngrid,-1,zv)
      do ig=1,m_ngknr(jk)
        s_wkit(ig,ispn,j,ikloc)=zv(igfft(m_igkignr(ig,jk)))
      enddo !ig
    enddo !ispn
  enddo !ikloc
enddo !j
deallocate(zv)
return
end subroutine
