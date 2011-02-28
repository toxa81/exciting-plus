subroutine sic_genfvprj
use modmain
use mod_sic
implicit none
integer ik,ikloc,ispn,ld,ierr,j1,j2,i
integer ias,ig,ist,ir,j
complex(8) zt1
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: ovlp(:,:)
complex(8), allocatable :: ovlp1(:,:,:)
complex(8), allocatable :: wb(:,:,:)
complex(8), allocatable :: wvb(:,:,:)
complex(8), allocatable :: expikr(:)
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), external :: zfinp_
!
sic_wb=zzero
sic_wvb=zzero
if (.not.tsic_wv) return 
call timer_start(t_sic_genfvprj)

!do ikloc=1,nkptloc
!  do j=1,sic_wantran%nwan
!    n=sic_wantran%iwan(j)
!    do ispn=1,nspinor
!      do ist=1,nstfv
!        istsv=ist+(ispn-1)*nstfv
!        do i=1,nstsv
!          sic_wb_tmp(j,ist,ispn,ikloc)=sic_wb_tmp(j,ist,ispn,ikloc)+&
!            dconjg(wann_c(n,i,ikloc)*evecsvloc(istsv,i,ikloc))
!        enddo !j
!      enddo !ispn
!    enddo !i
!  enddo !n
!enddo !ikloc  

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wffvmt(lmmaxvr,nufrmax,natmtot,nstfv))
allocate(ovlp(sic_wantran%nwan,sic_wantran%nwan))
allocate(ovlp1(sic_wantran%nwan,sic_wantran%nwan,nkpt))
allocate(wb(sic_wantran%nwan,nstfv,nspinor))
allocate(wvb(sic_wantran%nwan,nstfv,nspinor))

allocate(expikr(ngrtot))
allocate(wfmt_(lmmaxvr,nrmtmax))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))

! recompute muffin-tin integrals every n iterations
if (mod(iscl,5).eq.1) call sic_genmti

ovlp1=zzero
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),sfacgk(1,1,1,ikloc),&
    apwalm)

! == code #1 == Bloch sums of W_n and (W*V)_n
  do ir=1,ngrtot
    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir)))
  enddo
  do ist=1,nstfv
    wfmt=zzero
    wfir=zzero
! generate first-variational wave function
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
        evecfvloc(1,ist,1,ikloc),lmmaxvr,wfmt_)
! convert to spherical coordinates
!      call zgemm('N','N',lmmaxvr,nrmt(ias2is(ias)),lmmaxvr,zone,zbshtvr,lmmaxvr, &
!       wfmt_,lmmaxvr,zzero,wfmt(1,1,ias),lmmaxvr)
      wfmt(:,:,ias)=wfmt_(:,:)
    enddo
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig(ig,1,ikloc)))=evecfvloc(ig,ist,1,ikloc)/sqrt(omega)
    enddo
    call zfftifc(3,ngrid,1,wfir)
    do ir=1,ngrtot
      wfir(ir)=wfir(ir)*expikr(ir)
    enddo 
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
        !sic_wb(j,ist,ispn,ikloc)=zfinp_(s_wankmt(1,1,1,ispn,j,ikloc),wfmt,&
        !  s_wankir(1,ispn,j,ikloc),wfir)
        !sic_wvb(j,ist,ispn,ikloc)=zfinp_(s_wvkmt(1,1,1,ispn,j,ikloc),wfmt,&
        !  s_wvkir(1,ispn,j,ikloc),wfir)
      enddo
    enddo
  enddo

  call genwffvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfvloc(1,1,1,ikloc),&
    apwalm,wffvmt)

! == code #2 == expansion of fv function in the big sphere
!  do ist=1,nstfv
!    do ir=1,s_nr
!      do itp=1,s_ntp
!        vrc(:)=s_spx(:,itp)*s_r(ir)+wanpos <---!!!
!        call s_get_wffvval(ikloc,vrc,wffvmt(1,1,1,ist),&
!          evecfvloc(1,ist,1,ikloc),wffvtp(itp,ir))
!      enddo
!    enddo
!    call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,wffvtp,&
!      s_ntp,zzero,wffvlm,lmmaxwan)
!    do j=1,sic_wantran%nwan
!      n=sic_wantran%iwan(j)
!      do ispn=1,nspinor
!        zt1=zzero
!        zt2=zzero
!        do ir=1,s_nr
!          zt1=zt1+zdotc(lmmaxwan,s_wanlm(1,ir,ispn,j),1,wffvlm(1,ir),1)*s_rw(ir)
!          zt2=zt2+zdotc(lmmaxwan,s_wvlm(1,ir,ispn,j),1,wffvlm(1,ir),1)*s_rw(ir)
!        enddo
!        sic_wb(j,ist,ispn,ikloc)=zt1
!        sic_wvb(j,ist,ispn,ikloc)=zt2
!      enddo !ispn
!    enddo !j
!  enddo !istfv

!! == code #3 == optimized version of code #2
!  ld=lmmaxvr*nufrmax*natmtot
!  do ispn=1,nspinor
!! muffin-tin part of <W_n|\phi>
!    call zgemm('T','N',sic_wantran%nwan,nstfv,ld,zone,&
!      sic_wuy(1,1,1,1,ispn,ikloc),ld,wffvmt,ld,zzero,wb(1,1,ispn),&
!      sic_wantran%nwan)
!      write(*,*)"mt:",wb(1,1,1)
!! interstitial part of <W_n|\phi>
!    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
!      sic_wgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),&
!      nmatmax,zone,wb(1,1,ispn),sic_wantran%nwan)
!      write(*,*)"tot:",wb(1,1,1)     
!! muffin tin part of <(W*V)_n|\phi>
!    call zgemm('T','N',sic_wantran%nwan,nstfv,ld,zone,&
!      sic_wvuy(1,1,1,1,ispn,ikloc),ld,wffvmt,ld,zzero,wvb(1,1,ispn),&
!      sic_wantran%nwan)
!! interstitial part of <(W*V)_n|\phi>
!    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
!      sic_wvgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),&
!      nmatmax,zone,wvb(1,1,ispn),sic_wantran%nwan)
!  enddo !ispn
!  sic_wb(:,:,:,ikloc)=wb(:,:,:)
!  sic_wvb(:,:,:,ikloc)=wvb(:,:,:)

  wb=sic_wb(:,:,:,ikloc)
  wvb=sic_wvb(:,:,:,ikloc)
  write(*,*)"tot:",wb(1,1,1)

! compute overlap matrix
  ovlp=zzero
  do j1=1,sic_wantran%nwan
    do j2=1,sic_wantran%nwan
      do ispn=1,nspinor
        do i=1,nstfv
          ovlp(j1,j2)=ovlp(j1,j2)+wb(j1,i,ispn)*dconjg(wb(j2,i,ispn))
        enddo
      enddo
    enddo
  enddo
  ovlp1(:,:,ik)=ovlp(:,:)
!! compute O^{-1/2}
!  call isqrtzhe(sic_wantran%nwan,ovlp,ierr)
!  if (ierr.ne.0) then
!    write(*,'("Warning(sic_genfvprj): overlap matrix is degenerate")')
!    write(*,'("  iteration : ",I4)')iscl
!    do j1=1,sic_wantran%nwan
!      write(*,'(255F12.6)')(abs(ovlp1(j1,j2)),j2=1,sic_wantran%nwan)
!    enddo
!    sic_wb(:,:,:,ikloc)=wb(:,:,:)
!    sic_wvb(:,:,:,ikloc)=wvb(:,:,:)
!  else
!    do j1=1,sic_wantran%nwan
!      do j2=1,sic_wantran%nwan
!        sic_wb(j1,:,:,ikloc)=sic_wb(j1,:,:,ikloc)+ovlp(j2,j1)*wb(j2,:,:)
!        sic_wvb(j1,:,:,ikloc)=sic_wvb(j1,:,:,ikloc)+ovlp(j2,j1)*wvb(j2,:,:)
!      enddo
!    enddo
!  endif
enddo !ikloc
call mpi_grid_reduce(ovlp1(1,1,1),sic_wantran%nwan*sic_wantran%nwan*nkpt,dims=(/dim_k/))
if (mpi_grid_root()) then
  open(210,file="SIC_GENFVPRJ.OUT",form="formatted",status="replace")
  do ik=1,nkpt
    write(210,'(" ik : ",I4)')ik
    do j1=1,sic_wantran%nwan
      write(210,'(255F12.6)')(abs(ovlp1(j1,j2,ik)),j2=1,sic_wantran%nwan)
    enddo
    write(210,*)
  enddo
  close(210)
endif
deallocate(apwalm,wffvmt)
deallocate(ovlp,wb,wvb,ovlp1)
deallocate(expikr,wfmt_,wfmt,wfir)
call timer_stop(t_sic_genfvprj)
return
end
