subroutine sic_genfvprj
use modmain
use mod_sic
implicit none
!complex(8), allocatable :: expikr(:)
integer ik,ikloc,ir,n,ispn,ist,istsv,j,i,jas,is,ig,ld
real(8) x(3),t1
complex(8) expikt
complex(8), allocatable :: apwalm(:,:,:,:)
!complex(8), allocatable :: wfmt(:,:,:)
!complex(8), allocatable :: wfmt_(:,:)
!complex(8), allocatable :: wfir(:)
!complex(8), allocatable :: wffvtp(:,:)
!complex(8), allocatable :: wffvlm(:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
!complex(8), allocatable :: wantp(:,:,:)
!complex(8), allocatable :: wvtp(:,:,:)

!integer l,lm,io,it,ias,itp,i1,i2,i3
!real(8) d1,vrc(3)
!complex(8) zt1,zt2
!complex(8) ufrval(lmmaxvr,nufrmax)
complex(8), external :: zdotu,zdotc,zfinp_
!real(8), allocatable :: tp(:,:)
!
sic_wb=zzero
sic_wvb=zzero
call timer_start(t_sic_genfvprj)
! this code is mainly for tests 
!  on first SIC iteration Wannier functions are generated from LDA Hamiltonian
!  so we can compute overlap between Wannier states and first-variational
!  states analytically
if (.not.tsic_wv) then
  do ikloc=1,nkptloc
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        do ist=1,nstfv
          istsv=ist+(ispn-1)*nstfv
          do i=1,nstsv
! TODO: zgemm?
            sic_wb(j,ist,ispn,ikloc)=sic_wb(j,ist,ispn,ikloc)+&
              dconjg(wann_c(n,i,ikloc)*evecsvloc(istsv,i,ikloc))
          enddo !j
        enddo !ispn
      enddo !i
    enddo !n
  enddo !ikloc  
  call timer_stop(t_sic_genfvprj)
  return
endif

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
!allocate(expikr(ngrtot))
!allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
!allocate(wfmt_(lmmaxvr,nrmtmax))
!allocate(wfir(ngrtot))
!allocate(tp(2,lmmaxvr))
!call sphcover(lmmaxvr,tp)

!allocate(wffvtp(s_ntp,s_nr))
!allocate(wffvlm(lmmaxwan,s_nr))
allocate(wffvmt(lmmaxvr,nufrmax,natmtot,nstfv))
!allocate(wantp(s_ntp,s_nr,nspinor))
!allocate(wvtp(s_ntp,s_nr,nspinor))


! recompute muffin-tin integrals every n iterations
if (mod(iscl,5).eq.1) call sic_genmti

do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),sfacgk(1,1,1,ikloc),&
    apwalm)

! == code #1 == Bloch sums of W_n and (W*V)_n
!  do ir=1,ngrtot
!    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir)))
!  enddo
!  do ist=1,nstfv
!    wfmt=zzero
!    wfir=zzero
!! generate first-variational wave function
!    do ias=1,natmtot
!      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
!        evecfvloc(1,ist,1,ikloc),lmmaxvr,wfmt_)
!! convert to spherical coordinates
!      call zgemm('N','N',lmmaxvr,nrmt(ias2is(ias)),lmmaxvr,zone,zbshtvr,lmmaxvr, &
!       wfmt_,lmmaxvr,zzero,wfmt(1,1,ias),lmmaxvr)
!    enddo
!    do ig=1,ngk(1,ik)
!      wfir(igfft(igkig(ig,1,ikloc)))=evecfvloc(ig,ist,1,ikloc)/sqrt(omega)
!    enddo
!    call zfftifc(3,ngrid,1,wfir)
!    do ir=1,ngrtot
!      wfir(ir)=wfir(ir)*expikr(ir)
!    enddo 
!    do j=1,sic_wantran%nwan
!      do ispn=1,nspinor
!        sic_wb(j,ist,ispn,ikloc)=zfinp_(s_wankmt(1,1,1,ispn,j,ikloc),wfmt,&
!          s_wankir(1,ispn,j,ikloc),wfir)
!        sic_wvb(j,ist,ispn,ikloc)=zfinp_(s_wvkmt(1,1,1,ispn,j,ikloc),wfmt,&
!          s_wvkir(1,ispn,j,ikloc),wfir)
!      enddo
!    enddo
!  enddo

  call genwffvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfvloc(1,1,1,ikloc),&
    apwalm,wffvmt)

! == code #2 == expansion of fv function in the big sphere
!  do ist=1,nstfv
!    do ir=1,s_nr
!      do itp=1,s_ntp
!        vrc(:)=s_spx(:,itp)*s_r(ir)
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
!        !write(*,*)"diff = ",abs(sic_wb(j,ist,ispn,ikloc)-zt1),abs(sic_wvb(j,ist,ispn,ikloc)-zt2)
!      enddo !ispn
!    enddo !j
!  enddo !istfv

! == code #3 == optimized version of code #2
  ld=lmmaxvr*nufrmax*natmtot
  do ispn=1,nspinor
    call zgemm('T','N',sic_wantran%nwan,nstfv,ld,zone,&
      sic_wuy(1,1,1,1,ispn,ikloc),ld,wffvmt,ld,zzero,&
      sic_wb(1,1,ispn,ikloc),sic_wantran%nwan)
    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
      sic_wgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),nmatmax,&
      zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan)
    call zgemm('T','N',sic_wantran%nwan,nstfv,ld,zone,&
      sic_wvuy(1,1,1,1,ispn,ikloc),ld,wffvmt,ld,zzero,&
      sic_wvb(1,1,ispn,ikloc),sic_wantran%nwan)
    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
      sic_wvgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),nmatmax,&
      zone,sic_wvb(1,1,ispn,ikloc),sic_wantran%nwan)
  enddo
!  do ist=1,nstfv
!    do ispn=1,nspinor
!      do j=1,sic_wantran%nwan
!        zt1=zzero
!        zt2=zzero
!        do ias=1,natmtot
!          do io=1,nufrmax
!            do lm=1,lmmaxvr
!              zt1=zt1+sic_wuy(lm,io,ias,j,ispn,ikloc)*wffvmt(lm,io,ias,ist)
!              zt2=zt2+sic_wvuy(lm,io,ias,j,ispn,ikloc)*wffvmt(lm,io,ias,ist)
!            enddo
!          enddo
!        enddo
!        do ig=1,ngk(1,ik)
!          zt1=zt1+sic_wgk(ig,j,ispn,ikloc)*evecfvloc(ig,ist,1,ikloc)
!          zt2=zt2+sic_wvgk(ig,j,ispn,ikloc)*evecfvloc(ig,ist,1,ikloc)
!        enddo
!        !write(*,*)"diff = ",abs(sic_wb(j,ist,ispn,ikloc)-zt1),abs(sic_wvb(j,ist,ispn,ikloc)-zt2)
!        sic_wb(j,ist,ispn,ikloc)=zt1
!        sic_wvb(j,ist,ispn,ikloc)=zt2
!      enddo !j
!    enddo !ispn  
!  enddo !ist
enddo !ikloc
deallocate(apwalm,wffvmt)
call timer_stop(t_sic_genfvprj)
return
end
