subroutine sic_genfvprj
use modmain
use mod_sic
implicit none
!integer ngk1,h
!integer, allocatable :: igkig1(:)
!real(8), allocatable :: gkc1(:)
!real(8), allocatable :: tpgkc1(:,:)
!real(8), allocatable :: vgkl1(:,:)
!real(8), allocatable :: vgkc1(:,:)
!complex(8), allocatable :: sfacgk1(:,:)
!complex(8), allocatable :: evecfv1(:,:,:)
!complex(8), allocatable :: wfir(:)
!complex(8), allocatable :: wfmt_(:,:)
!complex(8), allocatable :: wfir_(:)
!complex(8), allocatable :: a(:,:,:)
!complex(8), allocatable :: b(:,:,:)
!complex(8), allocatable :: expikr(:)
integer ik,ikloc,ir,n,ispn,istfv,istsv,j,i
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: wu(:,:,:,:)
complex(8), allocatable :: wuk(:,:,:)
complex(8), allocatable :: wvu(:,:,:,:)
complex(8), allocatable :: wvuk(:,:,:)
integer l,lm,io,it,ias
real(8) d1
complex(8), external :: zdotu
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
        do istfv=1,nstfv
          istsv=istfv+(ispn-1)*nstfv
          do i=1,nstsv
! TODO: zgemm?
            sic_wb(j,istfv,ispn,ikloc)=sic_wb(j,istfv,ispn,ikloc)+&
              dconjg(wann_c(n,i,ikloc)*evecsvloc(istsv,i,ikloc))
          enddo !j
        enddo !ispn
      enddo !i
    enddo !n
  enddo !ikloc  
  call timer_stop(t_sic_genfvprj)
  return
endif
allocate(wffvmt(lmmaxvr,nufrmax,natmtot,nstfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
! experimental code
allocate(wu(lmmaxvr,nufrmax,natmtot,sic_orbitals%ntr))
allocate(wuk(lmmaxvr,nufrmax,natmtot))
allocate(wvu(lmmaxvr,nufrmax,natmtot,sic_orbitals%ntr))
allocate(wvuk(lmmaxvr,nufrmax,natmtot))
! muffin-tin contribution to <W_n|\phi> and <W_n|V_n|\phi>
do j=1,sic_wantran%nwan
  do ispn=1,nspinor
    wu=zzero
    wvu=zzero
! precompute <W_n|u_{l}> and <W_n|V_n|u_{l}>
    do it=1,sic_orbitals%ntr
      do i=1,nmtloc
        ias=(mtoffs+i-1)/nrmtmax+1
        ir=mod(mtoffs+i-1,nrmtmax)+1
        do lm=1,lmmaxvr
          l=lm2l(lm)
          do io=1,nufr(l,ias2is(ias))
            wu(lm,io,ias,it)=wu(lm,io,ias,it)+rmtwt(i)*&
              dconjg(sic_orbitals%wanmt(lm,i,it,ispn,j))*ufr(ir,l,io,ias2ic(ias))
            wvu(lm,io,ias,it)=wvu(lm,io,ias,it)+rmtwt(i)*&
              dconjg(sic_orbitals%wvmt(lm,i,it,ispn,j))*ufr(ir,l,io,ias2ic(ias))
          enddo !io
        enddo !lm
      enddo !i
    enddo !it
    call mpi_grid_reduce(wu(1,1,1,1),&
      lmmaxvr*nufrmax*natmtot*sic_orbitals%ntr,all=.true.)
    call mpi_grid_reduce(wvu(1,1,1,1),&
      lmmaxvr*nufrmax*natmtot*sic_orbitals%ntr,all=.true.)
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      wuk=zzero
      wvuk=zzero
      do it=1,sic_orbitals%ntr
        d1=dot_product(vkc(:,ik),sic_orbitals%vtc(:,it))
        wuk(:,:,:)=wuk(:,:,:)+wu(:,:,:,it)*exp(zi*d1)
        wvuk(:,:,:)=wvuk(:,:,:)+wvu(:,:,:,it)*exp(zi*d1)
      enddo
      call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
        sfacgk(1,1,1,ikloc),apwalm)
      call genwffvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfvloc(1,1,1,ikloc),&
        apwalm,wffvmt)
      do istfv=1,nstfv
        sic_wb(j,istfv,ispn,ikloc)=sic_wb(j,istfv,ispn,ikloc)+& 
          zdotu(lmmaxvr*nufrmax*natmtot,wuk,1,wffvmt(1,1,1,istfv),1)
        sic_wvb(j,istfv,ispn,ikloc)=sic_wvb(j,istfv,ispn,ikloc)+& 
          zdotu(lmmaxvr*nufrmax*natmtot,wvuk,1,wffvmt(1,1,1,istfv),1)          
!        do ias=1,natmtot
!          do lm=1,lmmaxvr
!            do io=1,nufr(lm2l(lm),ias2is(ias))
!              sic_wb(j,istfv,ispn,ikloc)=sic_wb(j,istfv,ispn,ikloc)+&
!                wuk(lm,io,ias)*wffvmt(lm,io,ias,istfv)
!              sic_wvb(j,istfv,ispn,ikloc)=sic_wvb(j,istfv,ispn,ikloc)+&
!                wvuk(lm,io,ias)*wffvmt(lm,io,ias,istfv)
!            enddo !io
!          enddo !lm
!        enddo !ias
      enddo !istfv
    enddo !ikloc
  enddo !ispn
enddo !j
deallocate(wu,wuk,wvu,wvuk)
! interstitial contribution to <W_n|\phi> and <W_n|V_n|\phi>
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do ispn=1,nspinor
    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
      sic_wgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),nmatmax,&
      zone,sic_wb(1,1,ispn,ikloc),sic_wantran%nwan)
    call zgemm('T','N',sic_wantran%nwan,nstfv,ngk(1,ik),zone,&
      sic_wvgk(1,1,ispn,ikloc),ngkmax,evecfvloc(1,1,1,ikloc),nmatmax,&
      zone,sic_wvb(1,1,ispn,ikloc),sic_wantran%nwan)
  enddo
enddo
! old working code
!allocate(igkig1(ngkmax))
!allocate(gkc1(ngkmax))
!allocate(tpgkc1(2,ngkmax))
!allocate(vgkl1(3,ngkmax))
!allocate(vgkc1(3,ngkmax))
!allocate(sfacgk1(ngkmax,natmtot))
!allocate(evecfv1(nmatmax,nstfv,nspnfv))
!allocate(wfir(ngrtot))
!allocate(wfmt_(lmmaxvr,nmtloc))
!allocate(wfir_(ngrloc))
!allocate(a(sic_wantran%nwan,nstfv,nspinor))
!allocate(b(sic_wantran%nwan,nstfv,nspinor))
!allocate(expikr(ngrloc))
!do ik=1,nkpt
!  call timer_start(t_sic_genfvprj_wfk)
!  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
!  if (mpi_grid_dim_pos(dim_k).eq.h) evecfv1(:,:,:)=evecfvloc(:,:,:,ikloc)
!  call mpi_grid_bcast(evecfv1(1,1,1),nmatmax*nstfv*nspnfv,dims=(/dim_k,dim2/),&
!    root=(/h,0,0/))
!  call gengpvec(vkl(1,ik),vkc(1,ik),ngk1,igkig1,vgkl1,vgkc1,gkc1,tpgkc1)
!  call gensfacgp(ngk1,vgkc1,ngkmax,sfacgk1)
!  call match(ngk1,gkc1,tpgkc1,sfacgk1,apwalm)
!  call genwffvmt(lmaxvr,lmmaxvr,ngk1,evecfv1,apwalm,wffvmt)
!  do ir=1,ngrloc
!    expikr(ir)=exp(zi*dot_product(vkc(:,ik),vgrc(:,ir+groffs)))/sqrt(omega)
!  enddo
!  call timer_stop(t_sic_genfvprj_wfk)
!  a=zzero
!  b=zzero
!! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!!  variational Bloch wave-functions
!  do istfv=1,nstfv
!    wfir=zzero
!    call timer_start(t_sic_genfvprj_wfmt)
!    call sic_wavefmt(wffvmt,istfv,wfmt_)
!    call timer_stop(t_sic_genfvprj_wfmt)
!    call timer_start(t_sic_genfvprj_wfir)
!    do ig=1,ngk(1,ik)
!      wfir(igfft(igkig1(ig)))=evecfv1(ig,istfv,1)
!    enddo
!    call zfftifc(3,ngrid,1,wfir)
!    call timer_stop(t_sic_genfvprj_wfir)
!    call timer_start(t_sic_genfvprj_dotp)
!    call sic_copy_ir_z(.true.,wfir,wfir_)
!    do ir=1,ngrloc
!      wfir_(ir)=wfir_(ir)*expikr(ir)
!    enddo
!    do j=1,sic_wantran%nwan
!      n=sic_wantran%iwan(j)
!      do ispn=1,nspinor
!        a(j,istfv,ispn)=sic_dot_lb(.false.,vkc(1,ik),sic_orbitals%wanmt(1,1,1,ispn,j),&
!          sic_orbitals%wanir(1,1,ispn,j),sic_orbitals%twanmt(1,1,n),wfmt_,wfir_)
!        b(j,istfv,ispn)=sic_dot_lb(.false.,vkc(1,ik),sic_orbitals%wvmt(1,1,1,ispn,j),&
!          sic_orbitals%wvir(1,1,ispn,j),sic_orbitals%twanmt(1,1,n),wfmt_,wfir_)
!      enddo
!    enddo !j
!    call timer_stop(t_sic_genfvprj_dotp)
!  enddo !istfv
!  call timer_start(t_sic_genfvprj_dotp)
!  call mpi_grid_reduce(a(1,1,1),sic_wantran%nwan*nstfv*nspinor,root=(/h,0,0/))
!  call mpi_grid_reduce(b(1,1,1),sic_wantran%nwan*nstfv*nspinor,root=(/h,0,0/))
!  if (mpi_grid_dim_pos(dim_k).eq.h) then
!    !sic_wb(:,:,:,ikloc)=a(:,:,:)
!    !sic_wvb(:,:,:,ikloc)=b(:,:,:)
!    write(*,*)sum(abs(sic_wb(:,:,:,ikloc)-a(:,:,:)))
!    write(*,*)sum(abs(sic_wvb(:,:,:,ikloc)-b(:,:,:)))   
!  endif
!  call timer_stop(t_sic_genfvprj_dotp)
!enddo !ik
!deallocate(igkig1,gkc1,tpgkc1,vgkl1,vgkc1,sfacgk1,evecfv1)
!deallocate(wfir,wfmt_,wfir_,a,b,expikr)
deallocate(wffvmt,apwalm)
call timer_stop(t_sic_genfvprj)
return
end
