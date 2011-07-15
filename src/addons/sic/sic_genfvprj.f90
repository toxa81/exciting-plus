subroutine sic_genfvprj(ikloc,evecfv,apwalm)
use modmain
use mod_sic
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
! local variables
integer ik,ispn,ierr
integer ias,ist,j,j1,is,ilo,ir,l,m,lm,io,i,ig
complex(8) zt2(2),z1,z2
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wb(:,:,:,:)
complex(8), allocatable :: om(:,:)
complex(8), allocatable :: zv1(:)
complex(8), allocatable :: zv2(:)
logical, parameter :: tsic_ort=.false.
!
sic_wb(:,:,:,ikloc)=zzero
sic_wvb(:,:,:,ikloc)=zzero
if (.not.tsic_wv) return 
if (debug_level.ge.4) then
  allocate(wb(3,sic_wantran%nwan,nstfv,nspinor))
  wb=zzero
  allocate(om(sic_wantran%nwan,sic_wantran%nwan))
endif
call timer_start(t_sic_genfvprj)
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

if (tsveqn)  then

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wfmt_,wfmt,wfmt2,ias,j,ispn,zt2)
allocate(wfmt_(lmmaxapw,nrmtmax))
!allocate(wfmt(mt_ntp,nrmtmax,natmtot))
allocate(wfmt(lmmaxapw,nrmtmax,natmtot))
allocate(wfmt2(lmmaxapw,nrmtmax,natmtot))
!$OMP DO
do ist=1,nstfv
  wfmt=zzero
! generate first-variational wave function
  do ias=1,natmtot
    call wavefmt(1,lmaxapw,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
      evecfv(1,ist,1),lmmaxapw,wfmt_)
    wfmt(:,:,ias)=wfmt_(:,:)
! convert to spherical coordinates
    !call zgemm('T','N',mt_ntp,nrmt(ias2is(ias)),lmmaxapw,zone,mt_ylmf,&
    !  lmmaxapw,wfmt_,lmmaxapw,zzero,wfmt(1,1,ias),mt_ntp)
  enddo
  do j=1,sic_wantran%nwan
    do ispn=1,nspinor
      do ias=1,natmtot
        do lm=1,lmmaxapw
          wfmt2(lm,:,ias)=s_wkmt(:,lm,ias,ispn,j,ikloc)
        enddo
      enddo
      sic_wb(j,ist,ispn,ikloc)=s_zfinp(.true.,.true.,lmmaxapw,ngk(1,ik),&
        wfmt2,wfmt,s_wkit(1,ispn,j,ikloc),&
        evecfv(1,ist,1),zt2)
      !sic_wb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
      !  s_wkmt(1,1,1,ispn,j,ikloc),wfmt,s_wkit(1,ispn,j,ikloc),&
      !  evecfv(1,ist,1),zt2)
      if (debug_level.ge.4) then  
        wb(1,j,ist,ispn)=sic_wb(j,ist,ispn,ikloc)
        wb(2:3,j,ist,ispn)=zt2(:)
      endif
      !sic_wvb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
      !  s_wvkmt(1,1,1,ispn,j,ikloc),wfmt,s_wvkit(1,ispn,j,ikloc),&
      !  evecfv(1,ist,1))
    enddo
  enddo
enddo !ist
!$OMP END DO
deallocate(wfmt_,wfmt,wfmt2)
!$OMP END PARALLEL

else

allocate(zv1(ngk(1,ik)))
allocate(zv2(ngk(1,ik)))
do j=1,sic_wantran%nwan
  do ispn=1,nspinor
! interstitial contribution from APW
    do ig=1,ngk(1,ik)
      zv1(ig)=dconjg(s_wkit(ig,ispn,j,ikloc))
      zv2(ig)=dconjg(s_wvkit(ig,ispn,j,ikloc))
    enddo
! muffin-tin contribution from APW
    do ias=1,natmtot
      is=ias2is(ias)
      do lm=1,lmmaxapw
        l=lm2l(lm)
        do io=1,apword(l,is)
          z1=zzero
          z2=zzero
          do ir=1,nrmt(is)
            z1=z1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*&
              apwfr(ir,1,io,l,ias)*mt_rw(ir,is)
            z2=z2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*&
              apwfr(ir,1,io,l,ias)*mt_rw(ir,is)
          enddo
          do ig=1,ngk(1,ik)
            zv1(ig)=zv1(ig)+z1*apwalm(ig,io,lm,ias)
            zv2(ig)=zv2(ig)+z2*apwalm(ig,io,lm,ias)
          enddo !ig
        enddo !io
      enddo !lm
! muffin-tin contribution from l.o.
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        do m=-l,l
          lm=idxlm(l,m)
          i=ngk(1,ik)+idxlo(lm,ilo,ias)
          z1=zzero
          z2=zzero
          do ir=1,nrmt(is)
            z1=z1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*&
              lofr(ir,1,ilo,ias)*mt_rw(ir,is)
            z2=z2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*&
              lofr(ir,1,ilo,ias)*mt_rw(ir,is)
          enddo
          sic_wb(j,i,ispn,ikloc)=z1
          sic_wvb(j,i,ispn,ikloc)=z2
        enddo !m
      enddo !ilo
    enddo !ias
    sic_wb(j,1:ngk(1,ik),ispn,ikloc)=zv1(1:ngk(1,ik))
    sic_wvb(j,1:ngk(1,ik),ispn,ikloc)=zv2(1:ngk(1,ik))
  enddo !ispn
enddo !j
deallocate(zv1,zv2)
endif

call timer_stop(t_sic_genfvprj)
if (debug_level.ge.4) then
  call dbg_open_file
  write(fdbgout,*)
  write(fdbgout,'("[sic_genfvprj] ik : ",I4,"    vkl : ",3F10.6)')ik,vkl(:,ik)
  do ist=1,nstfv
    write(fdbgout,'("  ist : ",I4)')ist
    do j=1,sic_wantran%nwan 
      write(fdbgout,'("    j  : ",I4)')j
      do ispn=1,nspinor
        write(fdbgout,'("     total : ",F12.8," (",2F12.8,")")')abs(wb(1,j,ist,ispn)),&
          dreal(wb(1,j,ist,ispn)),dimag(wb(1,j,ist,ispn))
        write(fdbgout,'("        mt : ",F12.8," (",2F12.8,")")')abs(wb(2,j,ist,ispn)),&
          dreal(wb(2,j,ist,ispn)),dimag(wb(2,j,ist,ispn))
        write(fdbgout,'("        it : ",F12.8," (",2F12.8,")")')abs(wb(3,j,ist,ispn)),&
          dreal(wb(3,j,ist,ispn)),dimag(wb(3,j,ist,ispn))
      enddo
    enddo
  enddo
  om=zzero
  do j=1,sic_wantran%nwan
    do j1=1,sic_wantran%nwan
      do ispn=1,nspinor
        do ist=1,nstfv
          om(j,j1)=om(j,j1)+wb(1,j,ist,ispn)*dconjg(wb(1,j1,ist,ispn))
        enddo
      enddo
    enddo
  enddo
  write(fdbgout,'("  overlap matrix")')
  do j=1,sic_wantran%nwan
    write(fdbgout,'(255F12.6)')(abs(om(j,j1)),j1=1,sic_wantran%nwan)
  enddo
  call dbg_close_file
  deallocate(wb,om)
endif
if (tsic_ort) then
  allocate(om(sic_wantran%nwan,sic_wantran%nwan))
  allocate(wb(2,sic_wantran%nwan,nstfv,nspinor))
  wb(1,:,:,:)=sic_wb(:,:,:,ikloc)
  wb(2,:,:,:)=sic_wvb(:,:,:,ikloc)
  om=zzero
  do j=1,sic_wantran%nwan
    do j1=1,sic_wantran%nwan
      do ispn=1,nspinor
        do ist=1,nstfv
          om(j,j1)=om(j,j1)+wb(1,j,ist,ispn)*dconjg(wb(1,j1,ist,ispn))
        enddo
      enddo
    enddo
  enddo
  call isqrtzhe(sic_wantran%nwan,om,ierr)
  if (ierr.ne.0) then
    write(*,'("Warning(sic_genfvprj): overlap matrix is degenerate")')
  else
    sic_wb(:,:,:,ikloc)=zzero
    sic_wvb(:,:,:,ikloc)=zzero
    do j=1,sic_wantran%nwan
      do j1=1,sic_wantran%nwan
        sic_wb(j,:,:,ikloc)=sic_wb(j,:,:,ikloc)+om(j1,j)*wb(1,j1,:,:)
        sic_wvb(j,:,:,ikloc)=sic_wvb(j,:,:,ikloc)+om(j1,j)*wb(2,j1,:,:)
      enddo
    enddo
  endif
  deallocate(om,wb)
endif
return
end
