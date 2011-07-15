subroutine sic_test_fvprj(fout)
use modmain
use mod_sic
use mod_nrkp
implicit none
integer, intent(in) :: fout
integer ikloc,ik,ist,ir,itp,itp1,i,j,j1,n,ispn,istsv,lm,io,ig,ias
integer io1,io2
real(8) vrc(3),t1
complex(8) zt1,zt2(2)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: evecfvnr(:,:,:)
complex(8), allocatable :: evecsvnr(:,:)
complex(8), allocatable :: wffvtp(:,:)
complex(8), allocatable :: wantp(:,:)
complex(8), allocatable :: wffvlm(:,:)
complex(8), allocatable :: wfmt(:,:)
complex(8), allocatable :: fvmt(:,:,:)
complex(8), allocatable :: fvir(:)
complex(8), allocatable :: wnkmt(:,:,:,:,:)
complex(8), allocatable :: wnkir(:,:,:)
complex(8), allocatable :: wb(:,:,:,:)      
complex(8), allocatable :: om(:,:)
!
if (debug_level.lt.4) return

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfvnr(nmatmax,nstfv,nspnfv))
allocate(evecsvnr(nstsv,nstsv))
allocate(fvmt(lmmaxvr,nrmtmax,natmtot))
allocate(fvir(ngrtot))
allocate(wnkmt(lmmaxvr,nrmtmax,natmtot,nspinor,sic_wantran%nwan))
allocate(wnkir(ngrtot,nspinor,sic_wantran%nwan))
allocate(wb(3,sic_wantran%nwan,nstfv,nspinor))
allocate(om(sic_wantran%nwan,sic_wantran%nwan))

do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnr)
  call getevecsv(vklnr(1,ik),evecsvnr)
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),&
    sfacgknr(1,1,ikloc),apwalm)
  wnkmt=zzero
  wnkir=zzero
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    do ispn=1,nspinor
      do ias=1,natmtot
        do ir=1,nrmt(ias2is(ias))
          do lm=1,lmmaxvr
            do io=1,nufr(lm2l(lm),ias2is(ias))
              wnkmt(lm,ir,ias,ispn,j)=wnkmt(lm,ir,ias,ispn,j)+&
                wann_unkmt(lm,io,ias,ispn,n,ikloc)*&
                ufr(ir,lm2l(lm),io,ias2ic(ias))
            enddo
          enddo !lm
        enddo !ir
      enddo !ias
      do ig=1,ngknr(ikloc)
        wnkir(igfft(igkignr(ig,ikloc)),ispn,j)=wann_unkit(ig,ispn,n,ikloc)/&
          sqrt(omega)
      enddo
      call zfftifc(3,ngrid,1,wnkir(1,ispn,j))
    enddo !ispn
  enddo !j
  wb=zzero
  do ist=1,nstfv
    fvmt=zzero
    fvir=zzero
! generate first-variational wave function
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngknr(ikloc),apwalm,&
        evecfvnr(1,ist,1),lmmaxvr,fvmt(1,1,ias))
    enddo
    do ig=1,ngknr(ikloc)
      fvir(igfft(igkignr(ig,ikloc)))=evecfvnr(ig,ist,1)/sqrt(omega)
    enddo
    call zfftifc(3,ngrid,1,fvir)
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
        zt2=zzero   
        wb(1,j,ist,ispn)=s_zfinp(.true.,.false.,lmmaxvr,ngrtot,wnkmt(1,1,1,ispn,j),&
          fvmt,wnkir(1,ispn,j),fvir,zt2)
        wb(2:3,j,ist,ispn)=zt2(:) 
      enddo
    enddo
  enddo !ist
  call dbg_open_file
  write(fdbgout,*)
  write(fdbgout,'("[sic_test_fvprj] ik : ",I4)')ik
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
enddo !ikloc
deallocate(apwalm)
deallocate(evecfvnr)
deallocate(evecsvnr)
deallocate(wnkmt,wnkir)
deallocate(wb)
deallocate(fvmt,fvir)
deallocate(om)
return
end

