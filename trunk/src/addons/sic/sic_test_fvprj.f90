subroutine sic_test_fvprj(fout)
use modmain
use mod_sic
use mod_nrkp
implicit none
integer, intent(in) :: fout
integer ikloc,ik,ist,ir,itp,i,j,n,ispn,istsv,lm,io,ig,ias
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
complex(8), allocatable :: zprod(:,:,:)      
!
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
!allocate(wffvmt(lmmaxvr,nufrmax,natmtot,nstfv))
allocate(evecfvnr(nmatmax,nstfv,nspnfv))
allocate(evecsvnr(nstsv,nstsv))
!allocate(wffvtp(s_ntp,s_nr))
!allocate(wantp(s_ntp,s_nr))
!allocate(wffvlm(lmmaxwan,s_nr))
allocate(wfmt(lmmaxvr,nrmtmax))
allocate(fvmt(lmmaxvr,nrmtmax,natmtot))
allocate(fvir(ngrtot))
allocate(wnkmt(lmmaxvr,nrmtmax,natmtot,nspinor,sic_wantran%nwan))
allocate(wnkir(ngrtot,nspinor,sic_wantran%nwan))
allocate(zprod(3,sic_wantran%nwan,nkptnr))
zprod=zzero

do ikloc=1,1 !nkptnrloc
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
        wfmt=zzero
        do ir=1,nrmt(ias2is(ias))
          do lm=1,lmmaxvr
            do io=1,nufr(lm2l(lm),ias2is(ias))
              wfmt(lm,ir)=wfmt(lm,ir)+wann_unkmt(lm,io,ias,ispn,n,ikloc)*&
                ufr(ir,lm2l(lm),io,ias2ic(ias))
            enddo
          enddo !lm
        enddo !ir
        wnkmt(:,:,ias,ispn,j)=wfmt(:,:)
      enddo !ias
      do ig=1,ngknr(ikloc)
        wnkir(igfft(igkignr(ig,ikloc)),ispn,j)=wann_unkit(ig,ispn,n,ikloc)/&
          sqrt(omega)
      enddo
      call zfftifc(3,ngrid,1,wnkir(1,ispn,j))
    enddo !ispn
    zt1=zzero
    zt2=zzero
    do ispn=1,nspinor
      zt1=zt1+s_zfinp(.true.,.false.,lmmaxvr,ngrtot,wnkmt(1,1,1,ispn,j),&
        wnkmt(1,1,1,ispn,j),wnkir(1,ispn,j),wnkir(1,ispn,j),zt2)
    enddo
    zprod(1,j,ik)=zt1
    zprod(2:3,j,ik)=zt2(1:2)
  enddo !j

  do ist=1,nstfv
    !write(220,'("  ist : ",I4)')ist
    fvmt=zzero
    fvir=zzero
! generate first-variational wave function
    do ias=1,natmtot
      call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngknr(ikloc),apwalm,&
        evecfvnr(1,ist,1),lmmaxvr,wfmt)
! convert to spherical coordinates
      !call zgemm('N','N',lmmaxvr,nrmt(ias2is(ias)),lmmaxvr,zone,zbshtvr,&
      !  lmmaxvr,wfmt,lmmaxvr,zzero,fvmt(1,1,ias),lmmaxvr)
      fvmt(:,:,ias)=wfmt(:,:)
    enddo
    do ig=1,ngknr(ikloc)
      fvir(igfft(igkignr(ig,ikloc)))=evecfvnr(ig,ist,1)/sqrt(omega)
    enddo
    call zfftifc(3,ngrid,1,fvir)
    zt2=zzero
    zt1=s_zfinp(.true.,.false.,lmmaxvr,ngrtot,fvmt,fvmt,fvir,fvir,zt2)
    !write(*,*)"ist:",ist,"norm:",zt1,"partial:",dreal(zt2)
    !write(220,'("  norm : ",2G18.10)')dreal(zt1),dimag(zt1) 
    !do j=1,sic_wantran%nwan
    !  n=sic_wantran%iwan(j)
    !  do ispn=1,nspinor
    !    zt1=zfinp_(wnkmt(1,1,1,ispn,j),fvmt,wnkir(1,ispn,j),fvir)
    !    zt2=zzero
    !    istsv=ist+(ispn-1)*nstfv
    !    do i=1,nstsv
    !      zt2=zt2+dconjg(wanncnrloc(n,i,ikloc)*evecsvnr(istsv,i))
    !    enddo !i
    !    !write(220,'("    n : ",I4,"  ispn : ",I4,&
    !    !  &"  numerical : ",2G18.10,"  analytical : ",2G18.10,&
    !    !  &"  diff : ",G18.10)')n,ispn,dreal(zt1),dimag(zt1),&
    !    !  dreal(zt2),dimag(zt2),abs(zt1-zt2)
    !  enddo
    !enddo
  enddo !ist
  !write(220,*)
  !call flushifc(220)
enddo !ikloc
call mpi_grid_reduce(zprod(1,1,1),3*sic_wantran%nwan*nkptnr,dims=(/dim_k/))
if (mpi_grid_root()) then
  open(210,file="sic_blochsum_exact.out",form="formatted",status="replace")
  do ik=1,nkptnr
    write(210,'(" ik : ",I4)')ik
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      write(210,'("  n : ",I4,6X," <W_nk|W_nk> : ",3G18.10)')&
        n,dreal(zprod(1,j,ik)),dreal(zprod(2,j,ik)),dreal(zprod(3,j,ik))
    enddo
    write(210,*)
  enddo !ik
  close(210)
endif



!! == test2: W_{nk} is analytical, integration is analytical ==
!write(220,'("test2")')
!do ikloc=1,1!nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  write(220,'("ik : ",I4,"  vklnr : ",3G18.10)')ik,vklnr(:,ik)
!  call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnr)
!  call getevecsv(vklnr(1,ik),evecsvnr)
!  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),&
!    sfacgknr(1,1,ikloc),apwalm)
!  call genwffvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvnr,apwalm,wffvmt)
!
!  do ist=1,nstfv
!    write(220,'("  ist : ",I4)')ist
!    fvir=zzero
!    do ig=1,ngknr(ikloc)
!      fvir(igfft(igkignr(ig,ikloc)))=evecfvnr(ig,ist,1)
!    enddo
!    call zfftifc(3,ngrid,1,fvir)
!    do ir=1,ngrtot
!      fvir(ir)=fvir(ir)*cfunir(ir)
!    enddo
!    call zfftifc(3,ngrid,-1,fvir)
!    do j=1,sic_wantran%nwan
!      n=sic_wantran%iwan(j)
!      do ispn=1,nspinor
!        zt1=zzero
!        do ias=1,natmtot
!          do lm=1,lmmaxvr
!            do io1=1,nufr(lm2l(lm),ias2is(ias))
!              do io2=1,nufr(lm2l(lm),ias2is(ias))
!                zt1=zt1+dconjg(wann_unkmt(lm,io1,ias,ispn,n,ikloc))*&
!                  wffvmt(lm,io2,ias,ist)*ufrp(lm2l(lm),io1,io2,ias2ic(ias))
!              enddo
!            enddo
!          enddo !lm
!        enddo
!        do ig=1,ngknr(ikloc)
!          zt1=zt1+dconjg(wann_unkit(ig,ispn,n,ikloc))*&
!            fvir(igfft(igkignr(ig,ikloc)))
!        enddo
!        zt2=zzero
!        istsv=ist+(ispn-1)*nstfv
!        do i=1,nstsv
!          zt2=zt2+dconjg(wanncnrloc(n,i,ikloc)*evecsvnr(istsv,i))
!        enddo !i
!        write(220,'("    n : ",I4,"  ispn : ",I4,&
!          &"  numerical : ",2G18.10,"  analytical : ",2G18.10,&
!          &"  diff : ",G18.10)')n,ispn,dreal(zt1),dimag(zt1),&
!          dreal(zt2),dimag(zt2),abs(zt1-zt2)
!      enddo
!    enddo
!  enddo !ist
!  write(220,*)
!  call flushifc(220)
!enddo !ikloc


!t1=0.d0
!
!
!
!
!do ikloc=1,1 !nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  write(220,'("ik : ",I4,"  vklnr : ",3G18.10)')ik,vklnr(:,ik)
!  call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnr)
!  call getevecsv(vklnr(1,ik),evecsvnr)
!  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),&
!    sfacgknr(1,1,ikloc),apwalm)
!  call genwffvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvnr,apwalm,wffvmt)
!  do ist=1,nstfv
!    write(220,'("  ist : ",I4)')ist
!    do ir=1,s_nr
!      do itp=1,s_ntp
!        vrc(:)=s_spx(:,itp)*s_r(ir)
!        call s_get_wffvval(vrc,ngknr(ikloc),vkcnr(1,ik),vgkcnr(1,1,ikloc), &
!          wffvmt(1,1,1,ist),evecfvnr(1,ist,1),wffvtp(itp,ir))
!      enddo
!    enddo
!    do ir=1,3000
!      vrc(:)=(/ir*10.d0/3000.d0,0.d0,0.d0/)
!      call s_get_wffvval(vrc,ngknr(ikloc),vkcnr(1,ik),vgkcnr(1,1,ikloc), &
!          wffvmt(1,1,1,ist),evecfvnr(1,ist,1),zt1)
!      write(200,*)vrc(1),dreal(zt1),dimag(zt1)
!    enddo
!    write(200,*)
!!    call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,wffvtp,&
!!      s_ntp,zzero,wffvlm,lmmaxwan)
!    do j=1,sic_wantran%nwan
!      n=sic_wantran%iwan(j)
!      do ispn=1,nspinor
!        call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,s_wanlm(1,1,ispn,j),&
!          lmmaxwan,zzero,wantp,s_ntp)
!
!        zt1=zzero
!        !do ir=1,s_nr
!        !  zt1=zt1+zdotc(lmmaxwan,s_wanlm(1,ir,ispn,j),1,wffvlm(1,ir),1)*s_rw(ir)
!        !enddo
!        do ir=1,s_nr
!          do itp=1,s_ntp
!            zt1=zt1+dconjg(wantp(itp,ir))*wffvtp(itp,ir)*s_tpw(itp)*s_rw(ir)
!          enddo
!        enddo
!        zt2=zzero
!        istsv=ist+(ispn-1)*nstfv
!        do i=1,nstsv
!          zt2=zt2+dconjg(wanncnrloc(n,i,ikloc)*evecsvnr(istsv,i))
!        enddo !i
!        t1=max(t1,abs(zt1-zt2))
!        write(220,'("    n : ",I4,"  ispn : ",I4,&
!          &"  numerical : ",2G18.10,"  analytical : ",2G18.10,&
!          &"  diff : ",G18.10)')n,ispn,dreal(zt1),dimag(zt1),&
!          dreal(zt2),dimag(zt2),abs(zt1-zt2)
!      enddo !ispn
!    enddo !j
!  enddo !ist
!  write(220,*)
!  call flushifc(220)
!enddo !ikloc
!call mpi_grid_reduce(t1,dims=(/dim_k/),op=op_max)
!if (wproc) then
!  write(fout,*)
!  write(fout,'("Maximum deviation for <W_n|\phi_{jk}> : ",G18.10)')t1
!endif
!close(220)
!call bstop
deallocate(apwalm)
!deallocate(wffvmt)
deallocate(evecfvnr)
deallocate(evecsvnr)
!deallocate(wffvtp)
!deallocate(wffvlm)
!deallocate(wantp)
deallocate(wnkmt,wnkir)
deallocate(wfmt)
deallocate(zprod)
deallocate(fvmt,fvir)
return
end

