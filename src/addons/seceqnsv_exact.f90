subroutine seceqnsv_exact(ikloc,apwalm,evalfv,evecfv,evecsv)
use modmain
use modldapu
use mod_sic
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer is,ias,ik,ic
integer ist,jst,i,j,ispn
integer ir,j1,j2,io1,io2,l1,l2,lm1,lm2,lm3
integer lwork,info,ig2,ig3,ifg1,ivg1(3)
real(8) cb
complex(8) zt1,zt2
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: wfmt(:,:,:,:)
complex(8), allocatable :: wfbfmt(:,:,:,:,:)
complex(8), allocatable :: zfft(:)
complex(8), allocatable :: wfbfit(:,:)
complex(8), allocatable :: zm1(:,:,:,:)
complex(8), allocatable :: zm2(:,:)
! external functions
complex(8) zdotc,zfmtinp
external zdotc,zfmtinp
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_seceqnsv)
call timer_start(t_seceqnsv_setup)
! no calculation of second-variational eigenvectors
if (.not.tevecsv) then  
  do i=1,nstsv
    evalsv(i,ik)=evalfv(i)
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
! generate first-variational wave functions
allocate(wfmt(lmmaxapw,nufrmax,natmtot,nstfv))
call genwffvmt(lmaxapw,lmmaxapw,ngk(1,ik),evecfv,apwalm,wfmt)
allocate(wfbfmt(lmmaxapw,nufrmax,natmtot,nstfv,nspinor))
allocate(zm1(lmmaxapw,lmmaxapw,nufrmax,nufrmax))
call timer_start(t_seceqnsv_setup_mt)
! multiply wave-function with magnetic field
wfbfmt=zzero
do ias=1,natmtot
  is=ias2is(ias)
  ic=ias2ic(ias)
  if (spinpol) then
    zm1=zzero
    j1=0
    do l1=0,lmaxapw
      do io1=1,nufr(l1,is)
        j1=j1+1
        do lm1=l1**2+1,(l1+1)**2
          j2=0
          do l2=0,lmaxapw
            do io2=1,nufr(l2,is)
              j2=j2+1
              do lm2=l2**2+1,(l2+1)**2
                zt2=zzero
                do lm3=1,lmmaxvr
                  zt2=zt2+gntyry(lm3,lm2,lm1)*sv_ubu(lm3,j2,j1,ias,1)
                enddo
                zm1(lm2,lm1,io2,io1)=zt2
              enddo !lm2
            enddo !io2
          enddo !l2
        enddo !lm1
      enddo !io1
    enddo !l1
    do ist=1,nstfv
      do l1=0,lmaxapw
        do io1=1,nufr(l1,is)
          do lm1=l1**2+1,(l1+1)**2
            zt1=zzero
            do l2=0,lmaxapw
              do io2=1,nufr(l2,is)
                do lm2=l2**2+1,(l2+1)**2
                  zt1=zt1+wfmt(lm2,io2,ias,ist)*zm1(lm2,lm1,io2,io1)
                enddo
              enddo
            enddo !l2
            wfbfmt(lm1,io1,ias,ist,1)=zt1
            wfbfmt(lm1,io1,ias,ist,2)=-zt1
          enddo
        enddo
      enddo !l1
    enddo
  endif !spinpol
  if ((ldapu.ne.0).and.(llu(is).ge.0)) then
    l1=llu(is)
    do ispn=1,nspinor
      do ist=1,nstfv
        do io1=1,nufr(l1,is)
          do lm1=l1**2+1,(l1+1)**2
            do io2=1,nufr(l1,is)
              do lm2=l1**2+1,(l1+1)**2
                wfbfmt(lm1,io1,ias,ist,ispn)=wfbfmt(lm1,io1,ias,ist,ispn)+&
                  ufrp(l1,io1,io2,ic)*vmatlu(lm1,lm2,ispn,ispn,ias)*&
                  wfmt(lm2,io2,ias,ist)
              enddo !lm2
            enddo !lm1
          enddo !lm1
        enddo !io1
      enddo !ist
    enddo !ispn
  endif
enddo !ias
deallocate(zm1)
call timer_stop(t_seceqnsv_setup_mt)
call timer_start(t_seceqnsv_setup_it)
allocate(zfft(ngrtot))
allocate(zm2(ngk(1,ik),ngk(1,ik)))
allocate(wfbfit(ngk(1,ik),nstfv))
cb=gfacte/(4.d0*solsc)
do ir=1,ngrtot
  zfft(ir)=zone*(bxcir(ir,1)+cb*bfieldc(3))*cfunir(ir)
enddo
call zfftifc(3,ngrid,-1,zfft)
do ig3=1,ngk(1,ik)
  do ig2=1,ngk(1,ik)
    ivg1(:)=ivg(:,igkig(ig3,1,ikloc))-ivg(:,igkig(ig2,1,ikloc))
    ifg1=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
    zm2(ig2,ig3)=zfft(ifg1)
  enddo
enddo
deallocate(zfft)
call zgemm('T','N',ngk(1,ik),nstfv,ngk(1,ik),zone,zm2,ngk(1,ik),&
  &evecfv,nmatmax,zzero,wfbfit,ngk(1,ik))
deallocate(zm2)
call timer_stop(t_seceqnsv_setup_it)
evecsv=zzero
do ist=1,nstfv
  do jst=1,nstfv
    i=ist
    j=jst
    if (i.le.j) then
      evecsv(i,j)=evecsv(i,j)+zdotc(lmmaxapw*nufrmax*natmtot,&
        &wfmt(1,1,1,ist),1,wfbfmt(1,1,1,jst,1),1)+zdotc(ngk(1,ik),evecfv(1,ist),1,wfbfit(1,jst),1)
    endif
    if (spinpol) then
      i=ist+nstfv
      j=jst+nstfv
      if (i.le.j) then
        evecsv(i,j)=evecsv(i,j)+zdotc(lmmaxapw*nufrmax*natmtot,&
          &wfmt(1,1,1,ist),1,wfbfmt(1,1,1,jst,2),1)-zdotc(ngk(1,ik),evecfv(1,ist),1,wfbfit(1,jst),1)
      endif
    endif
  enddo
enddo
deallocate(wfmt,wfbfmt,wfbfit)
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist)
  end do
end do
call timer_stop(t_seceqnsv_setup)
if (mpi_grid_root((/dim2/))) then
  if (sic) call sic_hunif(ikloc,evecsv)
  call timer_start(t_seceqnsv_diag)
! diagonalise second-variational Hamiltonian
  allocate(rwork(3*nstsv))
  lwork=2*nstsv
  allocate(work(lwork))
  if (ndmag.eq.1) then
! collinear: block diagonalise H
    call zheev('V','U',nstfv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
    if (info.ne.0) goto 20
    i=nstfv+1
    call zheev('V','U',nstfv,evecsv(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
    if (info.ne.0) goto 20
    do i=1,nstfv
      do j=1,nstfv
        evecsv(i,j+nstfv)=0.d0
        evecsv(i+nstfv,j)=0.d0
      end do
    end do
  else
! non-collinear or spin-unpolarised: full diagonalisation
    call zheev('V','U',nstsv,evecsv,nstsv,evalsv(:,ik),work,lwork,rwork,info)
    if (info.ne.0) goto 20
  end if
  deallocate(rwork,work)
  call timer_stop(t_seceqnsv_diag)
endif
call mpi_grid_bcast(evecsv(1,1),nstsv*nstsv,dims=(/dim2/))
call mpi_grid_bcast(evalsv(1,ik),nstsv,dims=(/dim2/))
call timer_stop(t_seceqnsv)
timesv=0.d0
return
20 continue
write(*,*)
write(*,'("Error(seceqnsv1): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
call pstop  
return
end subroutine
