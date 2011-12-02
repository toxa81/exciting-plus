subroutine rhomagk_exact(ikloc)
use modmain
use mod_seceqn
implicit none
! arguments
integer, intent(in) :: ikloc
! local variables
integer ist,j,ia
integer is,ias,ik,ispn,jst,i1
integer ivg3(3),ifg3,ir
integer ig1,ig2,io1,io2,l1,l2,j1,j2,lm1,lm2,lm3,l3
integer nstocc
real(8) t1
complex(8) zt1,zt2(nspinor)
complex(8) zv1(lmmaxapw,nufrmax)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: wfsvit_(:,:,:)
complex(8), allocatable :: gntmp(:,:)
complex(8), allocatable :: zdens(:,:,:,:,:,:)
real(8), allocatable :: wo(:)
integer, allocatable :: istocc(:)
complex(8), allocatable :: evecfd_(:,:)
complex(8), allocatable ::zfft(:)
real(8), allocatable :: rhotmp(:,:)
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
allocate(wo(nstsv))
allocate(istocc(nstsv))
nstocc=0
do ist=1,nstsv
  t1=wkpt(ik)*occsv(ist,ik)
  if (abs(t1).gt.epsocc) then
    nstocc=nstocc+1
    wo(nstocc)=t1
    istocc(nstocc)=ist
  endif
enddo
allocate(wfsvmt(lmmaxapw,nufrmax,natmtot,nspinor,nstocc))
allocate(wfsvit(ngkmax,nspinor,nstocc))
wfsvmt=zzero
wfsvit=zzero
if (tsveqn) then
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
    sfacgk(:,:,1,ikloc),apwalm)
  call genwfsvocc(lmaxapw,lmmaxapw,ngk(1,ik),nstocc,istocc,&
    evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),apwalm,wfsvmt,wfsvit)
else
  allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
  call genapwalm(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
    sfacgk(:,:,1,ikloc),apwalm)
  allocate(evecfd_(nspinor*nmatmax,nstocc))
  do j=1,nstocc
    evecfd_(:,j)=evecfdloc(:,istocc(j),ikloc)
  enddo
  call genwfsvc(lmaxapw,lmmaxapw,ngk(1,ik),nstocc,apwalm,evecfd_,wfsvmt,wfsvit) 
  deallocate(evecfd_)
endif
deallocate(apwalm)
call timer_start(t_rho_mag_mt)
allocate(zdens(lmmaxapw,lmmaxapw,nufrmax,nufrmax,natmtot,nspinor))
zdens=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(is,l1,io1,l2,io2,jst,ispn,zv1,lm1,lm2) 
do ias=1,natmtot
  is=ias2is(ias)
  do l1=0,lmaxapw
    do io1=1,nufr(l1,is)
      do l2=0,lmaxapw
        do io2=1,nufr(l2,is)
          do jst=1,nstocc
            do ispn=1,nspinor
              zv1(:,:)=wfsvmt(:,:,ias,ispn,jst)
              do lm2=l2**2+1,(l2+1)**2
                do lm1=l1**2+1,(l1+1)**2
                  zdens(lm1,lm2,io1,io2,ias,ispn)=zdens(lm1,lm2,io1,io2,ias,ispn)+&
                    dconjg(zv1(lm1,io1))*zv1(lm2,io2)*wo(jst)
                enddo
              enddo
            enddo
          enddo !jst
        enddo
      enddo
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(gntmp,l1,l2,l3,ias,is,j1,j2,io1,io2,zt2,ispn,zt1)
allocate(gntmp(lmmaxapw,lmmaxapw)) 
!$OMP DO
do lm3=1,lmmaxvr
  gntmp(:,:)=gntyry(lm3,:,:)
  l3=lm2l(lm3)
  do ias=1,natmtot
    is=ias2is(ias)
    j1=0
    do l1=0,lmaxapw
      do io1=1,nufr(l1,is)
        j1=j1+1
        j2=0
        do l2=0,lmaxapw
          do io2=1,nufr(l2,is)
            j2=j2+1
            if (mod(l1+l2+l3,2).eq.0) then
              zt2=zzero
              do ispn=1,nspinor
                zt1=zzero
                do lm2=l2**2+1,(l2+1)**2
                  do lm1=l1**2+1,(l1+1)**2
                    zt1=zt1+zdens(lm1,lm2,io1,io2,ias,ispn)*gntmp(lm1,lm2)
                  enddo
                enddo
                zt2(ispn)=zt1
              enddo
              rhomagmt(j1,j2,lm3,ias,:)=rhomagmt(j1,j2,lm3,ias,:)+dreal(zt2(:))
            endif
          enddo
        enddo
      enddo
    enddo
  enddo !ias
enddo !lm3
!$OMP END DO
deallocate(gntmp)
!$OMP END PARALLEL 
deallocate(wfsvmt)
deallocate(zdens)
call timer_stop(t_rho_mag_mt)
call timer_start(t_rho_mag_it)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfft,rhotmp,jst,ispn,ig1,ir)
allocate(zfft(ngrtot))
allocate(rhotmp(ngrtot,nspinor))
rhotmp=0.d0
!$OMP DO
do jst=1,nstocc
  do ispn=1,nspinor
    zfft=zzero
    do ig1=1,ngk(1,ik) 
      zfft(igfft(igkig(ig1,1,ikloc)))=wfsvit(ig1,ispn,jst)
    enddo
    call zfftifc(3,ngrid,1,zfft) 
    do ir=1,ngrtot
      rhotmp(ir,ispn)=rhotmp(ir,ispn)+(abs(zfft(ir))**2)*wo(jst)/omega
    enddo
  enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP CRITICAL
rhomagit=rhomagit+rhotmp
!$OMP END CRITICAL
deallocate(rhotmp)
!$OMP END PARALLEL

!allocate(wfsvit_(nstocc,ngkmax,nspinor))
!do jst=1,nstocc
!  wfsvit_(jst,:,:)=wfsvit(:,:,jst)
!enddo
!do ig1=1,ngk(1,ik)
!  do ig2=1,ngk(1,ik)
!    ivg3(:)=ivg(:,igkig(ig2,1,ikloc))-ivg(:,igkig(ig1,1,ikloc))
!    ifg3=igfft(ivgig(ivg3(1),ivg3(2),ivg3(3)))
!    do ispn=1,nspinor
!      do jst=1,nstocc
!        rhomagit(ifg3,ispn)=rhomagit(ifg3,ispn)+&
!          dconjg(wfsvit_(jst,ig1,ispn))*wfsvit_(jst,ig2,ispn)*wo(jst)/omega
!      enddo
!    enddo
!  enddo
!enddo
call timer_stop(t_rho_mag_it)
deallocate(wfsvit)
deallocate(wo,istocc)
return
end subroutine

