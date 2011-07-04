subroutine rhomagk_exact(ikloc,evecfv,evecsv)
use modmain
use mod_seceqn
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ist,j,ilo,ia,io,l,lm,m
integer is,ias,ik,ispn,jst,n,i1
integer ivg3(3),ifg3
integer ig1,ig2,ic,io1,io2,l1,l2,j1,j2,lm1,lm2,lm3,l3
real(8) t1
complex(8) zt1,zt2(nspinor)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)

complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: rhoitk(:,:)
complex(8), allocatable :: gntmp(:,:)
real(8), allocatable :: wo(:)
integer, allocatable :: istocc(:)
complex(8), allocatable :: evecfd_(:,:)
integer nstocc

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
! generate wave functions in muffin-tins
if (tsveqn) then
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
    sfacgk(:,:,1,ikloc),apwalm)
  call genwfsvocc(lmaxapw,lmmaxapw,ngk(1,ik),nstocc,istocc,evecfv,evecsv,&
    apwalm,wfsvmt,wfsvit)
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
allocate(gntmp(lmmaxapw,lmmaxapw))
do lm3=1,lmmaxvr
  gntmp(:,:)=gntyry(lm3,:,:)
  l3=lm2l(lm3)
  do ias=1,natmtot
    ic=ias2ic(ias)
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
              do jst=1,nstocc
                do ispn=1,nspinor
                  do lm1=l1**2+1,(l1+1)**2
                    zt1=zzero
                    do lm2=l2**2+1,(l2+1)**2
                      zt1=zt1+wfsvmt(lm2,io2,ias,ispn,jst)*gntmp(lm2,lm1)
                    enddo
                    zt2(ispn)=zt2(ispn)+dconjg(wfsvmt(lm1,io1,ias,ispn,jst))*zt1*wo(jst)
                  enddo
                enddo
              enddo !jst
              rhomagmt(j1,j2,lm3,ias,:)=rhomagmt(j1,j2,lm3,ias,:)+dreal(zt2(:))
            endif
          enddo
        enddo
      enddo
    enddo
  enddo !lm3
enddo !ias
deallocate(gntmp)
deallocate(wfsvmt)
call timer_stop(t_rho_mag_mt)
call timer_start(t_rho_mag_it)
do ig1=1,ngk(1,ik)
  do ig2=1,ngk(1,ik)
    ivg3(:)=ivg(:,igkig(ig2,1,ikloc))-ivg(:,igkig(ig1,1,ikloc))
    ifg3=igfft(ivgig(ivg3(1),ivg3(2),ivg3(3)))
    do jst=1,nstocc
      do ispn=1,nspinor
        rhomagit(ifg3,ispn)=rhomagit(ifg3,ispn)+&
          dconjg(wfsvit(ig1,ispn,jst))*wfsvit(ig2,ispn,jst)*wo(jst)/omega
      enddo
    enddo
  enddo
enddo
call timer_stop(t_rho_mag_it)
deallocate(wfsvit)
deallocate(wo,istocc)
return
end subroutine

