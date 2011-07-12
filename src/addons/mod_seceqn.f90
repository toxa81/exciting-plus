module mod_seceqn

contains

subroutine sethml(ngp,ld,vgpc,igpig,apwalm,h,ispn,jspn)
use modmain
use mod_timer
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: ld
real(8), intent(in) :: vgpc(3,ngkmax)
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: h(ld,*)
integer, optional, intent(in) :: ispn
integer, optional, intent(in) :: jspn
!
logical tbf
complex(8), allocatable :: zv(:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8) zsum,zt1
real(8) t1
integer is,ia,ias,ig,naa,n
integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,io1,io2
integer i,j,ilo1,ilo2
integer iv(3)

call timer_start(t_seceqnfv_setup_h)
call timer_start(t_seceqnfv_setup_h_mt)
tbf=.false.
if (spinpol.and..not.tsveqn) tbf=.true.
if (tbf.and..not.(present(ispn).and.present(jspn))) then
  write(*,*)
  write(*,'("Error(sethml): spin blocks are not specified")')
  call pstop
endif

allocate(zv(ngp))
allocate(zm1(lmmaxapw*apwordmax*natmtot,ngp))
allocate(zm2(lmmaxapw*apwordmax*natmtot,ngp))
naa=0
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do lm1=1,lmmaxapw
    l1=lm2l(lm1)
    do io1=1,apword(l1,is)
      zv=zzero
      do lm2=1,lmmaxapw
        l2=lm2l(lm2)
        do io2=1,apword(l2,is)
          zsum=0.d0
          do l3=0,lmaxvr
            if (mod(l1+l2+l3,2).eq.0) then
              if (.not.tbf) then
                do lm3=l3**2+1,(l3+1)**2
                  zsum=zsum+gntyry(lm3,lm2,lm1)*haa(lm3,io2,l2,io1,l1,ias)
                enddo !m3
              else
                if (ispn.eq.1.and.jspn.eq.1) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm3,lm2,lm1)*(haa(lm3,io2,l2,io1,l1,ias)+&
                      baa(lm3,io2,l2,io1,l1,ias,1))
                  enddo !m3
                endif
                if (ispn.eq.2.and.jspn.eq.2) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm3,lm2,lm1)*(haa(lm3,io2,l2,io1,l1,ias)-&
                      baa(lm3,io2,l2,io1,l1,ias,1))
                  enddo !m3
                endif
              endif !.not.tbf
            endif !mod(l1+l2+l3,2).eq.0
          enddo !l3
          if (abs(zsum).gt.1.d-14) then
            call zaxpy(ngp,dconjg(zsum),apwalm(1,io2,lm2,ias),1,zv,1)
          endif
        enddo !io2
      enddo !lm2
! kinetic surface contribution
      do io2=1,apword(l1,is)
        zt1=zone*(0.5d0*rmt(is)**2)*apwfr(nrmt(is),1,io1,l1,ias)*&
          apwdfr(io2,l1,ias)
        call zaxpy(ngp,zt1,apwalm(1,io2,lm1,ias),1,zv,1)
      enddo !io2
      naa=naa+1
      zm1(naa,:)=zv(:)
      zm2(naa,:)=apwalm(1:ngp,io1,lm1,ias)
    enddo !io1
  enddo !lm1
enddo !ias
!----------------------!
!     APW-APW term     !
!----------------------!
n=lmmaxapw*apwordmax*natmtot
call zgemm('C','N',ngp,ngp,naa,zone,zm1,n,zm2,n,zone,h,ld)
deallocate(zm1,zm2,zv)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
!---------------------!
!     APW-lo term     !
!---------------------!  
  do ilo2=1,nlorb(is)
    l2=lorbl(ilo2,is)
    do m2=-l2,l2
      lm2=idxlm(l2,m2)
      i=ngp+idxlo(lm2,ilo2,ias)
      do lm1=1,lmmaxapw
        l1=lm2l(lm1)
        do io1=1,apword(l1,is)
          zsum=0.d0
          do l3=0,lmaxvr
            if (mod(l1+l2+l3,2).eq.0) then
              if (.not.tbf) then
                do lm3=l3**2+1,(l3+1)**2
                  zsum=zsum+gntyry(lm3,lm1,lm2)*hloa(lm3,ilo2,io1,l1,ias)
                end do !m3
              else
                if (ispn.eq.1.and.jspn.eq.1) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm3,lm1,lm2)*(hloa(lm3,ilo2,io1,l1,ias)+&
                      bloa(lm3,ilo2,io1,l1,ias,1))
                  end do !m3
                end if
                if (ispn.eq.2.and.jspn.eq.2) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm3,lm1,lm2)*(hloa(lm3,ilo2,io1,l1,ias)-&
                      bloa(lm3,ilo2,io1,l1,ias,1))
                  end do !m3
                end if
              end if !.not.tbf
            end if
          end do !l3
          if (abs(zsum).gt.1.d-14) then
            do ig=1,ngp
              h(ig,i)=h(ig,i)+dconjg(apwalm(ig,io1,lm1,ias))*zsum
            enddo
          endif
        end do !io1
      end do !lm1
    end do !m2
  end do !ilo2
!--------------------!
!     lo-lo term     !
!--------------------!
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      i=ngp+idxlo(lm1,ilo1,ias)
      do ilo2=1,nlorb(is)
        l2=lorbl(ilo2,is)
        do m2=-l2,l2
          lm2=idxlm(l2,m2)
          j=ngp+idxlo(lm2,ilo2,ias)
          if (i.le.j) then
            zsum=0.d0
            do l3=0,lmaxvr
              if (mod(l1+l2+l3,2).eq.0) then
                if (.not.tbf) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm3,lm1,lm2)*hlolo(lm3,ilo1,ilo2,ias)
                  end do
                else
                  if (ispn.eq.1.and.jspn.eq.1) then
                    do lm3=l3**2+1,(l3+1)**2
                      zsum=zsum+gntyry(lm3,lm1,lm2)*(hlolo(lm3,ilo1,ilo2,ias)+&
                        blolo(lm3,ilo1,ilo2,ias,1))
                    end do
                  endif
                  if (ispn.eq.2.and.jspn.eq.2) then
                    do lm3=l3**2+1,(l3+1)**2
                      zsum=zsum+gntyry(lm3,lm1,lm2)*(hlolo(lm3,ilo1,ilo2,ias)-&
                        blolo(lm3,ilo1,ilo2,ias,1))
                    end do
                  endif
                end if
              end if
            end do
            h(i,j)=h(i,j)+zsum
          end if
        end do
      end do
    end do
  end do
enddo
call timer_stop(t_seceqnfv_setup_h_mt)
!---------------------!
!     interstitial    !
!---------------------!
call timer_start(t_seceqnfv_setup_h_it)
do j=1,ngp
  do i=1,j
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    t1=0.5d0*dot_product(vgpc(:,i),vgpc(:,j))
    if (.not.tbf) then
      h(i,j)=h(i,j)+veffig(ig)+t1*cfunig(ig)
    else
      if (ispn.eq.1.and.jspn.eq.1) then
        h(i,j)=h(i,j)+veffig(ig)+beffig(ig,1)+t1*cfunig(ig)
      endif
      if (ispn.eq.2.and.jspn.eq.2) then
        h(i,j)=h(i,j)+veffig(ig)-beffig(ig,1)+t1*cfunig(ig)
      endif
    endif
  end do
end do
call timer_stop(t_seceqnfv_setup_h_it)
call timer_stop(t_seceqnfv_setup_h)
return
end subroutine


subroutine setovl(ngp,nmatp,igpig,apwalm,o)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: nmatp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)
!
complex(8), allocatable :: zm1(:,:)
integer is,ia,ias,ig,naa
integer l1,m1,lm1,io1
integer i,j,ilo1,ilo2
integer iv(3)
!
call timer_start(t_seceqnfv_setup_o)
call timer_start(t_seceqnfv_setup_o_mt)
allocate(zm1(apwordmax*lmmaxapw*natmtot,ngkmax))
zm1=zzero
naa=0
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do l1=0,lmaxapw
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      do io1=1,apword(l1,is)
         naa=naa+1
         zm1(naa,:)=apwalm(:,io1,lm1,ias)
      enddo
    enddo
  enddo
enddo
!----------------------!
!     APW-APW term     !
!----------------------!
call zgemm('C','N',ngp,ngp,naa,zone,zm1,apwordmax*lmmaxapw*natmtot,&
  zm1,apwordmax*lmmaxapw*natmtot,zone,o(1,1),nmatp)
deallocate(zm1)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
!---------------------!
!     APW-lo term     !
!---------------------!  
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      i=ngp+idxlo(lm1,ilo1,ias)
      do io1=1,apword(l1,is)
        do ig=1,ngp
          o(ig,i)=o(ig,i)+dconjg(apwalm(ig,io1,lm1,ias))*oalo(io1,ilo1,ias)
        end do
      end do
    end do
  end do
!--------------------!
!     lo-lo term     !
!--------------------!
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do ilo2=1,nlorb(is)
      if (lorbl(ilo2,is).eq.l1) then
        do m1=-l1,l1
          lm1=idxlm(l1,m1)
          i=ngp+idxlo(lm1,ilo1,ias)
          j=ngp+idxlo(lm1,ilo2,ias)
          if (i.le.j) then
            o(i,j)=o(i,j)+dcmplx(ololo(ilo1,ilo2,ias),0.d0)
          end if
        end do
      end if
    end do
  end do
enddo
call timer_stop(t_seceqnfv_setup_o_mt)
call timer_start(t_seceqnfv_setup_o_it)
!---------------------!
!     interstitial    !
!---------------------!
do j=1,ngp
  do i=1,j
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    o(i,j)=o(i,j)+cfunig(ig)
  end do
end do
call timer_stop(t_seceqnfv_setup_o_it)
call timer_stop(t_seceqnfv_setup_o)
return
end subroutine

! transformation from second-variational to full diagonalization 
!  eigen-vector representaion
subroutine evecsvfd(evecfv,evecsv,evecfd,nbnd,ibnd)
use modmain
implicit none
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstfv,nstsv)
complex(8), intent(out) :: evecfd(nspinor*nmatmax,*)
integer, optional, intent(in) :: nbnd
integer, optional, intent(in) :: ibnd(*) 
!
integer ispn,i,j
logical tbnd
complex(8), allocatable :: evecsv_(:,:)
!
tbnd=.false.
if (present(nbnd)) tbnd=.true.
if (tbnd) then
  allocate(evecsv_(nstsv,nbnd))
  do i=1,nbnd
    evecsv_(:,i)=evecsv(:,ibnd(i))
  enddo
endif
do ispn=1,nspinor
  i=(ispn-1)*nstfv+1
  j=(ispn-1)*nmatmax+1
  if (tbnd) then
    call zgemm('N','N',nmatmax,nbnd,nstfv,zone,evecfv,nmatmax,&
      evecsv_(i,1),nstsv,zzero,evecfd(j,1),nspinor*nmatmax)
  else
    call zgemm('N','N',nmatmax,nstsv,nstfv,zone,evecfv,nmatmax,&
      evecsv(i,1),nstsv,zzero,evecfd(j,1),nspinor*nmatmax)
  endif
enddo
if (tbnd) deallocate(evecsv_)
return
end subroutine


subroutine genwfsvc(lmax,lmmax,ngp,nwf,apwalm,evecfd,wfsvmt,wfsvit)
use modmain
implicit none
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
integer, intent(in) :: nwf
complex(8), intent(in) :: apwalm(ngkmax,lmmaxapw,apwordmax,natmtot)  
complex(8), intent(in) :: evecfd(nspinor*nmatmax,nwf)
complex(8), intent(out) :: wfsvmt(lmmax,nufrmax,natmtot,nspinor,nwf)
complex(8), optional, intent(out) :: wfsvit(ngkmax,nspinor,nwf)
!
integer ispn,j,ias,is,l,m,io,ilo,lm,i1,n
integer ordl(0:lmax)  
complex(8), allocatable :: wfmt(:,:,:,:)
!
if (nwf.eq.0) return
wfsvmt=zzero
if (present(wfsvit)) wfsvit=zzero
allocate(wfmt(lmmaxapw,apwordmax,natmtot,nwf))
n=lmmaxapw*apwordmax*natmtot
do ispn=1,nspinor
  j=(ispn-1)*nmatmax
  call zgemm('T','N',n,nwf,ngp,zone,apwalm,ngkmax,evecfd(j+1,1),&
    nspinor*nmatmax,zzero,wfmt,n)
  do ias=1,natmtot
    ordl=0
    is=ias2is(ias)
    do l=0,lmax
      do io=1,apword(l,is)
        do m=-l,l
          lm=idxlm(l,m)
          wfsvmt(lm,io,ias,ispn,:)=wfmt(lm,io,ias,:)
        enddo !m
      enddo !io
      ordl(l)=apword(l,is)
    enddo !l
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      if (l.le.lmax) then
        ordl(l)=ordl(l)+1
        do m=-l,l
          lm=idxlm(l,m)
          i1=ngp+idxlo(lm,ilo,ias)+j
          wfsvmt(lm,ordl(l),ias,ispn,:)=evecfd(i1,:)
        enddo !m
      endif
    enddo !ilo    
  enddo
  if (present(wfsvit)) then
    wfsvit(1:ngp,ispn,:)=evecfd(j+1:j+ngp,:)
  endif
enddo !ispn
deallocate(wfmt)
return
end subroutine

subroutine genapwalm(ngp,gpc,tpgpc,sfacgp,apwalm)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
real(8), intent(in) :: gpc(ngkmax)
real(8), intent(in) :: tpgpc(2,ngkmax)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot)
complex(8), intent(out) :: apwalm(ngkmax,lmmaxapw,apwordmax,natmtot)
! local variables
integer np,is,ia,ias,omax
integer l,m,lm,io,jo
integer i,ir,igp,info
real(8) fpso,t1
complex(8) zt1,zt2
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: c(:)
real(8), allocatable :: djl(:,:,:)
complex(8), allocatable :: ylmgp(:,:)
complex(8), allocatable :: zd(:,:)
complex(8), allocatable :: zb(:,:)
! external functions
real(8) polynom
external polynom
! polynomial order
np=max(apwordmax+1,4)
allocate(ipiv(np))
allocate(c(np))
allocate(djl(0:lmaxapw,apwordmax,ngp))
allocate(ylmgp(lmmaxapw,ngp))
allocate(zd(apwordmax,apwordmax))
allocate(zb(apwordmax,ngp*(2*lmaxapw+1)))
! compute the spherical harmonics of the G+p-vectors
do igp=1,ngp
  call genylm(lmaxapw,tpgpc(:,igp),ylmgp(:,igp))
end do
fpso=fourpi/sqrt(omega)
! begin loops over atoms and species
do is=1,nspecies
! maximum APW order for this species
  omax=maxval(apword(1:lmaxapw,is))
! evaluate the spherical Bessel function derivatives for all G+p-vectors
  do igp=1,ngp
    t1=gpc(igp)*rmt(is)
    do io=1,omax
      call sbesseldm(io-1,lmaxapw,t1,djl(:,io,igp))
    end do
    t1=1.d0
    do io=2,omax
      t1=t1*gpc(igp)
      djl(:,io,igp)=t1*djl(:,io,igp)
    end do
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! begin loop over l
    do l=0,lmaxapw
      zt1=fpso*zil(l)
! set up matrix of derivatives
      do jo=1,apword(l,is)
        ir=nrmt(is)-np+1
        do io=1,apword(l,is)
          zd(io,jo)=polynom(io-1,np,spr(ir,is),apwfr(ir,1,jo,l,ias),c,rmt(is))
        end do
      end do
! set up target vectors
      i=0
      do igp=1,ngp
        zt2=zt1*sfacgp(igp,ias)
        do m=-l,l
          lm=idxlm(l,m)
          i=i+1
          do io=1,apword(l,is)
            zb(io,i)=djl(l,io,igp)*zt2*conjg(ylmgp(lm,igp))
          end do
        end do
      end do
! solve the general complex linear systems
      call zgesv(apword(l,is),i,zd,apwordmax,ipiv,zb,apwordmax,info)
      if (info.ne.0) then
        write(*,*)
        write(*,'("Error(genapwalm): could not find APW matching coefficients")')
        write(*,'(" for species ",I4)') is
        write(*,'(" and atom ",I4)') ia
        write(*,'(" ZGESV returned INFO = ",I8)') info
        write(*,*)
        call pstop
      end if
      i=0
      do igp=1,ngp
        do m=-l,l
          lm=idxlm(l,m)
          i=i+1
          do io=1,apword(l,is)
            apwalm(igp,lm,io,ias)=zb(io,i)
          end do
        end do
      end do
! end loop over l
    end do
! end loops over atoms and species
  end do
end do
deallocate(ipiv,c,djl,ylmgp,zd,zb)
return
end subroutine

end module
