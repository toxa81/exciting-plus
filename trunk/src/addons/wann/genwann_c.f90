subroutine genwann_c(ik,vpc,e,wfsvmt,wann_c_,ierr_)
use modmain
use mod_wannier
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: e(nstsv)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wann_c_(nwantot,nstsv)
integer, intent(out) :: ierr_
! local variables
complex(8), allocatable :: prjao(:,:)
integer ispn,j,ias,lm,itype,n,ilo
integer itr(3),i,iw,ierr
real(8) tr(3),d1
complex(8) zt1
logical, external :: bndint
real(8), external :: orbwt
integer, external :: hash

if (debug_level.gt.0) then
  if (debug_level.ge.1) then
    call mpi_grid_hash(wfsvmt(1,1,1,1,1),lmmaxvr*nufrmax*natmtot*nspinor*nstsv,&
      dim2,ierr)
    if (ierr.ne.0) call mpi_grid_msg("genwann_c","hash test of wfsvmt failed")
    call mpi_grid_hash(e(1),nstsv,dim2,ierr)
    if (ierr.ne.0) call mpi_grid_msg("genwann_c","hash test of e failed")
  endif
  if (debug_level.ge.5) then
    call dbg_open_file
    write(fdbgout,*)
    write(fdbgout,'("[genwann_c]")')
    write(fdbgout,'("  ik : ",I6)')ik  
    write(fdbgout,'("  hash(e) : ",I16)')hash(e,8*nstsv)
    write(fdbgout,'("  hash(wfsvmt) : ",I16)')&
      hash(wfsvmt,16*lmmaxvr*nufrmax*natmtot*nspinor*nstsv)
    call dbg_close_file
  endif
  if (debug_level.ge.6) then
    call dbg_open_file
    write(fdbgout,*)
    write(fdbgout,'("[genwann_c]")')
    write(fdbgout,'("  ik : ",I6)')ik  
    write(fdbgout,'("  e : ",6G12.6)')e(:)
    call dbg_close_file
  endif
endif

! compute <\psi|g_n>
allocate(prjao(nwantot,nstsv))
prjao=zzero
do n=1,nwantot
  if (.not.wannier_lc) then
    ias=wan_info(1,n)
    lm=wan_info(2,n)
    ispn=wan_info(3,n)
    itype=wan_info(4,n)
    ilo=wannier_prjlo(wan_info(7,n),wan_info(6,n))
    do j=1,nstsv
      if (bndint(j,e(j),wann_eint(1,itype),wann_eint(2,itype))) then
        call genprjao(ias,lm,ispn,j,wfsvmt,prjao(n,j),ilo)
        if (wannier_soft_eint) then
          prjao(n,j)=prjao(n,j)*orbwt(e(j),wannier_soft_eint_e1(itype), &
            wannier_soft_eint_e2(itype),wannier_soft_eint_w1(itype),&
            wannier_soft_eint_w2(itype))
        endif
      endif
    enddo
  else
    do i=1,wanlc_norb(n)
      d1=wanlc_iorb_alpha(i,n)
      iw=wanlc_iorb(1,i,n)
      itr(:)=wanlc_iorb(2:4,i,n)
      tr(:)=avec(:,1)*itr(1)+avec(:,2)*itr(2)+avec(:,3)*itr(3)
      ias=wan_info(1,iw)
      lm=wan_info(2,iw)
      ispn=wan_info(3,iw)
      itype=wan_info(4,iw)
      do j=1,nstsv
        if (bndint(j,e(j),wann_eint(1,itype),wann_eint(2,itype))) then
          call genprjao(ias,lm,ispn,j,wfsvmt,zt1)
          if (wannier_soft_eint) then
            zt1=zt1*orbwt(e(j),wannier_soft_eint_e1(itype),&
              wannier_soft_eint_e2(itype),wannier_soft_eint_w1(itype),&
              wannier_soft_eint_w2(itype))
          endif               
! <psi_k(r)|g(r-T)>=<psi(r+T)|g(r)>=e^{-ikT}<psi(r)|g(r)>
          prjao(n,j)=prjao(n,j)+zt1*d1*exp(-zi*dot_product(vpc,tr))
        endif !bndint
      enddo
    enddo !i
  endif
enddo !n

! remove small contribution
do j=1,nstsv
  d1=0.d0
  do n=1,nwantot
    d1=d1+abs(prjao(n,j))**2
  enddo
  if (d1.lt.wannier_min_prjao) prjao(:,j)=zzero
enddo

if (.false.) then
  write(*,'("Total contribution of projected orbitals : ")')
  do j=1,nstsv
    d1=0.d0
    do n=1,nwantot
      d1=d1+abs(prjao(n,j))**2
    enddo
    write(*,'("  band : ",I4,"  wt : ",F12.6)')j,d1
  enddo
endif

call wann_ort(ik,prjao,ierr_)
wann_c_=prjao
deallocate(prjao)
return
end

subroutine wann_ort(ik,wann_u_mtrx,ierr_)
use modmain
use mod_wannier
implicit none
integer, intent(in) :: ik
complex(8), intent(inout) :: wann_u_mtrx(nwantot,nstsv)
integer, intent(out) :: ierr_
complex(8), allocatable :: s(:,:)
complex(8), allocatable :: sdiag(:)
complex(8), allocatable :: wann_u_mtrx_ort(:,:)
integer m1,m2,j,ierr
integer, external :: hash

if (debug_level.gt.0) then
  if (debug_level.ge.5) then
    call dbg_open_file
    write(fdbgout,*)
    write(fdbgout,'("[wann_ort]")')
    write(fdbgout,'("  ik : ",I6)')ik  
    write(fdbgout,'("  hash(wann_u_mtrx) : ",I16)')hash(wann_u_mtrx,16*nwantot*nstsv)
    call dbg_close_file
  endif
endif

allocate(s(nwantot,nwantot))
allocate(sdiag(nwantot))
allocate(wann_u_mtrx_ort(nwantot,nstsv))

! compute ovelap matrix
s=zzero
do m1=1,nwantot
  do m2=1,nwantot
    do j=1,nstsv
      s(m1,m2)=s(m1,m2)+dconjg(wann_u_mtrx(m1,j))*wann_u_mtrx(m2,j)
    enddo
  enddo
  sdiag(m1)=s(m1,m1)
enddo
! compute S^{-1/2}
call isqrtzhe(nwantot,s,ierr)
ierr_=ierr
! compute Wannier function expansion coefficients
wann_u_mtrx_ort=zzero
if (ierr.eq.0) then
  do m1=1,nwantot
    do m2=1,nwantot
      wann_u_mtrx_ort(m1,:)=wann_u_mtrx_ort(m1,:)+wann_u_mtrx(m2,:)*s(m2,m1)
    enddo
  enddo
  wann_u_mtrx=wann_u_mtrx_ort
endif
deallocate(s,sdiag,wann_u_mtrx_ort)
return
end

real(8) function orbwt(e,e1,e2,w1,w2)
implicit none
real(8), intent(in) :: e
real(8), intent(in) :: e1
real(8), intent(in) :: e2
real(8), intent(in) :: w1
real(8), intent(in) :: w2
orbwt=1.d0/(exp((e-e2)/w2)+1.d0)+1.d0/(exp((-e+e1)/w1)+1.d0)-1.d0
if (orbwt.lt.1d-16) orbwt=0.d0 
return
end
