subroutine wann_ene_occ
use modmain
implicit none

real(8) wf_ene(wann_nmax,wann_nspin),wf_occ(wann_nmax,wann_nspin)
integer n,i,ik,ispn,ias,lm1,lm2,l,j,n1,n2,m1,m2
integer, external :: ikglob
real(8) w2
complex(8), allocatable :: wf_ene_mtrx(:,:,:,:,:)
complex(8), allocatable :: wf_occ_mtrx(:,:,:,:,:)
complex(8) z2

wf_ene=0.d0
wf_occ=0.d0

do ispn=1,wann_nspin
  do n=1,nwann(ispn)
    do ik=1,nkptloc(iproc)
      do i=1,nstfv
        w2=dreal(dconjg(wann_c(n,i,ispn,ik))*wann_c(n,i,ispn,ik))
        wf_ene(n,ispn)=wf_ene(n,ispn)+w2*evalsv(i+(ispn-1)*nstfv,ikglob(ik))*wkpt(ikglob(ik))
        wf_occ(n,ispn)=wf_occ(n,ispn)+w2*occsv(i+(ispn-1)*nstfv,ikglob(ik))*wkpt(ikglob(ik))
      enddo
    enddo
  enddo
enddo
call dsync(wf_ene,wann_nmax*wann_nspin,.true.,.false.)
call dsync(wf_occ,wann_nmax*wann_nspin,.true.,.false.)
if (iproc.eq.0) then
  write(60,*)
  do ispn=1,wann_nspin
    write(60,'("Wannier functions, spin ",I1)')ispn
    write(60,'(" wf  energy (Ha)  energy (eV)    occupancy ")')
    write(60,'("-------------------------------------------")')
    do n=1,nwann(ispn)
      write(60,'(1X,I2,1X,F12.6,1X,F12.6,1X,F12.6)')n,wf_ene(n,ispn),wf_ene(n,ispn)*ha2ev,wf_occ(n,ispn)
    enddo
  enddo
  if (nwann(ispn).gt.10) then
    write(60,*)
    write(60,'("E1 : ",F12.6)')sum(wf_ene(1:5,1))/5.d0
    write(60,'("E2 : ",F12.6)')sum(wf_ene(6:10,1))/5.d0
    write(60,'("N1 : ",F12.6)')sum(wf_occ(1:5,1))
    write(60,'("N2 : ",F12.6)')sum(wf_occ(6:10,1))
  endif
endif

allocate(wf_ene_mtrx(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(wf_occ_mtrx(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
wf_occ_mtrx=zzero
wf_ene_mtrx=zzero
do ispn=1,wann_nspin
  do i=1,wann_natom
    ias=wann_iatom(1,i)
    do l=0,lmaxlu
      do m1=-l,l
        do m2=-l,l
          lm1=idxlm(l,m1)
          lm2=idxlm(l,m2)
          n1=iasiwann(ias,lm1,ispn)
          n2=iasiwann(ias,lm2,ispn)
          if (n1.ne.-1.and.n2.ne.-1) then
            do ik=1,nkptloc(iproc)
              do j=1,nstfv
                z2=dconjg(wann_c(n1,j,ispn,ik))*wann_c(n2,j,ispn,ik)
                wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)=&
                  wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)+&
                  z2*occsv(j+(ispn-1)*nstfv,ikglob(ik))*wkpt(ikglob(ik))
                wf_ene_mtrx(lm1,lm2,ispn,ispn,ias)=&
                  wf_ene_mtrx(lm1,lm2,ispn,ispn,ias)+&
                  z2*evalsv(j+(ispn-1)*nstfv,ikglob(ik))*wkpt(ikglob(ik))
              enddo !j
            enddo !ik
          endif
        enddo !m2
      enddo !m1
    enddo !l
  enddo !i
enddo !ispn
call zsync(wf_occ_mtrx,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
call zsync(wf_ene_mtrx,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
! convert from Rlm to Ylm basis
do ias=1,natmtot
  call mtrxbas(lmmaxlu,yrlm_lcs(1,1,ias),wf_occ_mtrx(1,1,1,1,ias))
  call mtrxbas(lmmaxlu,yrlm_lcs(1,1,ias),wf_ene_mtrx(1,1,1,1,ias))
enddo
! symmetrise matrix
call symdmat(lmaxlu,lmmaxlu,wf_occ_mtrx)
call symdmat(lmaxlu,lmmaxlu,wf_ene_mtrx)
! convert back to Rlm
do ias=1,natmtot
  call mtrxbas(lmmaxlu,rylm_lcs(1,1,ias),wf_occ_mtrx(1,1,1,1,ias))
  call mtrxbas(lmmaxlu,rylm_lcs(1,1,ias),wf_ene_mtrx(1,1,1,1,ias))
enddo
! save energies and occupancies
do ispn=1,wann_nspin
  do i=1,wann_natom
    ias=wann_iatom(1,i)
    do l=0,lmaxlu
      do m1=-l,l
        lm1=idxlm(l,m1)
        lm2=idxlm(l,m2)
        n1=iasiwann(ias,lm1,ispn)
        if (n1.ne.-1) then
          wann_ene(n1,ispn)=wf_ene_mtrx(lm1,lm1,ispn,ispn,ias)
          wann_occ(n1,ispn)=wf_occ_mtrx(lm1,lm1,ispn,ispn,ias)
        endif
      enddo !m1
    enddo !l
  enddo !i
enddo !ispn

if (iproc.eq.0) then
  write(60,*)
  write(60,'("On-site matrices in WF basis")')
  do i=1,wann_natom
    ias=wann_iatom(1,i)
    write(60,*)
    write(60,'("ias : ",I4)')ias
    do ispn=1,wann_nspin
      write(60,'("  ispn : ",I1)')ispn
      do l=0,lmaxlu
        if (any(iasiwann(ias,idxlm(l,-l):idxlm(l,l),ispn).gt.0)) then
          write(60,'("     l : ",I1)')l
          write(60,'("  occupancy matrix")')
          write(60,'("    real part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dreal(wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          write(60,'("    imag part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dimag(wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          w2=0.d0
          do lm1=l**2+1,(l+1)**2
            w2=w2+dreal(wf_occ_mtrx(lm1,lm1,ispn,ispn,ias))
          enddo
          write(60,'("    occupancy : ",F12.6)')w2
          write(60,'("  energy matrix")')
          write(60,'("    real part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dreal(wf_ene_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          write(60,'("    imag part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dimag(wf_ene_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          w2=0.d0
          n=0
          do lm1=l**2+1,(l+1)**2
            if (iasiwann(ias,lm1,ispn).gt.0) n=n+1
            w2=w2+dreal(wf_ene_mtrx(lm1,lm1,ispn,ispn,ias))
          enddo
          if (n.ne.0)   write(60,'("    average energy : ",F12.6)')w2/n
        endif
      enddo
    enddo
  enddo
endif


deallocate(wf_ene_mtrx,wf_occ_mtrx)
return
end
  
