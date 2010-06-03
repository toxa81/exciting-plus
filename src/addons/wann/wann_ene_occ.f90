subroutine wann_ene_occ
use modmain
use modldapu
use mod_mpi_grid
implicit none
! local variables
real(8) wf_ene(nwann),wf_occ(nwann),t(2),w2
integer n,i,ik,ispn,ias,lm1,lm2,l,j,n1,n2,m1,m2,ispn1,ispn2,ikloc
!integer, external :: ikglob
complex(8), allocatable :: wf_ene_mtrx(:,:,:,:,:)
complex(8), allocatable :: wf_occ_mtrx(:,:,:,:,:)
complex(8) z2

wf_ene=0.d0
wf_occ=0.d0
do n=1,nwann
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    do i=1,nstsv
      w2=dreal(dconjg(wann_c(n,i,ikloc))*wann_c(n,i,ikloc))
      wf_ene(n)=wf_ene(n)+w2*evalsv(i,ik)*wkpt(ik)
      wf_occ(n)=wf_occ(n)+w2*occsv(i,ik)*wkpt(ik)
    enddo
  enddo
enddo
call mpi_grid_reduce(wf_ene(1),nwann,dims=(/dim_k/),side=.true.)
call mpi_grid_reduce(wf_occ(1),nwann,dims=(/dim_k/),side=.true.)
if (wproc.and.nosym) then
  write(60,*)
  write(60,'(" WF  energy (Ha)  energy (eV)    occupancy ")')
  write(60,'("-------------------------------------------")')
  do n=1,nwann
    write(60,'(1X,I2,1X,F12.6,1X,F12.6,1X,F12.6)')n,wf_ene(n),wf_ene(n)*ha2ev,wf_occ(n)
  enddo
endif

allocate(wf_ene_mtrx(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(wf_occ_mtrx(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
wf_occ_mtrx=zzero
wf_ene_mtrx=zzero
! generate occupancy matrix in WF basis for collinear (!!!) case only
! for noncollinear case the product <W_{n\sigma}|W_{n'\sigma'}> is required
!  (in collinear case it is diagonal in spin index)
do n1=1,nwann
  do n2=1,nwann
    if (iwann(1,n1).eq.iwann(1,n2)) then
      ias=iwann(1,n1)
      lm1=iwann(2,n1)
      lm2=iwann(2,n2)
      ispn1=iwann(3,n1)
      ispn2=iwann(3,n2)
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)        
        do j=1,nstsv
          z2=dconjg(wann_c(n1,j,ikloc))*wann_c(n2,j,ikloc)
          wf_occ_mtrx(lm1,lm2,ispn1,ispn2,ias)=&
            wf_occ_mtrx(lm1,lm2,ispn1,ispn2,ias)+&
            z2*occsv(j,ik)*wkpt(ik)
          wf_ene_mtrx(lm1,lm2,ispn1,ispn2,ias)=&
            wf_ene_mtrx(lm1,lm2,ispn1,ispn2,ias)+&
            z2*evalsv(j,ik)*wkpt(ik)
        enddo !j
      enddo !ikloc
    endif
  enddo
enddo 
call mpi_grid_reduce(wf_occ_mtrx(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),side=.true.,all=.true.)
call mpi_grid_reduce(wf_ene_mtrx(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),side=.true.,all=.true.)
! convert from Rlm to Ylm basis
do ias=1,natmtot
  call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),wf_occ_mtrx(1,1,1,1,ias))
  call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),wf_ene_mtrx(1,1,1,1,ias))
enddo
! symmetrise matrix
call symdmat(lmaxlu,lmmaxlu,wf_occ_mtrx)
call symdmat(lmaxlu,lmmaxlu,wf_ene_mtrx)
! convert back to Rlm
do ias=1,natmtot
  call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),wf_occ_mtrx(1,1,1,1,ias))
  call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),wf_ene_mtrx(1,1,1,1,ias))
  call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),wf_v_mtrx(1,1,1,1,ias))
enddo
if (wproc) then
  write(60,*)
  write(60,'("On-site matrices in WF basis")')
  do i=1,wann_natom
    ias=wann_iprj(1,i)
    write(60,*)
    write(60,'("ias : ",I4)')ias
    do l=0,lmaxlu
      if (sum(abs(wf_ene_mtrx(idxlm(l,-l):idxlm(l,l),idxlm(l,-l):idxlm(l,l),:,:,ias))).gt.1d-8) then
! occupancy matrix
        write(60,'("  occupancy matrix")')
        t=0.0
        do ispn=1,nspinor
          write(60,'("  ispn : ",I1)')ispn
          write(60,'("    real part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dreal(wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
            t(ispn)=t(ispn)+dreal(wf_occ_mtrx(lm1,lm1,ispn,ispn,ias))
          enddo
          write(60,'("    imag part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dimag(wf_occ_mtrx(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          write(60,'("    occupancy : ",F12.6)')t(ispn)
        enddo !ispn
        write(60,'("  total occupancy : ",F12.6)')sum(t)
        if (nspinor.eq.2) then
          write(60,'("  moment : ",F12.6)')t(1)-t(2)
        endif      
! energy matrix
        write(60,*)
        write(60,'("  energy matrix")')
        do ispn=1,nspinor
          write(60,'("  ispn : ",I1)')ispn
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
            if (abs(wf_ene_mtrx(lm1,lm1,ispn,ispn,ias)).gt.1d-12) then
              n=n+1
              w2=w2+dreal(wf_ene_mtrx(lm1,lm1,ispn,ispn,ias))
            endif
          enddo
          if (n.ne.0) write(60,'("    average energy : ",F12.6)')w2/n
        enddo !ispn
      endif
    enddo !l
  enddo !i
endif
deallocate(wf_ene_mtrx,wf_occ_mtrx)
return
end
  
