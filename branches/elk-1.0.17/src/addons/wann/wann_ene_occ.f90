subroutine wann_ene_occ
use modmain
use modldapu
use mod_mpi_grid
implicit none
complex(8), allocatable :: wann_ene_m(:,:,:,:,:)
complex(8), allocatable :: wann_occ_m(:,:,:,:,:)
allocate(wann_ene_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
allocate(wann_occ_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
call wann_ene_occ_(wann_ene_m,wann_occ_m)
deallocate(wann_ene_m,wann_occ_m)
return
end

subroutine wann_ene_occ_(wann_ene_m,wann_occ_m)
use modmain
use modldapu
use mod_mpi_grid
implicit none
complex(8), intent(inout) :: wann_ene_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot)
complex(8), intent(inout) :: wann_occ_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot)
! local variables
real(8) t(2),w2
integer n,i,ik,ispn,ias,lm1,lm2,l,j,n1,n2,m1,m2,ispn1,ispn2,ikloc
complex(8) z2

wann_ene=0.d0
wann_occ=0.d0
do n=1,nwann
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    do i=1,nstsv
      w2=dreal(dconjg(wann_c(n,i,ikloc))*wann_c(n,i,ikloc))
      wann_ene(n)=wann_ene(n)+w2*evalsv(i,ik)*wkpt(ik)
      wann_occ(n)=wann_occ(n)+w2*occsv(i,ik)*wkpt(ik)
    enddo
  enddo
enddo
call mpi_grid_reduce(wann_ene(1),nwann,dims=(/dim_k/))
call mpi_grid_reduce(wann_occ(1),nwann,dims=(/dim_k/))
if (wproc.and.nosym) then
  write(60,*)
  write(60,'(" WF  energy (Ha)  occupancy ")')
  write(60,'("----------------------------")')
  do n=1,nwann
    write(60,'(I4,2F12.6)')n,wann_ene(n),wann_occ(n)
  enddo
endif

wann_occ_m=zzero
wann_ene_m=zzero
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
          wann_occ_m(lm1,lm2,ispn1,ispn2,ias)=&
            wann_occ_m(lm1,lm2,ispn1,ispn2,ias)+&
            z2*occsv(j,ik)*wkpt(ik)
          wann_ene_m(lm1,lm2,ispn1,ispn2,ias)=&
            wann_ene_m(lm1,lm2,ispn1,ispn2,ias)+&
            z2*evalsv(j,ik)*wkpt(ik)
        enddo !j
      enddo !ikloc
    endif
  enddo
enddo 
call mpi_grid_reduce(wann_occ_m(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(wann_ene_m(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),all=.true.)
! convert from Rlm to Ylm basis
do ias=1,natmtot
  do ispn=1,nspinor
    call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),wann_occ_m(1,1,ispn,ispn,ias))
    call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),wann_ene_m(1,1,ispn,ispn,ias))
  enddo
enddo
! symmetrise matrix
call symdmat(lmaxlu,lmmaxlu,wann_occ_m)
call symdmat(lmaxlu,lmmaxlu,wann_ene_m)
! convert back to Rlm
do ias=1,natmtot
  do ispn=1,nspinor
    call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),wann_occ_m(1,1,ispn,ispn,ias))
    call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),wann_ene_m(1,1,ispn,ispn,ias))
  enddo
enddo
do n=1,nwann
  ias=iwann(1,n)
  lm1=iwann(2,n)
  ispn1=iwann(3,n)
  wann_ene(n)=wann_ene_m(lm1,lm1,ispn1,ispn1,ias)
enddo
if (wproc) then
  write(60,*)
  write(60,'("On-site matrices in WF basis")')
  do i=1,wann_natom
    ias=wann_iprj(1,i)
    write(60,*)
    write(60,'("ias : ",I4)')ias
    do l=0,lmaxlu
      if (sum(abs(wann_ene_m(idxlm(l,-l):idxlm(l,l),idxlm(l,-l):idxlm(l,l),:,:,ias))).gt.1d-8) then
! occupancy matrix
        write(60,'("  occupancy matrix")')
        t=0.0
        do ispn=1,nspinor
          write(60,'("  ispn : ",I1)')ispn
          write(60,'("    real part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dreal(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
            t(ispn)=t(ispn)+dreal(wann_occ_m(lm1,lm1,ispn,ispn,ias))
          enddo
          write(60,'("    imag part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dimag(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
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
            write(60,'(2X,7F12.6)')(dreal(wann_ene_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          write(60,'("    imag part")')
          do lm1=l**2+1,(l+1)**2
            write(60,'(2X,7F12.6)')(dimag(wann_ene_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
          enddo
          w2=0.d0
          n=0
          do lm1=l**2+1,(l+1)**2
            if (abs(wann_ene_m(lm1,lm1,ispn,ispn,ias)).gt.1d-12) then
              n=n+1
              w2=w2+dreal(wann_ene_m(lm1,lm1,ispn,ispn,ias))
            endif
          enddo
          if (n.ne.0) write(60,'("    average energy : ",F12.6)')w2/n
        enddo !ispn
      endif
    enddo !l
  enddo !i
endif
return
end
  
