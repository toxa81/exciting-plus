subroutine wann_ene_occ
use modmain
implicit none

real(8) wf_ene(wann_nmax,wann_nspin),wf_occ(wann_nmax,wann_nspin)
integer n,i,ik,ispn
integer, external :: ikglob
real(8) w2

wf_ene=0.d0
wf_occ=0.d0

do ispn=1,wann_nspin
  do n=1,nwann(ispn)
    do ik=1,nkptloc(iproc)
      do i=1,nstfv
        w2=dconjg(wann_c(n,i,ispn,ik))*wann_c(n,i,ispn,ik)
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


return
end
  
