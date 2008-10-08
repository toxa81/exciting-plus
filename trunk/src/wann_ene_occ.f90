subroutine wann_ene_occ
use modmain
use modwann
implicit none

real(8) wf_ene(wf_dim),wf_occ(wf_dim)
integer n,i,ik
integer, external :: ikglob

wf_ene=0.d0
wf_occ=0.d0

do n=1,wf_dim
  do ik=1,nkptloc(iproc)
    do i=1,nstfv
      wf_ene(n)=wf_ene(n)+dconjg(wfc(n,i,1,ik))*wfc(n,i,1,ik) * &
        evalsv(i,ikglob(ik))*wkpt(ikglob(ik))
      wf_occ(n)=wf_occ(n)+dconjg(wfc(n,i,1,ik))*wfc(n,i,1,ik) * &
        occsv(i,ikglob(ik))*wkpt(ikglob(ik))
    enddo
  enddo
enddo
call dsync(wf_ene,wf_dim,.true.,.false.)
call dsync(wf_occ,wf_dim,.true.,.false.)
if (iproc.eq.0) then
  write(60,*)
  write(60,'("Wannier functions")')
  write(60,'(" wf  energy (Ha)  energy (eV)    occupancy ")')
  write(60,'("-------------------------------------------")')
  do n=1,wf_dim
    write(60,'(1X,I2,1X,F12.6,1X,F12.6,1X,F12.6)')n,wf_ene(n),wf_ene(n)*ha2ev,wf_occ(n)
  enddo
  write(60,*)
endif

return
end
  
