subroutine wann_ene_occ
use modmain
use modwann
implicit none

real(8) wf_ene(wf_dim)
real(8) wf_occ(wf_dim)
integer n,i,ik

wf_ene=0.d0
wf_occ=0.d0

do n=1,wf_dim
  do ik=1,nkpt
    do i=1,nstfv
      wf_ene(n)=wf_ene(n)+dconjg(a_ort(n,i,1,ik))*a_ort(n,i,1,ik)*evalsv(i,ik)*wkpt(ik)
      wf_occ(n)=wf_occ(n)+dconjg(a_ort(n,i,1,ik))*a_ort(n,i,1,ik)*occsv(i,ik)*wkpt(ik)
    enddo
  enddo
enddo
write(60,*)
write(60,'("Wannier functions")')
write(60,'(" wf  energy (Ha)  energy (eV)    occupancy ")')
write(60,'("-------------------------------------------")')
do n=1,wf_dim
  write(60,'(1X,I2,1X,F12.6,1X,F12.6,1X,F12.6)')n,wf_ene(n),wf_ene(n)*ha2ev,wf_occ(n)
enddo
write(60,*)

return
end
  
