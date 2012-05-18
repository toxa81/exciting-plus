subroutine printmegqblh(iq)
use modmain
use mod_wannier
use mod_expigqr
implicit none
integer, intent(in) :: iq
integer i,ig,ikloc,ik,iboffs
integer ist1,ist2
character*100 fname
write(fname,'("MEGQBLH_iq_",I4.4,".OUT")')iq
if (iproc.eq.0) then
  open(200,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    write(200,'("ikloc : ",I4,"  idxkq : ",2I4)')ikloc,idxkq(:,ik)
    do ist2=1,nstsv
      write(200,'("  ist2 : ",I4)')ist2
      do ist1=1,nstsv
        do i=1,nmegqblh(ikloc)
          if (bmegqblh(1,i,ikloc).eq.ist1.and.&
              bmegqblh(2,i,ikloc).eq.ist2) then
            write(200,'("    ist1 : ",I4)')ist1
            do ig=1,ngvecme
              write(200,'("      ig : ",I4,"   ",3G18.10)')&
                ig,dreal(megqblh(i,ig,ikloc)),-dimag(megqblh(i,ig,ikloc)),&
                abs(megqblh(i,ig,ikloc))**2
            enddo
          endif
        enddo
      enddo
    enddo
  enddo !ikloc
  close(200)
endif
return
end