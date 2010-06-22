subroutine printmegqblh
use modmain
implicit none
integer i,ig,ikloc
if (iproc.eq.0) then
  open(200,file="MEGQBLH.OUT",form="FORMATTED",status="REPLACE")
  do ikloc=1,1 !nkptnrloc
    write(200,'("ikloc : ",I4)')ikloc
    do i=1,nmegqblh(ikloc)
      write(200,'("  ist1 : ",I4,"   ist2 : ",I4)')bmegqblh(1,i,ikloc),&
        bmegqblh(2,i,ikloc)
      do ig=1,ngvecme
        write(200,'("    ig : ",I4,"   megqblh(i,ig,ikloc) : ",2G18.10)')&
          ig,dreal(megqblh(i,ig,ikloc)),dimag(megqblh(i,ig,ikloc))
      enddo
    enddo
  enddo
  close(200)
endif
return
end