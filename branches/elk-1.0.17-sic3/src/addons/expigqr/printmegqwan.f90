subroutine printmegqwan(iq)
use modmain
implicit none
integer, intent(in) :: iq
integer j,ig
integer n1,n2
character*100 fname
write(fname,'("MEGQWAN_iq_",I4.4,".OUT")')iq
if (iproc.eq.0) then
  open(200,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
  do j=1,nmegqwan
    n1=imegqwan(1,j)
    n2=imegqwan(2,j)
    do ig=1,ngvecme
      write(200,'(I4,"    ",I4,"    ",3I4,"    ",3G18.10)')n1,n2,&
        imegqwan(3:5,j),dreal(megqwan(j,ig)),dimag(megqwan(j,ig)),&
        abs(megqwan(j,ig))**2
    enddo !ig      
  enddo !j
  close(200)
endif
return
end
