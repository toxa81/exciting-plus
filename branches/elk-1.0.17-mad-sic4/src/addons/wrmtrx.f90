subroutine wrmtrx(fname,nr,nc,a,lda)
implicit none
character*(*), intent(in) :: fname
integer, intent(in) :: nr
integer, intent(in) :: nc
integer, intent(in) :: lda
complex(8), intent(in) :: a(lda,*)
integer i,j,fout
if (fname.eq."") then
  fout=6
else
  fout=153
  open(fout,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
endif
do j=1,nc
  write(fout,'("j : ",I10)')j
  do i=1,nr
    write(fout,'(4X,"i : ",I10,"   a(i,j) : ",3G20.12)')i,dreal(a(i,j)),&
      dimag(a(i,j)),abs(a(i,j))
  enddo
enddo
if (fout.ne.6) close(fout)
return
end
