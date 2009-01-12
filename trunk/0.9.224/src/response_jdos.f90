subroutine response_jdos(occsvnr,evalsvnr)
use modmain
implicit none
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
real(8), intent(in) :: evalsvnr(nstsv,nkptnr)
! energy mesh
complex(8), allocatable :: w(:)
integer ik,i,n1,n2,ist1,ist2,ispn1,ispn2
complex(8), allocatable :: f(:)
character*100 fname
integer, external :: iknrglob

! setup energy mesh
nepts=1+maxomega/domega
allocate(w(nepts))
do i=1,nepts
  w(i)=dcmplx(domega*(i-1),lr_eta)/ha2ev
enddo

allocate(f(nepts))

do ispn1=1,nspinor
do ispn2=1,nspinor
  f=dcmplx(0.d0,0.d0)
  do i=1,nepts
    do n1=1,nstfv
      do n2=1,nstfv
        ist1=n1+(ispn1-1)*nstfv
	ist2=n2+(ispn2-1)*nstfv
        do ik=1,nkptnrloc(iproc)
          f(i)=f(i)+(occsvnr(ist1,iknrglob(ik))-occsvnr(ist2,iknrglob(ik)))/(evalsvnr(ist1,iknrglob(ik))-evalsvnr(ist2,iknrglob(ik))+w(i))
        enddo
      enddo
    enddo
  enddo
  call zsync(f,nepts,.true.,.false.)
  if (iproc.eq.0) then
    f=f/nkptnr/ha2ev/(au2ang)**3/omega
    write(fname,'("JDOS_S1=",I1,"_S2=",I1,".OUT")')ispn1,ispn2
    open(50,file=trim(fname),form='formatted',status='replace')
    do i=1,nepts
      write(50,'(2F12.6)')dreal(w(i))*ha2ev,-dimag(f(i))
    enddo
    close(50)
  endif
enddo
enddo

return
end
