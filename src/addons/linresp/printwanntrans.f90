subroutine printwanntrans(fout,dat)
use modmain
use mod_wannier
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: fout
complex(8), intent(in) :: dat(megqwantran%nwt)
integer ias,jas,lm1,lm2,n,n1,vtl(3)
real(8) vtc(3),vrc(3)
character*20 c1,c2,c3,c4
character, parameter :: orb(4)=(/'s','p','d','f'/)
integer i
character*200 str

do i=1,megqwantran%nwt
  n=megqwantran%iwt(1,i)
  n1=megqwantran%iwt(2,i)
  vtl=megqwantran%iwt(3:5,i)
  ias=wan_info(1,n)
  lm1=wan_info(2,n)
  jas=wan_info(1,n1)
  lm2=wan_info(2,n1)
  vtc(:)=avec(:,1)*vtl(1)+avec(:,2)*vtl(2)+avec(:,3)*vtl(3)
  vrc(:)=vtc(:)+atposc(:,ias2ia(jas),ias2is(jas))-&
                atposc(:,ias2ia(ias),ias2is(ias))       
  write(c1,'(I6)')ias2ia(ias)
  write(c2,'(I1)')lm2m(lm1)+lm2l(lm1)+1
  c3="("//trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
    orb(lm2l(lm1)+1)//trim(adjustl(c2))//")"
  write(c1,'(I6)')ias2ia(jas)
  write(c2,'(I1)')lm2m(lm2)+lm2l(lm2)+1
  c4="("//trim(spsymb(ias2is(jas)))//trim(adjustl(c1))//"-"//&
    orb(lm2l(lm2)+1)//trim(adjustl(c2))//")"
  write(c1,'("i : ",I4)')i
  str=trim(adjustl(c1))
  write(c1,'(I6)')n
  str=trim(str)//"   n="//trim(adjustl(c1))//" "//trim(adjustl(c3))
  write(c1,'(I6)')n1
  str=trim(str)//" ->  n'="//trim(adjustl(c1))//" "//trim(adjustl(c4))
  write(c1,'(I6)')vtl(1)
  write(c2,'(I6)')vtl(2)
  write(c3,'(I6)')vtl(3)
  str=trim(str)//" T=("//trim(adjustl(c1))//" "//trim(adjustl(c2))//" "//&
    trim(adjustl(c3))//")"
  write(c1,'(F12.4)')vrc(1)
  write(c2,'(F12.4)')vrc(2)
  write(c3,'(F12.4)')vrc(3)
  str=trim(str)//" R=("//trim(adjustl(c1))//" "//trim(adjustl(c2))//" "//&
    trim(adjustl(c3))//")"
  write(c1,'(F12.4)')sqrt(sum(vrc(:)**2))
  str=trim(str)//" D="//trim(adjustl(c1))
  write(c1,'(G18.10)')dreal(dat(i))
  write(c2,'(G18.10)')dimag(dat(i))
  write(c3,'(G18.10)')abs(dat(i))
  write(fout,'(A)')trim(str)
  str="           Re="//trim(adjustl(c1))//"  Im="//trim(adjustl(c2))//&
    "  Abs="//trim(adjustl(c3))
  write(fout,'(A)')trim(str)
  
!  write(fout,'("i : ",I4," n=",I4," ",A,"  -> ",I4,"T=",3I3," ",A,"    R=",3F12.6,&
!    &"   D=",F12.6,"  |me(G0q)|=",G18.10)')i,megqwantran%iwt(1,i),&
!    trim(c3),megqwantran%iwt(2,i),trim(c4),vtc,sqrt(sum(vtc(:)**2)),dat(i)
enddo
call flushifc(fout)
return
end
