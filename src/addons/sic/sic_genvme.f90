subroutine sic_genvme(fout)
use modmain
use mod_sic
use mod_hdf5
implicit none
integer, intent(in) :: fout
!
integer i,j,n,i1,j1,n1,vl(3)
logical texist
real(8) t1,t2,t3,pos1(3),pos2(3)
complex(8) me1,me2
!
if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("matrix elements")')
  write(fout,'(80("="))')
endif
! measure time for one me
if (sic_me_cutoff.gt.1.d0) then
  pos1=0.d0
  pos2=1.d0 
  call timer_start(t_sic_me,reset=.true.)
  me1=s_spinor_dotp(pos1,pos2,s_wvlm(1,1,1,1),s_wlm(1,1,1,1))
  call timer_stop(t_sic_me)
  if (wproc) then
    write(fout,*)
    write(fout,'("time for one matrix element : ",F8.3," sec.")')&
      timer_get_value(t_sic_me)
    call flushifc(fout)
  endif
endif
call timer_start(t_sic_me,reset=.true.)
call sic_genvme_dotp(.false.)
! check localization criterion 
t2=-1.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  t1=abs(sic_vme(i)-dconjg(sic_vme(j)))
  if (t1.ge.t2) then
    t2=t1
    i1=i
    j1=j
    me1=sic_vme(i)
    me2=sic_vme(j)
  endif
enddo
! symmetrize the potential matrix elements
!do i=1,sic_wantran%nwt
!  n=sic_wantran%iwt(1,i)
!  n1=sic_wantran%iwt(2,i)
!  vl(:)=sic_wantran%iwt(3:5,i)
!  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
!  z1=0.5d0*(sic_vme(i)+dconjg(sic_vme(j)))
!  sic_vme(i)=z1
!  sic_vme(j)=dconjg(z1)
!enddo
! read old matrix elements
!allocate(vme_old(sic_wantran%nwt))
!vme_old=zzero
!inquire(file="sic.hdf5",exist=texist)
!if (texist) then
!  call hdf5_read("sic.hdf5","/","nwt",i)
!  if (i.eq.sic_wantran%nwt) then
!    call hdf5_read("sic.hdf5","/","vme",vme_old(1),(/sic_wantran%nwt/))
!  endif
!endif
! compute RMS difference
t3=0.d0
do i=1,sic_wantran%nwt
  t3=t3+abs(sic_vme(i)-sic_vme_old(i))**2
enddo
t3=sqrt(t3/sic_wantran%nwt)
call timer_stop(t_sic_me)
if (wproc) then
  write(fout,*)
  write(fout,'("maximum deviation from ""localization criterion"" : ",G18.10)')t2
  write(fout,'("matrix elements with maximum difference : ",2I6)')i1,j1
  write(fout,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,i1),dreal(me1),dimag(me1)
  write(fout,'("  n : ",I4,"    n'' : ",I4,"    T : ",3I4,8X,2G18.10)')&
    sic_wantran%iwt(:,j1),dreal(me2),dimag(me2)
  write(fout,*)
  write(fout,'("diagonal matrix elements (<(W*V)_n|W_n>) :")')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    i=sic_wantran%iwtidx(n,n,0,0,0)
    write(fout,'("  n : ",I4,8X,2G18.10)')n,dreal(sic_vme(i)),dimag(sic_vme(i))
  enddo  
  write(fout,*)
  write(fout,'("matrix elements RMS difference :",G18.10)')t3
  write(fout,*)
  write(fout,'("done in : ",F8.3," sec.")')timer_get_value(t_sic_me)
  write(fout,*)
  call flushifc(fout)
endif
xml_info%sic_vme_rms=t3
xml_info%sic_vme_err=t2
return
end
