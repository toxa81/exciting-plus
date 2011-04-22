subroutine sic_me(fout)
use modmain
use mod_sic
use mod_hdf5
implicit none
integer, intent(in) :: fout
integer i,j,n,i1,j1,n1,i2,j2,n2,ispn,vl(3),nwtloc,iloc,ik
logical texist
real(8) t1,t2,t3,pos1(3),pos2(3),vtrc(3)
complex(8) me1,me2,z1
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: vwank(:,:)
!
if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("matrix elements")')
  write(fout,'(80("="))')
endif
! measure time fo one me
pos1=0.d0
pos2=1.d0 
call timer_start(t_sic_me,reset=.true.)
me1=s_spinor_dot_ll(pos1,pos2,s_wvlm(1,1,1,1),s_wanlm(1,1,1,1))
call timer_stop(t_sic_me)
if (wproc) then
  write(fout,*)
  write(fout,'("time for one matrix element : ",F8.3," sec.")')&
    timer_get_value(t_sic_me)
endif
call timer_start(t_sic_me,reset=.true.)
! read old matrix elements
allocate(vwanme_old(sic_wantran%nwt))
vwanme_old=zzero
inquire(file="sic.hdf5",exist=texist)
if (texist) then
  call hdf5_read("sic.hdf5","/","nwt",i)
  if (i.eq.sic_wantran%nwt) then
    call hdf5_read("sic.hdf5","/","vwanme",vwanme_old(1),(/sic_wantran%nwt/))
  endif
endif
! compute matrix elements of SIC potential
!  vwanme = <(W*V)_n|W_{n1,T}>
vwanme=zzero
nwtloc=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/))
do iloc=1,nwtloc
  i=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/),loc=iloc)
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  !do ispn=1,nspinor
  !  vwanme(i)=vwanme(i)+s_dot_ll(pos1,pos2,s_wvlm(1,1,ispn,j),s_wanlm(1,1,ispn,j1))
  !enddo
  vwanme(i)=s_spinor_dot_ll(pos1,pos2,s_wvlm(1,1,1,j),s_wanlm(1,1,1,j1))
enddo
call mpi_grid_reduce(vwanme(1),sic_wantran%nwt,all=.true.)
! check localization criterion 
t2=-1.d0
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  t1=abs(vwanme(i)-dconjg(vwanme(j)))
  if (t1.ge.t2) then
    t2=t1
    i1=i
    j1=j
    me1=vwanme(i)
    me2=vwanme(j)
  endif
enddo
! symmetrize the potential matrix elements
do i=1,sic_wantran%nwt
  n=sic_wantran%iwt(1,i)
  n1=sic_wantran%iwt(2,i)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n1,n,-vl(1),-vl(2),-vl(3))
  z1=0.5d0*(vwanme(i)+dconjg(vwanme(j)))
  vwanme(i)=z1
  vwanme(j)=dconjg(z1)
enddo
! compute RMS difference
t3=0.d0
do i=1,sic_wantran%nwt
  t3=t3+abs(vwanme(i)-vwanme_old(i))**2
enddo
deallocate(vwanme_old)
call timer_stop(t_sic_me)
if (wproc) then
  write(fout,*)
  write(fout,'("maximum deviation from ""localization criterion"" : ",F12.6)')t2
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
    write(fout,'("  n : ",I4,8X,2G18.10)')n,dreal(vwanme(i)),dimag(vwanme(i))
  enddo  
  t3=sqrt(t3/sic_wantran%nwt)
  write(fout,*)
  write(fout,'("matrix elements RMS difference :",G18.10)')t3
  write(fout,*)
  write(fout,'("done in : ",F8.3," sec.")')timer_get_value(t_sic_me)
  write(fout,*)
  call flushifc(fout)
endif
!! check hermiticity of V_nn'(k)
!allocate(vwank(sic_wantran%nwan,sic_wantran%nwan))
!do ik=1,nkpt
!  vwank=zzero
!  do i=1,sic_wantran%nwt
!    n1=sic_wantran%iwt(1,i)
!    j1=sic_wantran%idxiwan(n1)
!    n2=sic_wantran%iwt(2,i)
!    j2=sic_wantran%idxiwan(n2)
!    vl(:)=sic_wantran%iwt(3:5,i)
!    vtrc(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
!    z1=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
!    vwank(j1,j2)=vwank(j1,j2)+z1*vwanme(i)
!  enddo
!  t1=0.d0
!  do j1=1,sic_wantran%nwan
!    do j2=1,sic_wantran%nwan
!      t1=max(t1,abs(vwank(j1,j2)-dconjg(vwank(j2,j1))))
!    enddo
!  enddo
!  if (wproc) then
!    write(fout,'("ik : ",I4,"   max.herm.err : ",G18.10 )')ik,t1
!  endif
!enddo
!deallocate(vwank)
return
end
