program sicvlm
use mod_hdf5
implicit none
integer lmmaxwan,s_nr_min,nspinor,s_nr,n
integer nrtot,nrxyz(3),recsz,i,ir,i1,i2,i3
real(8) orig(3),bound3d(3,3),zero3d(3),val,vrc(3)
real(4), allocatable :: vyz(:)
real(8), allocatable :: vlm(:,:,:)
real(8), allocatable :: s_r(:)
complex(8), allocatable :: wlm(:,:)
complex(8) zval
character*100 fname
call hdf5_initialize
fname="vlm__0001.hdf5"
call hdf5_read(trim(fname),"/","lmmaxwan",lmmaxwan)
call hdf5_read(trim(fname),"/","s_nr_min",s_nr_min)
call hdf5_read(trim(fname),"/","nspinor",nspinor)
call hdf5_read(trim(fname),"/","s_nr",s_nr)
allocate(vlm(lmmaxwan,s_nr_min,nspinor))
allocate(s_r(s_nr))
call hdf5_read(trim(fname),"/","vlm",vlm(1,1,1),&
  (/lmmaxwan,s_nr_min,nspinor/))
call hdf5_read(trim(fname),"/","s_r",s_r(1),(/s_nr/))

allocate(wlm(lmmaxwan,s_nr))
fname="sic.hdf5"
call hdf5_read(trim(fname),"/wann/n0001/s0001","wlm",wlm(1,1),(/lmmaxwan,s_nr/))

nrxyz=(/1200,1200,5/)
zero3d=0.d0
bound3d(:,1)=(/10.d0,0.d0,0.d0/)
bound3d(:,2)=(/0.d0,10.d0,0.d0/)
bound3d(:,3)=(/0.d0,0.d0,10.d0/)

nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
n=1
recsz=4*nrxyz(2)*nrxyz(3)

allocate(vyz(nrxyz(2)*nrxyz(3)))

write(fname,'("sicv",I3.3,".dx")')n
open(70,file=trim(fname),status="REPLACE",form="FORMATTED")
write(70,'("object 1 class gridpositions counts",3I6)')nrxyz(1),nrxyz(2),nrxyz(3)
write(70,'("origin ",3G18.10)')orig(:)
do i=1,3
  write(70,'("delta ",3G18.10)')bound3d(:,i)/nrxyz(i)
enddo
write(70,'("object 2 class gridconnections counts",3I6)')nrxyz(1),nrxyz(2),nrxyz(3)
write(70,'("object 3 class array type float rank 1 shape 1 items ",&
  &I8," lsb ieee data file sicv",I3.3,".bin,0")')nrxyz(1)*nrxyz(2)*nrxyz(3),n
write(70,'("object ""regular positions regular connections"" class field")')
write(70,'("component ""positions"" value 1")')
write(70,'("component ""connections"" value 2")')
write(70,'("component ""data"" value 3")')
write(70,'("end")')
close(70)
write(fname,'("sicv",I3.3,".bin")')n
open(70,file=trim(fname),form="unformatted",status="replace",access="direct",recl=recsz)
do i1=0,nrxyz(1)-1
  write(*,'("slab ",I4," out of ",I4)')i1+1,nrxyz(1)
  ir=1
  do i2=0,nrxyz(2)-1
    do i3=0,nrxyz(3)-1
      vrc(:)=orig(:)+i1*bound3d(:,1)/nrxyz(1)+&
                     i2*bound3d(:,2)/nrxyz(2)+&
                     i3*bound3d(:,3)/nrxyz(3)
      !call func_val(vrc,lmmaxwan,s_nr_min,s_nr,s_r,vlm(1,1,1),val)
      call func_zval(vrc,lmmaxwan,s_nr_min,s_nr,s_r,wlm(1,1),zval)
      vyz(ir)=sngl(abs(zval)**2)
      ir=ir+1
    enddo
  enddo
  write(70,rec=i1+1)vyz
enddo
close(70)


return
end

subroutine func_val(x,lmmaxwan,s_nr_min,s_nr,s_r,vlm,val)
implicit none
real(8), intent(in) :: x(3)
integer, intent(in) :: lmmaxwan
integer, intent(in) :: s_nr_min
integer, intent(in) :: s_nr
real(8), intent(in) :: s_r(s_nr)
real(8), intent(in) :: vlm(lmmaxwan,s_nr_min)
real(8), intent(out) :: val
!
integer lmax,ir1,ir,lm
real(8) rlm(lmmaxwan) 
real(8) x0,tp(2),dx
!
lmax=int(sqrt(dble(lmmaxwan)+0.1d0))-1
if (sum(x(:)**2).gt.(s_r(s_nr_min)**2)) then
  val=0.d0
  return
endif

call sphcrd(x,x0,tp)
call genrlm(lmax,tp,rlm)

ir1=0
do ir=s_nr_min-1,1,-1
  if (s_r(ir).le.x0) then
    ir1=ir
    exit
  endif
enddo
if (ir1.eq.0) then
  ir1=1
  dx=0.d0
else
  dx=(x0-s_r(ir1))/(s_r(ir1+1)-s_r(ir1))
endif
val=0.d0
do lm=1,lmmaxwan
  val=val+(vlm(lm,ir1)+dx*(vlm(lm,ir1+1)-vlm(lm,ir1)))*rlm(lm)
enddo
return
end subroutine

subroutine func_zval(x,lmmaxwan,s_nr_min,s_nr,s_r,zflm,zval)
implicit none
real(8), intent(in) :: x(3)
integer, intent(in) :: lmmaxwan
integer, intent(in) :: s_nr_min
integer, intent(in) :: s_nr
real(8), intent(in) :: s_r(s_nr)
complex(8), intent(in) :: zflm(lmmaxwan,s_nr)
complex(8), intent(out) :: zval
!
integer lmax,ir1,ir,lm
complex(8) ylm(lmmaxwan) 
real(8) x0,tp(2),dx
!
lmax=int(sqrt(dble(lmmaxwan)+0.1d0))-1
zval=dcmplx(0.d0,0.d0)
if (sum(x(:)**2).gt.(s_r(s_nr)**2)) then
  return
endif

call sphcrd(x,x0,tp)
call genylm(lmax,tp,ylm)

ir1=0
do ir=s_nr-1,1,-1
  if (s_r(ir).le.x0) then
    ir1=ir
    exit
  endif
enddo
if (ir1.eq.0) then
  ir1=1
  dx=0.d0
else
  dx=(x0-s_r(ir1))/(s_r(ir1+1)-s_r(ir1))
endif
do lm=1,lmmaxwan
  zval=zval+(zflm(lm,ir1)+dx*(zflm(lm,ir1+1)-zflm(lm,ir1)))*ylm(lm)
enddo
return
end subroutine


