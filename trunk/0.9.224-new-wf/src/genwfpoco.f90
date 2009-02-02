subroutine genwfpoco(ik,ngp,igpig,evecfv,evecsv,apwalm,wann_poco)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wann_poco(nstsv,nstsv)

integer i,j,n,ispn,istfv,itype
complex(8) zt1
complex(8), allocatable :: wffvwann(:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)

! calculate first-fariatinal states
allocate(wffvmt(nstfv,nrfmax,lmmaxvr,natmtot))
call genwffvmt(lmaxvr,lmmaxvr,ngp,evecfv,apwalm,wffvmt)
! calculate <u_{n,k,\sigma}|\phi_{i}>, where |u_{n,k,\sigma}> are Fourier-transforms of 
!  Wannier functions and |\phi_{i}> are first-variational states
allocate(wffvwann(nstfv,wann_nmax,wann_nspin))
call genwffvwann(ik,ngp,igpig,wffvmt,evecfv,wffvwann)

do i=1,nstfv
  do j=1,nstfv
    do ispn=1,wann_nspin
      zt1=dcmplx(0.d0,0.d0)
      do n=1,nwann(ispn)
        itype=iwann(n,ispn,4)
        zt1=zt1+wffvwann(i,n,ispn)*dconjg(wffvwann(j,n,ispn))*wann_v(itype)
      enddo
      wann_poco(i+(ispn-1)*nstfv,j+(ispn-1)*nstfv)=zt1
    enddo
  enddo
enddo
deallocate(wffvmt,wffvwann)


return
end

!subroutine writewfpoco
!use modmain
!use modwann
!implicit none
!
!integer ik,i,j1,j2
!character*4 c4
!integer, external :: ikglob
!
!if (iproc.eq.0) then
!  call write_integer(nkpt,1,'/dimensions','nkpt')
!  call write_integer(nstsv,1,'/dimensions','nstsv')
!endif
!  
!do i=0,nproc-1
!  if (i.eq.iproc) then
!    do ik=1,nkptloc(iproc)
!      write(c4,'(I4.4)')ikglob(ik)
!      call write_complex16(wfpoco(1,1,ik),nstsv*nstsv,'/kpoints/'//c4,'wfpoco')
!      call write_real8(vkl(1,ikglob(ik)),3,'/kpoints/'//c4,'vkl')
!    enddo
!  endif
!  call barrier
!enddo
!return
!end
!
!subroutine readwfpoco
!use modmain
!use modwann
!implicit none
!
!integer ik,i,j1,j2,nkpt_,nstsv_
!integer, external :: ikglob
!character*4 c4
!
!do i=0,nproc-1
!  if (i.eq.iproc) then
!    do ik=1,nkptloc(iproc)
!      write(c4,'(I4.4)')ikglob(ik)
!      call read_complex16(wfpoco(1,1,ik),nstsv*nstsv,'/kpoints/'//c4,'wfpoco')
!    enddo
!  endif
!  call barrier
!enddo
!return
!end
!
!
!
