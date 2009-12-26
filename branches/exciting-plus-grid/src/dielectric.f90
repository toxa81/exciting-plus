
! Copyright (C) 2002-2008 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dielectric
use modmain
!use modtest
implicit none
! local variables
integer ik,jk,isym
integer ist,jst,iw,i,j,l
integer recl,iostat
real(8) eji,wplas,t1,t2
real(8) v1(3),v2(3),v3(3)
complex(8) zv(3),eta,zt1,zt2
character(256) fname
! allocatable arrays
integer, allocatable :: lspl(:)
real(8), allocatable :: w(:)
real(8), allocatable :: delta(:,:,:)
complex(8), allocatable :: pmat(:,:,:,:)
complex(8), allocatable :: sigma(:)
integer ikloc
! external functions
real(8) sdelta
external sdelta
! initialise universal variables
call init0
call init1
! read Fermi energy from file
if (mpi_grid_root()) then
  call readfermi
  do ik=1,nkpt
! get the eigenvalues and occupancies from file
    call getevalsv(vkl(:,ik),evalsv(:,ik))
    call getoccsv(vkl(:,ik),occsv(:,ik))
  end do
endif
call mpi_grid_bcast(efermi)
call mpi_grid_bcast(evalsv(1,1),nstsv*nkpt)
call mpi_grid_bcast(occsv(1,1),nstsv*nkpt)
! allocate local arrays
allocate(lspl(nkptnr))
allocate(w(nwdos))
if (usegdft) allocate(delta(nstsv,nstsv,nkpt))
allocate(pmat(3,nstsv,nstsv,nkptnrloc))
allocate(sigma(nwdos))
! compute generalised DFT correction
if (usegdft) then
  call readstate
  call poteff
  call linengy
  call genapwfr
  call genlofr
  do ik=1,nkpt
    call gdft(ik,delta(:,:,ik))
  end do
end if
! generate energy grid (starting from zero)
t1=wdos(2)/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)
end do
! find crystal symmetries which map non-reduced k-points to reduced equivalents
do ik=1,nkptnr
  call findkpt(vklnr(:,ik),isym,jk)
  lspl(ik)=lsplsymc(isym)
end do
if (mpi_grid_root()) then
! find the record length for momentum matrix element file
  inquire(iolength=recl) pmat(:,:,:,1)
endif
call mpi_grid_bcast(recl)
open(50,file='PMAT.OUT',action='READ',form='UNFORMATTED',access='DIRECT', &
 recl=recl,iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(dielectric): error opening PMAT.OUT")')
  write(*,*)
  call pstop
end if
! i divided by the complex relaxation time
eta=cmplx(0.d0,swidth)
! read pmat
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (mpi_grid_x(dim_k).eq.i) then
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! equivalent reduced k-point
        jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! read momentum matrix elements from direct-access file
        read(50,rec=jk) pmat(:,:,:,ikloc)
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif
! loop over dielectric tensor components
do l=1,noptcomp
  i=optcomp(1,l)
  j=optcomp(2,l)
  wplas=0.d0
  sigma(:)=0.d0
! loop over non-reduced k-points
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! equivalent reduced k-point
    jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! valance states
    do ist=1,nstsv
! conduction states
      do jst=1,nstsv
! rotate the matrix elements from the reduced to non-reduced k-point
! (note that the inverse operation is used)
        v1(:)=dble(pmat(:,ist,jst,ikloc))
        call r3mv(symlatc(:,:,lspl(ik)),v1,v2)
        v1(:)=aimag(pmat(:,ist,jst,ikloc))
        call r3mv(symlatc(:,:,lspl(ik)),v1,v3)
        zv(:)=cmplx(v2(:),v3(:),8)
        zt1=zv(i)*conjg(zv(j))
        eji=evalsv(jst,jk)-evalsv(ist,jk)
        if ((evalsv(ist,jk).le.efermi).and.(evalsv(jst,jk).gt.efermi)) then
! scissors correction
          eji=eji+scissor
! generalised DFT correction
          if (usegdft) eji=eji+delta(jst,ist,jk)
        end if
        t1=occsv(ist,jk)*(1.d0-occsv(jst,jk)/occmax)
        if (abs(t1).gt.epsocc) then
          t2=t1/(eji+swidth)
          do iw=1,nwdos
            sigma(iw)=sigma(iw)+t2*(zt1/(w(iw)-eji+eta) &
             +conjg(zt1)/(w(iw)+eji+eta))
          end do
        end if
! add to the plasma frequency
        if (intraband) then
          if (i.eq.j) then
            if (abs(eji).gt.1.d-8) then
              t2=occsv(jst,jk)-occsv(ist,jk)
              wplas=wplas+wkptnr(ik)*dble(zt1)*t2/eji
            end if
          end if
        end if
      end do
    end do
  end do
  call mpi_grid_reduce(sigma(1),nwdos,dims=(/dim_k/),side=.true.)
  if (mpi_grid_root()) then
    zt1=zi/(omega*dble(nkptnr))
    sigma(:)=zt1*sigma(:)
! intraband contribution
    if (intraband) then
      if (i.eq.j) then
        wplas=sqrt(abs(wplas)*fourpi/omega)
! write the plasma frequency to file
        write(fname,'("PLASMA_",2I1,".OUT")') i,j
        open(60,file=trim(fname),action='WRITE',form='FORMATTED')
        write(60,'(G18.10," : plasma frequency")') wplas
        close(60)
! add the intraband contribution to sigma
        t1=wplas**2/fourpi
        do iw=1,nwdos
          sigma(iw)=sigma(iw)+t1/(swidth-zi*w(iw))
        end do
      end if
    end if
! write the optical conductivity to file
    write(fname,'("SIGMA_",2I1,".OUT")') i,j
    open(60,file=trim(fname),action='WRITE',form='FORMATTED')
    do iw=1,nwdos
      write(60,'(2G18.10)') w(iw)*ha2ev,dble(sigma(iw))*ha2ev
    end do
    write(60,'("     ")')
    do iw=1,nwdos
      write(60,'(2G18.10)') w(iw)*ha2ev,aimag(sigma(iw))*ha2ev
    end do
    close(60)
! write the dielectric function to file
    write(fname,'("EPSILON_",2I1,".OUT")') i,j
    open(60,file=trim(fname),action='WRITE',form='FORMATTED')
    t1=0.d0
    if (i.eq.j) t1=1.d0
    do iw=1,nwdos
      if (w(iw).gt.1.d-8) then
        t2=t1-fourpi*aimag(sigma(iw)/(w(iw)+eta))
        write(60,'(2G18.10)') w(iw)*ha2ev,t2
      end if
    end do
    write(60,'("     ")')
    do iw=1,nwdos
      if (w(iw).gt.1.d-8) then
        t2=fourpi*dble(sigma(iw)/(w(iw)+eta))
        write(60,'(2G18.10)') w(iw)*ha2ev,t2
      end if
    end do
    close(60)
  endif
end do
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("Info(dielectric):")')
  write(*,'(" dielectric tensor written to EPSILON_ij.OUT")')
  write(*,'(" optical conductivity written to SIGMA_ij.OUT")')
  if (intraband) then
    write(*,'(" plasma frequency written to PLASMA_ij.OUT")')
  end if
  write(*,'(" for components")')
  do l=1,noptcomp
    write(*,'("  i = ",I1,", j = ",I1)') optcomp(1:2,l)
  end do
  write(*,*)
endif
deallocate(lspl,w,pmat,sigma)
if (usegdft) deallocate(delta)
return
end subroutine

