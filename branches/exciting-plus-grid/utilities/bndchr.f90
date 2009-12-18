program main
implicit none
      
integer lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
integer lm,ias,ispn,ist,ik,spin
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: dpp1d(:),e(:,:),w(:,:),lines(:,:)
integer, allocatable :: orb(:,:),tmp(:),wf1(:),wf2(:)
integer i,j,ist1,n,iin
real(8) scale,emin,emax,wmax,wmax1,efermi
real(8), parameter :: ha2ev = 27.21138386d0
logical wannier,l1,unitsev
integer wf_dim
real(8), allocatable :: wfc(:,:,:,:)
character*100 str
      
open(50,file='BNDCHR.OUT',form='formatted',status='old')
read(50,*)lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
read(50,*)efermi
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(dpp1d(nkpt))
allocate(e(nstsv,nkpt))
allocate(w(nstsv,nkpt))
allocate(orb(lmmax,natmtot))
allocate(tmp(lmmax))
      
do ik=1,nkpt
  read(50,*)dpp1d(ik)
  read(50,*)(e(ist,ik),ist=1,nstsv)
  read(50,*)((((bndchr(lm,ias,ispn,ist,ik),lm=1,lmmax), &
                ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
enddo
read(50,*)wannier
if (wannier) then
  read(50,*)wf_dim
  allocate(wfc(wf_dim,nstfv,nspinor,nkpt))
  do ik=1,nkpt
    read(50,*)((wfc(n,i,1,ik),n=1,wf_dim),i=1,nstfv)
  enddo
endif
close(50)

35 continue
write(*,'("Put Fermi level to zero? (1-Yes,2-No)")')
read(*,*)iin
if (.not.(iin.eq.1.or.iin.eq.2)) goto 35
if (iin.eq.1) then
  e=e-efermi
  efermi=0.d0
endif

40 continue
write(*,'("Energy units? (1-Ha,2-eV)")')
read(*,*)iin
if (.not.(iin.eq.1.or.iin.eq.2)) goto 40
if (iin.eq.2) e=e*ha2ev
 
emin=minval(e(:,:))
emax=maxval(e(:,:))
write(*,'("Band energy range : ",2F8.2)')emin,emax
write(*,'("Input energy interval (emin emax)")')
read(*,*)emin,emax
      
spin=1
if (nspinor.eq.2) then
  write(*,'("Input spin direction (1 or 2)")')
  read(*,*)spin
endif

      
      allocate(lines(4,nlines))
      
      l1=.false.
      if (wannier) then
        write(*,*)'Do you want band character of WFs (T) or lm- orbitals (F)'
        read(*,*)l1
      endif
      if (l1) goto 20  
      
      
      orb = 0
      write(*,*)'Number of atoms: ',natmtot
      do while (.true.)
        i = -1
        write(*,*)'Input atom number (/ to proceed):'
        read(*,*)i
        if (i.eq.-1) goto 10
        if (i.ge.1.and.i.le.natmtot) then
          tmp = 0
          write(*,*)'Input orbitals:'
          read(*,*)tmp
          do j = 1, lmmax
            if (tmp(j).ge.1.and.tmp(j).le.lmmax) then
              orb(tmp(j),i) = 1
            endif
          enddo
        endif
      enddo
      
  10  continue
      do ias = 1, natmtot
        write(*,*)(orb(lm,ias),lm=1,lmmax)
      enddo
      
      
      w = 0.d0
      do ik = 1, nkpt
      do ist = 1, nstsv
        do lm = 1, lmmax
        do ias = 1, natmtot
         do ispn=1,nspinor
            w(ist,ik) = w(ist,ik)+bndchr(lm,ias,ispn,ist,ik)
          enddo
        enddo
        enddo
      enddo
      enddo
      wmax = maxval(w(:,:))
      
      w = 0.d0
      do ik = 1, nkpt
      do ist = 1, nstsv
        do lm = 1, lmmax
        do ias = 1, natmtot
          do ispn=1,nspinor
            w(ist,ik) = w(ist,ik)+orb(lm,ias)*bndchr(lm,ias,ispn,ist,ik)
          enddo
        enddo
        enddo
      enddo
      enddo
      
      goto 30
      
      
  20  continue
  
      allocate(wf1(wf_dim),wf2(wf_dim))
      wf1=0
      wf2=0
      write(*,*)'Input WFs'
      read(*,*)wf2
      do i=1,wf_dim
        if (wf2(i).gt.0.and.wf2(i).le.wf_dim) wf1(wf2(i))=1
      enddo
      write(*,*)(wf1(i),i=1,wf_dim)
      
      ispn=1
      w=0.d0
      do ik = 1, nkpt
      do ist = 1, nstfv
      do n=1,wf_dim
        w(ist,ik)=w(ist,ik)+wfc(n,ist,ispn,ik)**2
      enddo
      enddo
      enddo
      wmax = maxval(w(:,:))

      ispn=1
      w=0.d0
      do ik = 1, nkpt
      do ist = 1, nstfv
      do n=1,wf_dim
        w(ist,ik)=w(ist,ik)+wfc(n,ist,ispn,ik)**2*wf1(n)
      enddo
      enddo
      enddo
      
    30 continue
      
      
      
      scale = 1.d0
      write(*,*)'Input maximum line width [eV]'
      read(*,*)scale
      
      scale = 0.5d0*scale/wmax
      
      open(50,file='BNDS.DAT',form='formatted',status='replace')
      open(51,file='BNDS1.DAT',form='formatted',status='replace')
      open(52,file='BNDS2.DAT',form='formatted',status='replace')
      open(53,file='BNDS3.DAT',form='formatted',status='replace')
      open(54,file='BNDS_xydy.DAT',form='formatted',status='replace')
      do ist = 1,nstfv
        ist1=ist+(spin-1)*nstfv
        do ik = 1, nkpt
          write(50,*)dpp1d(ik),e(ist1,ik)
          write(51,*)dpp1d(ik),e(ist1,ik)+scale*w(ist1,ik)
          write(52,*)dpp1d(ik),e(ist1,ik)-scale*w(ist1,ik)
          write(53,*)dpp1d(ik),e(ist1,ik)+scale*w(ist1,ik)
          write(53,*)dpp1d(ik),e(ist1,ik)-scale*w(ist1,ik)
          write(53,*)
          write(54,*)dpp1d(ik),e(ist1,ik),2*scale*w(ist1,ik)
        enddo
        write(50,*)
        write(51,*)
        write(52,*)
        write(54,*)
      enddo
      close(50) 
      close(51) 
      close(52) 
      close(53) 
      close(54)
      
      open(50,file='BNDS.GNU',form='formatted',status='replace')
      write(50,*)'set term postscript portrait'
      write(50,*)'set noxzeroaxis'
      write(50,*)'set tics out'
      write(50,*)'set noxtics'
      write(50,*)'set nokey'
      write(50,'(" set yrange [",F12.6,":",F12.6,"]")')emin,emax      
      write(50,'(" set xrange [",F12.6,":",F12.6,"]")')0.d0,dpp1d(nkpt)
      write(50,*)"set ylabel 'Energy (eV)'"
      write(50,*)"set output 'bnds.ps'"      
      write(50,'("plot ''BNDS.DAT'' with line 1, ''BNDS1.DAT'' with line 1, &
                 & ''BNDS2.DAT'' with line 1, ''BNDS3.DAT'' with line 1, &
                 & ''BNDS4.DAT'' with line 2")')
      close(50)
      
      open(50,file='BANDLINES.OUT',form='formatted',status='old')
      do i = 1, nlines
        read(50,*)lines(1,i),lines(2,i)
        read(50,*)lines(3,i),lines(4,i)
        read(50,*)
      enddo
      close(50)
      
      open(50,file='BNDS4.DAT',form='formatted',status='replace')
      do i = 1, nlines
        write(50,*)lines(1,i),emin+1.d-4
        write(50,*)lines(3,i),emax-1.d-4
        write(50,*)
      enddo
      write(50,*)0.d0,efermi
      write(50,*)dpp1d(nkpt),efermi
      close(50)
     
      
      end
      