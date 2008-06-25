      program  ebands
      implicit none
      
      integer              :: lmmax,natmtot,nspinor,nstsv,nkpt,nlines
      integer              :: lm,ias,ispn,ist,ik
      real*8  ,allocatable :: bndchr(:,:,:,:,:)
      real*8  ,allocatable :: dpp1d(:),e(:,:),w(:,:),lines(:,:)
      integer ,allocatable :: orb(:,:),tmp(:)
      integer              :: i,j
      real*8               :: scale,emin,emax,wmax,wmax1
      real*8, parameter    :: ha2ev = 27.21138386d0            
      
      open(50,file='BANDS.OUT',form='formatted',status='old')
      read(50,*)lmmax,natmtot,nspinor,nstsv,nkpt,nlines
      allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
      allocate(dpp1d(nkpt))
      allocate(e(nstsv,nkpt))
      allocate(w(nstsv,nkpt))
      allocate(orb(lmmax,natmtot))
      allocate(tmp(lmmax))
      
      do ik = 1, nkpt
        read(50,*)dpp1d(ik)
        read(50,*)(e(ist,ik),ist=1,nstsv)
        read(50,*)((((bndchr(lm,ias,ispn,ist,ik),lm=1,lmmax), &
                      ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
      enddo
      close(50)
      
      allocate(lines(4,nlines))
      
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
      
      e = e*ha2ev
      
      emin = minval(e(:,:))
      emax = maxval(e(:,:))
      write(*,*)'Input energy interval (emin emax) [eV]'
      read(*,*)emin,emax
      
      ispn = 1
      
      w = 0.d0
      do ik = 1, nkpt
      do ist = 1, nstsv
        do lm = 1, lmmax
        do ias = 1, natmtot
          w(ist,ik) = w(ist,ik) + bndchr(lm,ias,ispn,ist,ik)
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
          w(ist,ik) = w(ist,ik) + orb(lm,ias)*bndchr(lm,ias,ispn,ist,ik)
        enddo
        enddo
      enddo
      enddo
      
      !wmax = maxval(w(:,:))
      !wmax1 = 0.d0
      !do ik = 1, nkpt
      !do ist = 1, nstsv
      !    if (e(ist,ik).ge.emin.and.e(ist,ik).le.emax) wmax1 = max(wmax1,w(ist,ik))
      !enddo
      !enddo
      !write(*,*)'Info: ratio of character weights:',wmax1/wmax
      
      scale = 1.d0
      write(*,*)'Input maximum line width [eV]'
      read(*,*)scale
      
      scale = 0.5d0*scale/wmax
      
      open(50,file='BNDS.DAT',form='formatted',status='replace')
      open(51,file='BNDS1.DAT',form='formatted',status='replace')
      open(52,file='BNDS2.DAT',form='formatted',status='replace')
      open(53,file='BNDS3.DAT',form='formatted',status='replace')
      do ist = 1, nstsv
        do ik = 1, nkpt
          write(50,*)dpp1d(ik),e(ist,ik)
          write(51,*)dpp1d(ik),e(ist,ik)+scale*w(ist,ik)
          write(52,*)dpp1d(ik),e(ist,ik)-scale*w(ist,ik)
          write(53,*)dpp1d(ik),e(ist,ik)+scale*w(ist,ik)
          write(53,*)dpp1d(ik),e(ist,ik)-scale*w(ist,ik)
          write(53,*)
        enddo
        write(50,*)
        write(51,*)
        write(52,*)
      enddo
      close(50) 
      close(51) 
      close(52) 
      close(53) 
      
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
      write(50,*)0.d0,0.d0
      write(50,*)dpp1d(nkpt),0.d0
      close(50)
     
      
      end
      