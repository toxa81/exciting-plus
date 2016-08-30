subroutine writewann
use modmain
use modldapu
use mod_nrkp
implicit none
integer i,j,ik,ikloc,n1,n2,j1,j2,nwan,ias,n,is,ia
logical lpmat
integer values(8), seconds, hash, isp, hdim, ntype, itype, istart, l, norb, lnext
real(8)  reh, imh, sc, sc1, sc2, sc3, alat(3,3)
character*1 xx
character(256) block
complex(8), allocatable :: zm(:,:)
real(8), allocatable :: dm(:,:),eval(:)
integer, allocatable   :: basis_desc(:,:)
real*8,            parameter :: zero=0.d0
call init0
call init1

wproc=mpi_grid_root()
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(writewann_h) : WF generation is switched off")')
  write(*,*)
  call pstop
endif
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
lpmat=.false.
if (task.eq.808) lpmat=.true.
call genwfnr(6,.false.)
if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(nwantot,nwantot,nkptnr))
wann_h=zzero
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(nwantot,nkptnr))
wann_e=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call genwann_h(.true.,evalsvnr(1,ik),wanncnrloc(1,1,ikloc),&
    &wann_h(1,1,ik),wann_e(1,ik))
enddo

!allocate(wann_ene_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!allocate(wann_occ_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!call wann_ene_occ_(wann_ene_m,wann_occ_m)

call mpi_grid_reduce(wann_h(1,1,1),nwantot*nwantot*nkptnr,dims=(/dim_k/),side=.true.)
call mpi_grid_reduce(wann_e(1,1),nwantot*nkptnr,dims=(/dim_k/),side=.true.)
!call mpi_grid_reduce(wann_p(1,1,1,1),3*nwantot*nwantot*nkpt,dims=(/dim_k/),side=.true.)
if (mpi_grid_root().and.task.eq.807) then
  call readfermi
  open(200,file="WANN_H.OUT",form="FORMATTED",status="REPLACE")
  write(200,'("# units of energy are Hartree, 1 Ha=",F18.10," eV")')ha2ev
  write(200,'("# fermi energy")')
  write(200,'(G18.10)')efermi
  write(200,'("# lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')avec(:,i)
  enddo
  write(200,'("# reciprocal lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')bvec(:,i)
  enddo
  write(200,'("# k-grid size")')
  write(200,'(3I6)')ngridk
  write(200,'("# number of k-points")')
  write(200,'(I8)')nkptnr
  write(200,'("# number of Wannier functions")')
  write(200,'(I8)')nwantot
  write(200,'("# number of atoms")')
  write(200,'(I8)')natmtot 
  write(200,'("# number of species")')
  write(200,'(I8)')nspecies  
  write(200,'("# wf -> atom mapping")')
  do n=1,nwantot
    write(200,'(2I8)')n,wan_info(wi_atom,n)
  enddo
  write(200,'("# atom -> species mapping")')
  do ias=1,natmtot
    is=ias2is(ias)
    write(200,'(2I8)')ias,is
  enddo
  do ias=1,natmtot
    ia=ias2ia(ias)
    is=ias2is(ias)
    write(200,'("# atom : ",I8)')ias
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')atposc(:,ia,is)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')atposl(:,ia,is)
  enddo
!  write(200,'("# occupancy matrix")')
!  do i=1,wann_natom
!    ias=wann_iprj(1,i)
!    write(200,'("#   atom : ",I4)')ias
!    do l=0,lmaxlu
!      if (sum(abs(wann_occ_m(idxlm(l,-l):idxlm(l,l),idxlm(l,-l):idxlm(l,l),:,:,ias))).gt.1d-8) then
!        t=0.0
!        do ispn=1,nspinor
!          write(200,'("#     ispn : ",I1)')ispn
!          write(200,'("#     real part")')
!          do lm1=l**2+1,(l+1)**2
!            write(200,'("#",1X,7F12.6)')(dreal(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
!            t(ispn)=t(ispn)+dreal(wann_occ_m(lm1,lm1,ispn,ispn,ias))
!          enddo
!          write(200,'("#     imag part")')
!          do lm1=l**2+1,(l+1)**2
!            write(200,'("#",1X,7F12.6)')(dimag(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
!          enddo
!          write(200,'("#     spin occupancy : ",F12.6)')t(ispn)
!        enddo !ispn
!        write(200,'("#   total occupancy : ",F12.6)')sum(t)
!      endif
!    enddo
!  enddo
  do ik=1,nkptnr
    write(200,'("# k-point : ",I8)')ik
    write(200,'("# weight")')
    write(200,'(G18.10)')wkptnr(ik)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')vklnr(:,ik)
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')vkcnr(:,ik)
    write(200,'("# real part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dreal(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# imaginary part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dimag(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# eigen-values of H")')
    write(200,'(255G18.10)')(wann_e(j,ik),j=1,nwantot)
  enddo
  close(200)
endif

if (mpi_grid_root().and.task.eq.809) then
  call readfermi

  call date_and_time(values=values)
  ! Inaccurate but easy calculation of timestamp.
  hash = (values(1) - 1970) * 365 * 24 * 60 * 60  &  
     + (values(2)-1) * 31 * 24 * 60 * 60  &
     + values(3 + 11 ) * 24 * 60 * 60  & 
     + (values(5)-values(4)) * 60 * 60  &
     + values(6) * 60  + values(7)  &    
     + int(values(8)/1000)               

  open(200, file="hamilt.am", form="FORMATTED", status="REPLACE")

  write(200,'("# This file was written on: ",I2.2,".",I2.2,".",I4.4,"  &
                                           ",I2.2,":",I2.2,":",I2.2)') &
             values(3),values(2),values(1),values(5),values(6),values(7)
  write(200,*)

  write(200,'(a5)') '&hash'
  write(200,*) hash
  write(200,*)
  
  write(200,'(a10)') '&Codestamp'
  write(200,*) 'exciting-plus'
  write(200,*)
  
  
  write(200,'(a6)') '&nspin'
  write(200,*) nspinor
  write(200,*)
  
  write(200,'(a4)') '&nkp'
  write(200,*) nkptnr
  write(200,*)

  hdim = int(nwantot/nspinor)
  write(200,'(a4)') '&dim'
  write(200,*) hdim !nwantot
  write(200,*)

  open(50,file='TOTENERGY.OUT',action='READ',form='FORMATTED',status='OLD')
  do while (.true.)
    read(50, '(G22.12)', iostat=i) engytot
    if (i /= 0) exit
  end do
  close(50)
  write(200,'(a5)') '&etot'
  write(200,*) engytot * ha2ev
  write(200,*)

  write(200,'(a7)') '&fermi'
  write(200,*) efermi * ha2ev
  write(200,*)


  write(200,'(a11)') '&crystcoord'
  write(200,*) 'true'
  write(200,*)

  write(200,'(a8)') '&kpoints'
  do ik=1,nkptnr
!   write(200,'(f15.12,3f9.5)') wkptnr(ik), vkcnr(:,ik)
    write(200,'(f15.12,3f9.5)') wkptnr(ik), vklnr(:,ik)
  end do
  write(200,*) 
  
  hdim = int(nwantot/nspinor)
  write(200,'(a12)') '&hamiltonian'
  do isp=1,nspinor 
	do ik = 1, nkptnr
	  do i = 1, nwantot
		do j = i, nwantot
		  if (wan_info(wi_spin,i) /= isp .or. wan_info(wi_spin,j) /= isp) cycle
		  reh=dreal(wann_h(i,j,ik)) * ha2ev
		  imh=dimag(wann_h(i,j,ik)) * ha2ev
		  if(dabs(reh)<1.d-12 .and. dabs(imh)<1.d-12) then
			write(200,'(6x,f2.0,18x,f2.0)')zero,zero
		  elseif(dabs(reh)<1.d-12) then
			write(200,'(6x,f2.0,12x,f20.12)')zero,imh
		  elseif(dabs(imh)<1.d-12) then
			write(200,'(f20.12,6x,f2.0)')reh,zero
		  else
			write(200,'(2f20.12)')reh,imh
		  end if
		end do !n2
	  end do   !n1
	end do     !ikp
  end do       !isp

  close(200)
  
  open(200, file="system.am", form="FORMATTED", status="REPLACE")

  write(200,'("# This file was written on: ",I2.2,".",I2.2,".",I4.4,"  &
                                           ",I2.2,":",I2.2,":",I2.2)') &
             values(3),values(2),values(1),values(5),values(6),values(7)
  write(200,*)

  write(200,'(a5)') '&hash'
  write(200,*) hash
  write(200,*)
  
  write(200,'(a10)') '&Codestamp'
  write(200,*) 'exciting-plus'
  write(200,*)

!print*,alat  
  sc  = 1.0
  sc1 = 1.0
  sc2 = 1.0
  sc3 = 1.0
  open(50,file='elk.in',action='READ',status='OLD',form='FORMATTED')
10 continue
  read(50,*,end=30) block
  if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
  select case(trim(block))
  case('avec')
    read(50,*,err=20) alat(:,1)
    read(50,*,err=20) alat(:,2)
    read(50,*,err=20) alat(:,3)
  case('scale')
    read(50,*,err=20) sc
  case('scale1')
    read(50,*,err=20) sc1
  case('scale2')
    read(50,*,err=20) sc2
  case('scale3')
    read(50,*,err=20) sc3
  case('')
    goto 10
  end select
  goto 10
20 continue
  write(*,*)
  write(*,'("Error(readinput): error reading from elk.in")')
  write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
  write(*,'("Check input convention in manual")')
  write(*,*)
  stop
30 continue
  close(50)

  alat(:,1)=sc1*alat(:,1)
  alat(:,2)=sc2*alat(:,2)
  alat(:,3)=sc3*alat(:,3)

  write(200,'(a5)') '&cell'
  write(200,'(f12.9)') sc !1.0 
  do i=1,3
    write(200,'(3f9.5)') alat(:,i) !avec(:,i) !/tmp 
  end do
  write(200,*)

  write(200,'(a7)') '&fermi'
  write(200,*) efermi * ha2ev
  write(200,*)

  write(200,'(a11)') '&crystcoord'
  write(200,*) 'true'
  write(200,*)

  write(200,'(a6)') '&atoms'
  write(200,*) natmtot
  do ias=1,natmtot
    is=ias2is(ias)
    ia=ias2ia(ias)
!   write(200,'(a4,x,3f9.5)') trim(spsymb(is)), atposc(:,ia,is)/sc
    write(200,'(a4,x,3f9.5)') trim(spsymb(is)), atposl(:,ia,is)
  end do
!  write(200,*) natmtot
!  do ias=1,natmtot
!    is=ias2is(ias)
!    ia=ias2ia(ias)
!    write(200,'(a4,x,3f9.5)') trim(spsymb(is)), atposc(:,ia,is)
!  end do
  write(200,*)

  ntype = int(wann_natom/nspinor)
  allocate (basis_desc(11,ntype))
  itype = 1
  istart = 1
  do i = 1,hdim-1
    basis_desc(1,itype) = wan_info(wi_atom,i)
    l = wan_info(wi_lm,i)
    if (l < 2)                 basis_desc(2,itype) = 0
    if (l < 5 .and. l >= 2 )   basis_desc(2,itype) = 1
    if (l < 10 .and. l >= 5 )  basis_desc(2,itype) = 2
    if (l < 17 .and. l >= 10 ) basis_desc(2,itype) = 3
    basis_desc(3,itype) = istart
    l = basis_desc(2,itype)
    if (basis_desc(2,itype) == 0) norb = 1
    if (basis_desc(2,itype) == 1) norb = 3
    if (basis_desc(2,itype) == 2) norb = 5
    if (basis_desc(2,itype) == 3) norb = 7
    do j=0,norb-1
      basis_desc(4+j,itype) = l*l+1+j
    end do

    lnext = wan_info(wi_lm,i+1)
    if (lnext < 2)                     lnext = 0
    if (lnext < 5 .and.  lnext >= 2 )  lnext = 1
    if (lnext < 10 .and. lnext >= 5 )  lnext = 2
    if (lnext < 17 .and. lnext >= 10 ) lnext = 3
    
    if ((wan_info(wi_atom,i) .ne. wan_info(wi_atom,i+1)) .or. &
          l .ne. lnext) then 
      itype = itype +1
      istart = istart + norb
    end if
  end do

  write(200,'(a20)') '# Basis description:'
  write(200,'(a14)') '# dim, nblocks'
  write(200,'(a74)') '# atom_sym, atom_num, l_sym, block_dim, block_start, orbitals(1:block_dim)'
  write(200,'(a6)') '&basis'
  write(200,'(i2,i4)') hdim, ntype
  do i = 1,ntype
    select case(basis_desc(2,i))
    case(0)
      xx='s'
      norb = 1
    case(1)
      xx='p'
      norb = 3
    case(2)
      xx='d'
      norb = 5
    case(3)
      xx='f'
      norb = 7
    end select
print'(11i3)', basis_desc(1,itype), basis_desc(2,itype), basis_desc(3,itype), basis_desc(4,itype), basis_desc(5,itype), basis_desc(6,itype),basis_desc(7,itype), basis_desc(8,itype), basis_desc(9,itype),basis_desc(10,itype), basis_desc(11,itype)
    is=ias2is(basis_desc(1,i))
    write(200,'(a3,i3,a2,i2,i4,4x,16i2)')trim(spsymb(is)), basis_desc(1,i), &
               xx, norb, basis_desc(3,i), (basis_desc(4+j,i), j=0,norb-1)
  end do

 write(200,*)

  deallocate(basis_desc)
  close(200)

endif

if (mpi_grid_root()) then
  do ias=1,natmtot
    nwan=nwannias(ias)
    if (nwan.ne.0) then
      write(*,*)"ias : ",ias,"  nwan : ",nwan
      allocate(zm(nwan,nwan))
      allocate(dm(nwan,nwan))
      allocate(eval(nwan))
      j1=0
      do n1=1,nwantot
        if (wan_info(wi_atom,n1).eq.ias) then
          j1=j1+1
          j2=0
          do n2=1,nwantot
            if (wan_info(wi_atom,n2).eq.ias) then
              j2=j2+1
              zm(j1,j2)=wann_h(n1,n2,1)
              dm(j1,j2)=dreal(zm(j1,j2))
            endif
          enddo
        endif
      enddo
      write(*,*)"Hamiltonian at gamma point :"
      do j1=1,nwan
        write(*,'(255F12.6)')(dreal(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)
      do j1=1,nwan
        write(*,'(255F12.6)')(dimag(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)"eigen-vectors : "
      call diagdsy(nwan,dm,eval)
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j1,j2),j2=1,nwan)
      enddo
      write(*,*)
      write(*,'(2X,7G18.10)')(eval(j1),j1=1,nwan)
      write(*,*)
      write(*,*)"transpose of eigen-vectors : "
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j2,j1),j2=1,nwan)
      enddo
      write(*,*)
      deallocate(zm,dm,eval)
    endif
  enddo
  write(*,*)
  do ias=1,natmtot
    nwan=nwannias(ias)
    if (nwan.ne.0) then
      write(*,*)"ias : ",ias,"  nwan : ",nwan
      allocate(zm(nwan,nwan))
      zm=zzero
      allocate(dm(nwan,nwan))
      allocate(eval(nwan))
      j1=0
      do n1=1,nwantot
        if (wan_info(wi_atom,n1).eq.ias) then
          j1=j1+1
          j2=0
          do n2=1,nwantot
            if (wan_info(wi_atom,n2).eq.ias) then
              j2=j2+1
              do ik=1,nkptnr
                zm(j1,j2)=zm(j1,j2) + wann_h(n1,n2,ik)/nkptnr
              enddo
              dm(j1,j2)=dreal(zm(j1,j2))
            endif
          enddo
        endif
      enddo
      write(*,*)"Average Hamiltonian :"
      do j1=1,nwan
        write(*,'(255F12.6)')(dreal(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)
      do j1=1,nwan
        write(*,'(255F12.6)')(dimag(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)"eigen-vectors : "
      call diagdsy(nwan,dm,eval)
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j1,j2),j2=1,nwan)
      enddo
      write(*,*)
      write(*,'(2X,7G18.10)')(eval(j1),j1=1,nwan)
      write(*,*)
      write(*,*)"transpose of eigen-vectors : "
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j2,j1),j2=1,nwan)
      enddo
      write(*,*)
      deallocate(zm,dm,eval)
    endif
  enddo
endif

!deallocate(wann_ene_m,wann_occ_m)
return
end
