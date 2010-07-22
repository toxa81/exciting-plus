subroutine genchi0(iq)
use modmain
use mod_addons_q
use mod_nrkp
implicit none
! arguments
integer, intent(in) :: iq
! local variables
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0wan_k(:,:,:)
complex(8), allocatable :: chi0wan(:,:)
complex(8), allocatable :: mexp(:,:,:)
!complex(8), allocatable :: megqwan1(:,:)
complex(8), allocatable :: megqwan_tmp(:,:)
integer, allocatable :: imegqwan_tmp(:,:)
character*10 c1,c2,c3,c4
character, parameter :: orb(4)=(/'s','p','d','f'/)
integer lm1,lm2,ias,jas

integer i,iw,i1,i2,ikloc,n,j
integer ist1,ist2,nfxcloc,nwloc,jwloc,iwloc
integer it1(3),it2(3),it(3)
character*100 qnm,qdir,fout,fstat
integer ik,n1,n2
real(8) vtc(3)

real(8) t1,t2,t3,t4,t5,t6,t7

call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI0.OUT"
  open(150,file=trim(fout),form="FORMATTED",status="REPLACE")
  fstat=trim(qnm)//"_chi0_stat.txt"
endif

if (wproc) then
  write(150,*)
  write(150,'("Calculation of chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  write(150,'("  energy interval [eV] : ", 2F9.4)')lr_w0,lr_w1
  write(150,'("  energy step     [eV] : ", F9.4)')lr_dw
  write(150,'("  eta             [eV] : ", F9.4)')lr_eta
  write(150,*)  
  write(150,'("Included band interval (Ha)        : ",2F8.2)')&
    chi0_include_bands(1),chi0_include_bands(2)
  write(150,'("Excluded band interval (Ha)        : ",2F8.2)')&
    chi0_exclude_bands(1),chi0_exclude_bands(2) 
  call flushifc(150)
endif
  
! cutoff small matrix elements
if (wannier_chi0_chi) then
  allocate(megqwan_tmp(nmegqwan,ngvecme))
  megqwan_tmp=zzero
  allocate(imegqwan_tmp(5,nmegqwanmax))
  imegqwan_tmp=0
  j=0
  do i=1,nmegqwan
    if (abs(megqwan(i,iig0q)).ge.megqwan_cutoff1.and.&
        abs(megqwan(i,iig0q)).le.megqwan_cutoff2) then
      j=j+1
      imegqwan_tmp(:,j)=imegqwan(:,i)
      megqwan_tmp(j,:)=megqwan(i,:)
    endif
  enddo
  nmegqwan=j
  megqwan(:,:)=megqwan_tmp(:,:)
  imegqwan(:,:)=imegqwan_tmp(:,:)
  deallocate(megqwan_tmp)
  deallocate(imegqwan_tmp)
endif

! for response in Wannier bais
if (wannier_chi0_chi) then
  if (wproc) then
    write(150,*)
    write(150,'("Wannier chi0 AFM : ",L1)')wannier_chi0_afm
    write(150,*)
    write(150,'("megqwan value cutoff (min,max) : ",2F12.6)')&
      megqwan_cutoff1,megqwan_cutoff2
    write(150,'("megqwan distance cutoff (min,max) : ",2F12.6)')&
      megqwan_mindist,megqwan_maxdist
    write(150,'("Number of Wannier transitions after cutoff : ",I6)')nmegqwan
    write(150,'("List of Wannier transitions")')
    do i=1,nmegqwan
      ias=iwann(1,imegqwan(1,i))
      lm1=iwann(2,imegqwan(1,i))
      jas=iwann(1,imegqwan(2,i))
      lm2=iwann(2,imegqwan(2,i))
      vtc(:)=avec(:,1)*imegqwan(3,i)+avec(:,2)*imegqwan(4,i)+&
             avec(:,3)*imegqwan(5,i)  
      vtc(:)=vtc(:)+atposc(:,ias2ia(jas),ias2is(jas))-&
        atposc(:,ias2ia(ias),ias2is(ias))       
      write(c1,'(I6)')ias2ia(ias)
      write(c2,'(I1)')lm2m(lm1)+lm2l(lm1)+1
      c3="("//trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
        orb(lm2l(lm1)+1)//trim(adjustl(c2))//")"
      write(c1,'(I6)')ias2ia(jas)
      write(c2,'(I1)')lm2m(lm2)+lm2l(lm2)+1
      c4="("//trim(spsymb(ias2is(jas)))//trim(adjustl(c1))//"-"//&
        orb(lm2l(lm2)+1)//trim(adjustl(c2))//")"
      write(150,'("i : ",I4,"   ",I4," ",A,"  -> ",I4," ",A,"    R=",3F12.6,&
        &"   D=",F12.6,"  |me(G0q)|=",G18.10)')i,imegqwan(1,i),&
        trim(c3),imegqwan(2,i),trim(c4),vtc,sqrt(sum(vtc(:)**2)),&
        abs(megqwan(i,iig0q))
    enddo
    call flushifc(150)
  endif
  allocate(chi0wan(nmegqwan,nmegqwan))
  allocate(chi0wan_k(nmegqwan,nmegqwan,nkptnrloc))
  allocate(mexp(nmegqwan,nmegqwan,nkptnrloc))
  do i1=1,nmegqwan
    do i2=1,nmegqwan
      it1(:)=imegqwan(3:5,i1)
      it2(:)=imegqwan(3:5,i2)
      it(:)=it1(:)-it2(:)
      vtc(:)=avec(:,1)*it(1)+avec(:,2)*it(2)+avec(:,3)*it(3)
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! phase e^{i(k+q)T}
        mexp(i1,i2,ikloc)=exp(dcmplx(0.d0,dot_product(vkcnr(:,ik)+vqc(:,iq),&
          vtc(:))))   
      enddo
    enddo
  enddo
! arrangement for zgemm  
  allocate(wann_cc(nmegqblhwanmax,nmegqwan,nkptnrloc))
  allocate(wann_cc2(nmegqblhwanmax,nmegqwan))
  wann_cc=zzero
  do ikloc=1,nkptnrloc
    do i1=1,nmegqblhwan(ikloc)
      i=imegqblhwan(i1,ikloc)
      ist1=bmegqblh(1,i,ikloc)
      ist2=bmegqblh(2,i,ikloc)
      do n=1,nmegqwan
        n1=imegqwan(1,n)
        n2=imegqwan(2,n)
        wann_cc(i1,n,ikloc)=wanncnrloc(n1,ist1,ikloc)*dconjg(wann_c_jk(n2,ist2,ikloc))
      enddo
    enddo !i1
  enddo !ikloc
endif !wannier_chi0_chi

allocate(chi0(ngvecme,ngvecme))
allocate(megqblh2(nmegqblhlocmax,ngvecme))

! distribute nfxca between 3-rd dimension 
nfxcloc=mpi_grid_map(nfxca,dim_b)
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)

allocate(chi0loc(ngvecme,ngvecme,nwloc))
if (wannier_chi0_chi) allocate(chi0wanloc(nmegqwan,nmegqwan,nwloc))

call timer_start(1,reset=.true.)
call timer_reset(2)
call timer_reset(3)
call timer_reset(4)
call timer_reset(5)
call timer_reset(6)
call timer_reset(7)
call timer_reset(8)
! loop over energy points
do iw=1,lr_nw
  chi0=zzero
  if (wannier_chi0_chi) chi0wan_k=zzero
! sum over fraction of k-points
  call timer_start(2)
  do ikloc=1,nkptnrloc
    if (nmegqblhloc(1,ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call genchi0blh(ikloc,lr_w(iw),chi0)
    endif
  enddo
! find the processor j which will get the full chi0 and chi0wan matrices
  jwloc=mpi_grid_map(lr_nw,dim_k,glob=iw,x=j)
! sum over k-points and band transitions
  call mpi_grid_reduce(chi0(1,1),ngvecme*ngvecme,dims=(/dim_k,dim_b/),&
    root=(/j,0/))
  chi0=chi0/nkptnr/omega
! processor j saves chi0 to local array  
  if (mpi_grid_x(dim_k).eq.j) chi0loc(:,:,jwloc)=chi0(:,:)
  call timer_stop(2)
! for response in Wannier basis
  if (wannier_chi0_chi) then
    call timer_start(3)
    chi0wan(:,:)=zzero
    do ikloc=1,nkptnrloc
      if (nmegqblhwan(ikloc).gt.0) then
        call genchi0wan_k(ikloc,lr_w(iw),chi0wan_k(1,1,ikloc))
      endif
      chi0wan(:,:)=chi0wan(:,:)+mexp(:,:,ikloc)*chi0wan_k(:,:,ikloc)
    enddo !ikloc
! sum chi0wan over all k-points
    call mpi_grid_reduce(chi0wan(1,1),nmegqwan*nmegqwan,dims=(/dim_k/),&
      root=(/j/))
    chi0wan(:,:)=chi0wan(:,:)/nkptnr/omega
    if (wannier_chi0_afm) chi0wan(:,:)=chi0wan(:,:)*2.d0
! processor j saves chi0wan to local array  
    if (mpi_grid_x(dim_k).eq.j) chi0wanloc(:,:,jwloc)=chi0wan(:,:)
!    if (iw.eq.120.and..true.) then
!      call chi0wan_diag(chi0wan)
!    endif
    call timer_stop(3)
  endif !wannier_chi0_chi
  if (wproc) then
    open(160,file=trim(fstat),status="REPLACE",form="FORMATTED")
    write(160,'(I8)')iw
    close(160)
  endif
enddo !iw
call timer_stop(1)
t1=timer_get_value(1)
t2=timer_get_value(2)
t3=timer_get_value(3)
if (wproc) then
  write(150,*)
  write(150,'("Total time per frequency point   : ",F8.2)')t1/lr_nw
  write(150,'("  Bloch basis part (chi0)        : ",F8.2)')t2/lr_nw
  write(150,'("  Wannier basis part (chi0)      : ",F8.2)')t3/lr_nw
  call flushifc(150)
endif

call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(chi0)
deallocate(megqblh2)
if (wannier_chi0_chi) then
  deallocate(chi0wan)
  deallocate(chi0wan_k)
  deallocate(mexp)
  deallocate(wann_cc)
  deallocate(wann_cc2)
endif
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif
return
end
