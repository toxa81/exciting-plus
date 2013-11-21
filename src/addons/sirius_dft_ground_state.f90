#ifdef _SIRIUS_
subroutine sirius_dft_ground_state
use modmain
use modldapu
implicit none
integer is,ia,ik,i,mi,n,nwork,ias,rank,ispn,jspn,l,m1,m2,lm1,lm2,lm,kset_id,kset_fbz_id,i1,i2,i3
real(8), allocatable :: work(:),v(:)
real(8) dv,d1,d2,etot,etot_old


call sirius_elk_init(kset_id)

! WARNING: symmetrization routines for the interstitial part were changed to compute complex phase factors for translation
!          in a different way
!   old code:
!     t1=-dot_product(vgc(:,ig),vtc(:))
!     zt1=cmplx(cos(t1),sin(t1),8)
!   new code:
!     zt1=exp(dcmplx(0.d0,-twopi*(ivg(1,ig)*vtlsymc(1,isym)+ivg(2,ig)*vtlsymc(2,isym)+ivg(3,ig)*vtlsymc(3,isym))))
!  
!   The reason for this is to get rid of the dependency on vgc array (Cartesian coordinates of G-vectors)
!
! WARNING: output of core eigen-values is commented
!
! WARNING: varible ncmag is repalced by ndmag.eq.3 in the symmetrization routines
!
! Test of the ground state calculation
!

call sirius_ground_state_initialize(kset_id)

! generate initial density
call sirius_initial_density
call sirius_generate_effective_potential

! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngrtot
if (spinpol) n=n*(1+ndmag)
if (ldapu.ne.0) n=n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
allocate(v(n))
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,n,v,dv,nwork,v)
allocate(work(nwork))

call sirius_platform_mpi_rank(rank)
etot_old = 0d0
etot = 1d100

do iscl=1,maxscl
  call sirius_start_timer("elk::iteration")

!---------------!
! mix potential !
!---------------!
  call sirius_start_timer("elk::mixer")
! pack interstitial and muffin-tin effective potential and field into one array
  call mixpack(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)
! unpack potential and field
  call mixpack(.false.,n,v)
  call sirius_stop_timer("elk::mixer")

! check if all ranks can exit the scf loop
  i = 0   
  if ((dv.lt.epspot).and.abs(etot-etot_old).lt.epsengy) i = 1
  call sirius_global_set_sync_flag(i)
  call sirius_global_get_sync_flag(i)
  if (i.eq.1) exit

! solve H\psi = E\psi
  call sirius_start_timer("elk::find_eigen_states")
  call sirius_find_eigen_states(kset_id, 1)
  call sirius_stop_timer("elk::find_eigen_states")

  call sirius_start_timer("elk::band_occ")
! find Fermi level and band occupancies 
  call sirius_find_band_occupancies(kset_id)

  do ik=1,nkpt
    call sirius_get_band_energies(kset_id, ik, evalsv(1,ik))
    call sirius_get_band_occupancies(kset_id, ik, occsv(1,ik))
  enddo
!  call occupy

  if (rank.eq.0) call writeeval
  call sirius_stop_timer("elk::band_occ")

!  do ik=1,nkpt 
!    call sirius_density_set_band_occupancies(ik,occsv(1,ik))
!  enddo
  
! generate density and magnetization
  call sirius_start_timer("elk::generate_density")
  call sirius_generate_density(kset_id)
  call sirius_stop_timer("elk::generate_density")

  if (ldapu.ne.0) then
    dmatlu=zzero
    do ias=1,natmtot
      call sirius_get_occupation_matrix(ias, dmatlu(1,1,1,1,ias))
    enddo
! generate the LDA+U potential matrix
    call genvmatlu
! write the LDA+U matrices to file
    call writeldapu
    do is=1,nspecies
      if (llu(is).ge.0) then
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          call sirius_set_uj_correction_matrix(ias, llu(is), vmatlu(1,1,1,1,ias)) 
        enddo
      endif
    enddo
  endif

  call sirius_start_timer("elk::symmetrization")
  call symrf(1,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(1,magmt,magir)
  call sirius_stop_timer("elk::symmetrization")
 
! check the total charge
  !!== call sirius_start_timer("elk::integrate_density")
  !!== call sirius_density_integrate
  !!== call sirius_stop_timer("elk::integrate_density")

! generate potential
  call sirius_start_timer("elk::potential")
  call sirius_generate_effective_potential
  call sirius_stop_timer("elk::potential")

  call sirius_start_timer("elk::symmetrization")
! symmetrize potential
  call symrf(1,veffmt,veffir)
! symmetrise magnetic field
  if (spinpol) call symrvf(1,bxcmt,bxcir)
  call sirius_stop_timer("elk::symmetrization")
  
  call sirius_ground_state_print_info
  etot_old = etot
  call sirius_get_total_energy(etot)

  write(*,'("charge RMS : ",G18.10," total energy difference : ",G18.10)')dv,abs(etot - etot_old)
  call sirius_platform_barrier()
 
  call sirius_stop_timer("elk::iteration")  
enddo

!!== nkptnr = ngridk(1) * ngridk(2) * ngridk(3)
!!== if (allocated(vklnr)) deallocate(vklnr)
!!== allocate(vklnr(3, nkptnr))
!!== if (allocated(wkptnr)) deallocate(wkptnr)
!!== allocate(wkptnr(nkptnr))
!!== 
!!== ik = 0
!!== do i1 = 0, ngridk(1) - 1
!!==   do i2 = 0, ngridk(2) - 1
!!==     do i3 = 0, ngridk(3) - 1
!!==       ik = ik + 1
!!==       vklnr(:, ik) = (/ dble(i1) / ngridk(1), dble(i2) / ngridk(2), dble(i3) / ngridk(3) /)
!!==       wkptnr = 1.d0 / nkptnr
!!==     enddo
!!==   enddo
!!== enddo
!!== 
!!== call sirius_create_kset(nkptnr, vklnr, wkptnr, 1, kset_fbz_id)
!!== call sirius_find_eigen_states(kset_fbz_id, 1)
!!== call sirius_find_band_occupancies(kset_fbz_id)


call sirius_create_storage_file
!!== call sirius_save_kset(kset_fbz_id)
call sirius_save_potential
call sirius_print_timers

call sirius_write_json_output
call sirius_ground_state_clear
call sirius_platform_barrier
call sirius_delete_kset(kset_id)
!!== call sirius_delete_kset(kset_fbz_id)
call sirius_clear

return
end subroutine

!!== subroutine test_sirius_band
!!== use modmain
!!== implicit none
!!== integer is,ia,ik,i,mi,n,nwork,ist
!!== real(8), allocatable :: work(:),v(:)
!!== real(8) dv
!!== 
!!== call sirius_elk_init
!!== 
!!== !
!!== ! WARNING: symmetrization routines for the interstitial part were changed to compute complex phase factors for translation
!!== !          in a different way
!!== !   old code:
!!== !     t1=-dot_product(vgc(:,ig),vtc(:))
!!== !     zt1=cmplx(cos(t1),sin(t1),8)
!!== !   new code:
!!== !     zt1=exp(dcmplx(0.d0,-twopi*(ivg(1,ig)*vtlsymc(1,isym)+ivg(2,ig)*vtlsymc(2,isym)+ivg(3,ig)*vtlsymc(3,isym))))
!!== !  
!!== !   The reason for this is to get rid of the dependency on vgc array (Cartesian coordinates of G-vectors)
!!== !
!!== ! WARNING: output of core eigen-values is commented
!!== !
!!== ! WARNING: varible ncmag is repalced by ndmag.eq.3 in the symmetrization routines
!!== !
!!== ! Test of the ground state calculation
!!== !
!!== 
!!== ! set pointers to effective potential
!!== call sirius_set_effective_potential_ptr(veffmt,veffir)
!!== if (spinpol) then
!!== ! set pointer to effective magnetic field
!!==   call sirius_set_effective_magnetic_field_ptr(bxcmt,bxcir)
!!== endif
!!== 
!!== call sirius_bands()
!!== 
!!== call sirius_print_info()
!!== call sirius_print_timers()
!!== call sirius_platform_barrier()
!!== call sirius_clear()
!!== 
!!== return
!!== end subroutine


subroutine sirius_elk_init(kset_id)
use modmain
use modldapu
implicit none
integer, intent(out) :: kset_id
integer ist,l,m,lm,ncls,ia1,ia2
integer, allocatable :: icls(:)
real(8) sum
logical lsym(48)
integer isym,is,ia,ias
integer ik,io,ilo,iv(3)
integer i1,i2,i3,ispn
integer i
real(8) vl(3),vc(3)
real(8) boxl(3,4),t1

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!
lmmaxvr=(lmaxvr+1)**2
lmmaxapw=(lmaxapw+1)**2
! index to (l,m) pairs
if (allocated(idxlm)) deallocate(idxlm)
allocate(idxlm(0:50,-50:50))
lm=0
do l=0,50
  do m=-l,l
    lm=lm+1
    idxlm(l,m)=lm
  end do
end do
! array of i**l values
if (allocated(zil)) deallocate(zil)
allocate(zil(0:50))
do l=0,50
  zil(l)=zi**l
end do

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
! check if the system is an isolated molecule
if (molecule) then
  primcell=.false.
  tshift=.false.
end if
! find primitive cell if required
if (primcell) call findprim
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias

!------------------------!
!     spin variables     !
!------------------------!
! spin-orbit coupling or fixed spin moment implies spin-polarised calculation
if ((spinorb).or.(fixspin.ne.0).or.(spinsprl)) spinpol=.true.
! number of spinor components and maximum allowed occupancy
if (spinpol) then
  nspinor=2
  occmax=1.d0
else
  nspinor=1
  occmax=2.d0
end if
! number of spin-dependent first-variational functions per state
if (spinsprl) then
  nspnfv=2
else
  nspnfv=1
end if
! spin-polarised calculations require second-variational eigenvectors
if (spinpol) tevecsv=.true.
! set the magnetic fields to the initial values
bfieldc(:)=bfieldc0(:)
bfcmt(:,:,:)=bfcmt0(:,:,:)
! check for collinearity in the z-direction and set the dimension of the
! magnetisation and exchange-correlation vector fields
if (spinpol) then
  ndmag=1
  if ((abs(bfieldc(1)).gt.epslat).or.(abs(bfieldc(2)).gt.epslat)) ndmag=3
  do is=1,nspecies
    do ia=1,natoms(is)
      if ((abs(bfcmt(1,ia,is)).gt.epslat).or.(abs(bfcmt(2,ia,is)).gt.epslat)) ndmag=3
    end do
  end do
! source-free fields and spin-spirals are non-collinear in general
  if ((nosource).or.(spinsprl)) ndmag=3
! spin-orbit coupling is non-collinear in general
  if (spinorb) ndmag=3
else
  ndmag=0
end if
! spin-polarised cores
if (.not.spinpol) spincore=.false.
! set fixed spin moment effective field to zero
bfsmc(:)=0.d0
! set muffin-tin FSM fields to zero
bfsmcmt(:,:,:)=0.d0

!----------------------------------!
!     crystal structure set up     !
!----------------------------------!
! generate the reciprocal lattice vectors and unit cell volume
call reciplat
! compute the inverse of the lattice vector matrix
call r3minv(avec,ainv)
! compute the inverse of the reciprocal vector matrix
call r3minv(bvec,binv)
! Cartesian coordinates of the spin-spiral vector
call r3mv(bvec,vqlss,vqcss)
do is=1,nspecies
  do ia=1,natoms(is)
! map atomic lattice coordinates to [0,1) if not in molecule mode
    if (.not.molecule) call r3frac(epslat,atposl(:,ia,is),iv)
! determine atomic Cartesian coordinates
    call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
  end do
end do

!-------------------------------!
!     vector fields E and A     !
!-------------------------------!
efieldpol=.false.
if ((abs(efieldc(1)).gt.epslat).or.&
   &(abs(efieldc(2)).gt.epslat).or.&
   &(abs(efieldc(3)).gt.epslat)) then
  efieldpol=.true.
  tshift=.false.
! electric field vector in lattice coordinates
  call r3mv(ainv,efieldc,efieldl)
end if
afieldpol=.false.
if ((abs(afieldc(1)).gt.epslat).or.&
   &(abs(afieldc(2)).gt.epslat).or.&
   &(abs(afieldc(3)).gt.epslat)) then
  afieldpol=.true.
! vector potential added in second-variational step
  tevecsv=.true.
end if

!---------------------------------!
!     crystal symmetry set up     !
!---------------------------------!
! find Bravais lattice symmetries
call findsymlat
! use only the identity if required
if (nosym) nsymlat=1
! find the crystal symmetries and shift atomic positions if required
call findsymcrys
! find the site symmetries
call findsymsite
! check if fixed spin moments are invariant under the symmetry group
call checkfsm

!-------------------------!
!     LDA+U variables     !
!-------------------------!
if ((ldapu.ne.0).or.(task.eq.17)) then
! LDA+U requires second-variational eigenvectors
  tevecsv=.true.
! density matrices
  if (allocated(dmatlu)) deallocate(dmatlu)
  allocate(dmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! potential matrix elements
  if (allocated(vmatlu)) deallocate(vmatlu)
  allocate(vmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! zero the potential
  vmatlu(:,:,:,:,:)=0.d0
! energy for each atom
  if (allocated(engyalu)) deallocate(engyalu)
  allocate(engyalu(natmtot))
! interpolation constants (alpha)
  if (allocated(alphalu)) deallocate(alphalu)
  allocate(alphalu(natmtot))
end if

!---------------------!
!     k-point set     !
!---------------------!
! check if the system is an isolated molecule
if (molecule) then
  ngridk(:)=1
  vkloff(:)=0.d0
  autokpt=.false.
end if
! store the point group symmetries for reducing the k-point set
if (reducek.eq.0) then
  nsymkpt=1
  symkpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    if (reducek.eq.2) then
! check symmetry is symmorphic if required
      t1=abs(vtlsymc(1,isym))+abs(vtlsymc(2,isym))+abs(vtlsymc(3,isym))
      if (t1.gt.epslat) goto 10
! check also that the spin rotation is the same as the spatial rotation
      if (spinpol) then
        if (lspnsymc(isym).ne.lsplsymc(isym)) goto 10
      end if
    end if
    lsym(lsplsymc(isym))=.true.
10 continue
  end do
  nsymkpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymkpt=nsymkpt+1
      symkpt(:,:,nsymkpt)=symlat(:,:,isym)
    end if
  end do
end if

! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))

if (task.eq.2002) then
! for band structure plots generate k-points along a line
  call connect(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
  nkpt=npp1d
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,nkpt))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,nkpt))
  do ik=1,nkpt
    vkl(:,ik)=vplp1d(:,ik)
    call r3mv(bvec,vkl(:,ik),vkc(:,ik))
  end do
else
! setup the default k-point box
  boxl(:,1)=vkloff(:)/dble(ngridk(:))
  boxl(:,2)=boxl(:,1); 
  boxl(:,3)=boxl(:,1); 
  boxl(:,4)=boxl(:,1)
  boxl(1,2)=boxl(1,2)+1.d0
  boxl(2,3)=boxl(2,3)+1.d0
  boxl(3,4)=boxl(3,4)+1.d0
! allocate the reduced k-point set arrays
  if (allocated(ivk)) deallocate(ivk)
  allocate(ivk(3,ngridk(1)*ngridk(2)*ngridk(3)))
  if (allocated(vkl)) deallocate(vkl)
  allocate(vkl(3,ngridk(1)*ngridk(2)*ngridk(3)))
  if (allocated(vkc)) deallocate(vkc)
  allocate(vkc(3,ngridk(1)*ngridk(2)*ngridk(3)))
  if (allocated(wkpt)) deallocate(wkpt)
  allocate(wkpt(ngridk(1)*ngridk(2)*ngridk(3)))
  if (allocated(ikmap)) deallocate(ikmap)
  allocate(ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))
! generate the reduced k-point set
  call genppts(.false.,nsymkpt,symkpt,ngridk,epslat,bvec,boxl,nkpt,ikmap,ivk, &
              &vkl,vkc,wkpt)
endif

!
! set the basis parameters of the calculation: lattice vectors, number of atoms, maximum l-values
!   and plane-wave cutoffs
!
call sirius_clear

call sirius_set_lattice_vectors(avec(1,1), avec(1,2), avec(1,3))

call sirius_set_pw_cutoff(gmaxvr)

call sirius_set_aw_cutoff(rgkmax)

if (autormt) call sirius_set_auto_rmt(1)

do is=1,nspecies
  call sirius_add_atom_type(is, trim(spfname(is)))
enddo

do is=1,nspecies
  do ia=1,natoms(is)
    call sirius_add_atom(is, atposl(1,ia,is), bfcmt0(1,ia,is))
  enddo
enddo

! set equivalent atoms
ncls=0
allocate(icls(natmtot))
icls=0
do is=1,nspecies
  do ia1=1,natoms(is)
    if (icls(idxas(ia1,is)).eq.0) then
      ncls=ncls+1
      do ia2=1,natoms(is) 
        if (eqatoms(ia1,ia2,is)) then
          icls(idxas(ia2,is))=ncls
        endif
      enddo
    endif
  enddo
enddo
call sirius_set_equivalent_atoms(icls)
deallocate(icls)

if (ldapu.ne.0) call sirius_set_uj_correction(1)

call sirius_global_initialize(lmaxapw,lmaxvr,lmaxvr,ndmag)

call sirius_create_kset(nkpt,vkl,wkpt,1,kset_id)

! now it's time to get dimensions of important arrays that are allocated on the Fortran side
call sirius_get_num_grid_points(ngrtot)
call sirius_get_max_num_mt_points(nrmtmax)
call sirius_get_num_bands(nstsv)
call sirius_get_fft_grid_size(ngrid)
do i=1,3
  call sirius_get_fft_grid_limits(i,intgv(i,1),intgv(i,2))
enddo
call sirius_get_num_gvec(ngvec)

! get muffin-tin grids
if (allocated(spr)) deallocate(spr)
allocate(spr(nrmtmax,nspecies))
do is=1,nspecies
  call sirius_get_num_mt_points(is,nrmt(is))
  call sirius_get_mt_points(is,spr(1,is))
enddo

! get FFT related arrays
if (allocated(igfft)) deallocate(igfft)
allocate(igfft(ngrtot))
call sirius_get_fft_index(igfft)
if (allocated(ivg)) deallocate(ivg)
allocate(ivg(3,ngrtot))
call sirius_get_gvec(ivg)
if (allocated(ivgig)) deallocate(ivgig)
allocate(ivgig(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3,2)))
call sirius_get_index_by_gvec(ivgig)

! allocate memory for charge density 
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(lmmaxvr,nrmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngrtot))
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(magir(ngrtot,ndmag))
end if
!  allocate memory for effective potential
if (allocated(veffmt)) deallocate(veffmt)
allocate(veffmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(veffir)) deallocate(veffir)
allocate(veffir(ngrtot))
! XC magnetic field
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (spinpol) then
  allocate(bxcmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(bxcir(ngrtot,ndmag))
end if

call sirius_density_initialize(rhomt, rhoir, magmt, magir)
call sirius_potential_initialize(veffmt, veffir, bxcmt, bxcir)

if (allocated(occsv)) deallocate(occsv)
allocate(occsv(nstsv,nkpt))

if (allocated(evalsv)) deallocate(evalsv)
allocate(evalsv(nstsv,nkpt))

call sirius_get_num_electrons(chgtot)
call sirius_get_num_valence_electrons(chgval)
call sirius_get_num_core_electrons(chgcr)

! effective Wigner radius
rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)

return
end subroutine
#endif

