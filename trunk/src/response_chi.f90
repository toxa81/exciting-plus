subroutine response_chi(ivq0m,gvecchi1,gvecchi2)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)
! number of G-vectors for chi
integer, intent(inout) :: gvecchi1
integer, intent(inout) :: gvecchi2

integer gvecme1
integer gvecme2
integer ngvecchi
! number of G-vectors for matrix elements
integer ngvecme
! number of energy-mesh points
integer nepts
! energy mesh
complex(8), allocatable :: w(:)
! q-vector in lattice coordinates
real(8) vq0l(3)
! q-vector in Cartesian coordinates
real(8) vq0c(3)
! reduced q-vector in lattice coordinates
real(8) vq0rl(3)
! reduced q-vector in Cartesian coordinates
real(8) vq0rc(3)
! index of G-vector which brings q to first BZ
integer igq0
integer igq0s
! Kohn-Sham polarisability
complex(8), allocatable :: chi0(:,:,:,:)
complex(8), allocatable :: chi0s(:,:,:)
! true polarisability
complex(8), allocatable :: chi(:,:)
! G+q vectors in Cartesian coordinates
real(8), allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: gq0(:)
! Coulomb potential 
real(8), allocatable :: vc(:)

complex(8), allocatable :: epsilon_GqGq(:)
complex(8), allocatable :: chi_scalar(:)
complex(8), allocatable :: epsilon_eff(:)
complex(8), allocatable :: mtrx(:,:)
integer nspin_chi0

! allocatable arrays
real(8), allocatable :: func(:,:)
complex(8), allocatable :: epsilon(:,:)
integer, allocatable :: ipiv(:)

integer ie,ig,ngsh_me_,info,i,j,ig1,ig2
character*100 fname
character*4 name1,name2,name3

if (iproc.eq.0) then

  write(150,*)
  write(150,'("Calculation of charge polarizability chi")')  

  write(fname,'("CHI0[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
    ivq0m(1),ivq0m(2),ivq0m(3)
  write(150,'("Reading file ",A40)')trim(fname)
  open(160,file=trim(fname),form='unformatted',status='old')
  read(160)nepts,igq0
  read(160)gshme1,gshme2,gvecme1,gvecme2,ngvecme
  allocate(w(nepts))
  read(160)w(1:nepts)
  read(160)vq0l(1:3)
  read(160)vq0rl(1:3)
  read(160)vq0c(1:3)
  read(160)vq0rc(1:3)
  read(160)spin_me,nspin_chi0
  allocate(chi0(ngvecme,ngvecme,nepts,nspin_chi0))
  do ie=1,nepts
    read(160)chi0(1:ngvecme,1:ngvecme,ie,1:nspin_chi0)
  enddo
  if (task.eq.402) close(160)
  if (task.eq.403) close(160,status='delete')
  write(150,'("chi0 was calculated for ")')
  write(150,'("  G-shells  : ",I4," to ",I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ",I4)')gvecme1,gvecme2
  if (spinpol) then 
    if (spin_me.eq.1) write(150,'("chi0 was calculated for spin up")')
    if (spin_me.eq.2) write(150,'("chi0 was calculated for spin dn")')
    if (spin_me.eq.3) write(150,'("chi0 was calculated for both spins")')
  endif
  if (gvecchi1.lt.gvecme1.or.gvecchi1.gt.gvecme2) then
    write(150,*)
    write(150,'("Warning: minimum number of G-vectors was changed from ",&
      &I4," to ",I4)')gvecchi1,gvecme1
    gvecchi1=gvecme1 
  endif
  if (gvecchi2.lt.gvecme1.or.gvecchi2.gt.gvecme2) then
    write(150,*)
    write(150,'("Warning: maximum number of G-vectors was changed from ",&
      &I4," to ",I4)')gvecchi2,gvecme2
    gvecchi2=gvecme2 
  endif
  if (igq0.lt.gvecchi1.or.igq0.gt.gvecchi2) then
    write(*,*)
    write(*,'("Error(response_chi): not enough G-vectors for calculation of &
      &chi")')
    write(*,*)
    call pstop
  endif
  if (spin_me.le.2.and.spin_chi.ne.spin_me) then
    write(*,*)
    write(*,'("Error(response_chi): spin_chi != spin_me")')
    write(*,'("  spin_chi = ",I4)')spin_chi
    write(*,'("  spin_me = ",I4)')spin_me
    write(*,*)
    call pstop
  endif
  if (spin_me.le.2.and.spin_chi.eq.3) then
    write(*,*)
    write(*,'("Error(response_chi): chi0 is calculated only for one spin")')
    write(*,'("  can''t make sum of chi0(up) and chi0(dn)")')
    write(*,*)
    call pstop
  endif
  
  if (spinpol.and.spin_me.eq.3) then
    if (spin_chi.eq.1) write(150,'("using chi0(up)")')
    if (spin_chi.eq.2) then
      write(150,'("using chi0(dn)")')
      chi0(:,:,:,1)=chi0(:,:,:,2)
    endif
    if (spin_chi.eq.3) then
      write(150,'("using chi0(up)+chi0(dn)")')
      chi0(:,:,:,1)=chi0(:,:,:,1)+chi0(:,:,:,2)
    endif
  endif
  
  if (spinpol.and.afmchi0.and.(spin_chi.le.2)) then
    write(150,'("AFM case: chi0 is multiplied by 2")')
    chi0=chi0*2.d0
  endif  
  
  ngvecchi=gvecchi2-gvecchi1+1  
  write(150,*)
  write(150,'("Minimum and maximum G-vectors for chi : ",2I4)')gvecchi1,gvecchi2
  write(150,'("Number of G-vectors : ",I4)')ngvecchi
  
  allocate(chi0s(ngvecchi,ngvecchi,nepts))
  do ie=1,nepts
    ig1=gvecchi1-gvecme1+1
    ig2=ig1+ngvecchi-1
    chi0s(1:ngvecchi,1:ngvecchi,ie)=chi0(ig1:ig2,ig1:ig2,ie,1)
  enddo
  igq0s=igq0-gvecchi1+1
  
  allocate(chi(ngvecchi,nepts))
  allocate(vgq0c(3,ngvecchi))
  allocate(vc(ngvecchi))
  allocate(gq0(ngvecchi))
  allocate(epsilon(ngvecchi,ngvecchi))
  allocate(epsilon_GqGq(nepts))
  allocate(epsilon_eff(nepts))
  allocate(chi_scalar(nepts))
  allocate(mtrx(ngvecchi,ngvecchi))
  
  write(150,*)
  write(150,'("Coulomb potential matrix elements:")')
  write(150,'("   ig        |G+q|        Vc    ")')
  write(150,'(" ------------------------------ ")')
  do ig=1,ngvecchi
! generate G+q vectors  
    vgq0c(:,ig)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
    gq0(ig)=sqrt(vgq0c(1,ig)**2+vgq0c(2,ig)**2+vgq0c(3,ig)**2)
    vc(ig)=fourpi/gq0(ig)**2 
    write(150,'(1X,I4,2X,2F12.6)')ig,gq0(ig),vc(ig)
  enddo
  
  allocate(ipiv(ngvecchi))
  do ie=1,nepts
! compute epsilon=1-chi0*V
    do i=1,ngvecchi
      do j=1,ngvecchi
        epsilon(i,j)=-chi0s(i,j,ie)*vc(j)
      enddo
      epsilon(i,i)=dcmplx(1.d0,0.d0)+epsilon(i,i)
    enddo
    epsilon_GqGq(ie)=epsilon(igq0s,igq0s)
    chi_scalar(ie)=chi0s(igq0s,igq0s,ie)/epsilon_GqGq(ie)
! compute epsilon effective
    mtrx(:,:)=epsilon(:,:)
    call invzge(mtrx,ngvecchi)
    epsilon_eff(ie)=1.d0/mtrx(igq0s,igq0s)
! solve [1-chi0*V]^{-1}*chi=chi0
    chi(1:ngvecchi,ie)=chi0s(1:ngvecchi,igq0s,ie)
    call zgesv(ngvecchi,1,epsilon,ngvecchi,ipiv,chi(1,ie),ngvecchi,info)
    if (info.ne.0) then
      write(*,*)'Error solving linear equations'
      write(*,*)'info=',info
      stop
    endif
  enddo !ie
  deallocate(ipiv)
  
  chi0s=chi0s/ha2ev/(au2ang)**3
  chi=chi/ha2ev/(au2ang)**3
  chi_scalar=chi_scalar/ha2ev/(au2ang)**3
  
  write(fname,'("response[",I4.3,",",I4.3,",",I4.3,"].dat")') &
    ivq0m(1),ivq0m(2),ivq0m(3)

  open(160,file=trim(fname),form='formatted',status='replace')
  write(160,'("# k-mesh division                    : ",3I4)')ngridk(1),ngridk(2),ngridk(3)
  write(160,'("# Energy mesh parameters             : ")')
  write(160,'("#   maximum energy [eV]              : ", F7.2)')maxomega
  write(160,'("#   energy step    [eV]              : ", F7.2)')domega
  write(160,'("#   eta            [eV]              : ", F7.2)')eta
  write(160,'("# q-vector information               : ")')
  write(160,'("#   q-vector (mesh coord.)           : ",3I4)')ivq0m
  write(160,'("#   q-vector (lat. coord.)           : ",3F18.10)')vq0l
  write(160,'("#   q-vector (Cart. coord.) [a.u.]   : ",3F18.10)')vq0c
  write(160,'("#   q-vector length         [a.u.]   : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
  write(160,'("#   q-vector (Cart. coord.) [1/A]    : ",3F18.10)')vq0c/au2ang
  write(160,'("#   q-vector length         [1/A]    : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
  write(160,'("# G-vector information               : ")')
  write(160,'("#   G-shells                         : ",2I4)')gshchi1,gshchi2
  write(160,'("#   G-vectors                        : ",2I4)')gvecchi1,gvecchi2
  write(160,'("#   index of Gq vector               : ",I4)')igq0
  write(160,'("#")')
  write(160,'("# Definition of columns")')
  write(160,'("#   1: energy            [eV]")')
  write(160,'("#   2: -Re chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   3: -Im chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   4: -Re chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   5: -Im chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   6:  S(q,w)           [1/eV/A^3]")')
  write(160,'("#   7: -Re chi_scalar    [1/eV/A^3]")')
  write(160,'("#   8: -Im chi_scalar    [1/eV/A^3]")')
  write(160,'("#   9:  Re epsilon_eff             ")')
  write(160,'("#  10:  Im epsilon_eff             ")')
  write(160,'("#  11:  Re epsilon_GqGq            ")')
  write(160,'("#  12:  Im epsilon_GqGq            ")')
  write(160,'("#")')
  allocate(func(12,nepts))
  do ie=1,nepts
    func(1,ie)=dreal(w(ie))*ha2ev
    func(2,ie)=-dreal(chi0s(igq0s,igq0s,ie))
    func(3,ie)=-dimag(chi0s(igq0s,igq0s,ie))
    func(4,ie)=-dreal(chi(igq0s,ie))
    func(5,ie)=-dimag(chi(igq0s,ie))
    func(6,ie)=-2.d0*dimag(chi(igq0s,ie))
    func(7,ie)=-dreal(chi_scalar(ie))
    func(8,ie)=-dimag(chi_scalar(ie))
    func(9,ie)=dreal(epsilon_eff(ie))
    func(10,ie)=dimag(epsilon_eff(ie))
    func(11,ie)=dreal(epsilon_GqGq(ie))
    func(12,ie)=dimag(epsilon_GqGq(ie))
    write(160,'(16F16.8)')func(1:12,ie)
  enddo
  deallocate(func)
  close(160)
 
  deallocate(w)
  deallocate(chi0)
  deallocate(chi0s)
  deallocate(chi)
  deallocate(vgq0c)
  deallocate(vc)
  deallocate(gq0)
  deallocate(epsilon)
  deallocate(epsilon_GqGq)
  deallocate(epsilon_eff)
  deallocate(chi_scalar)
  deallocate(mtrx)
  
  write(150,*)
  write(150,'("Done.")')

endif

return
end  
