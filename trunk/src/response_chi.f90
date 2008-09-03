subroutine response_chi(ivq0m,ngvec_chi)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)
! number of G-vectors for chi
integer, intent(in) :: ngvec_chi

! number of G-vectors for matrix elements
integer ngvec_me
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
! Kohn-Sham polarisability
complex(8), allocatable :: chi0(:,:,:,:)
! true polarisability
complex(8), allocatable :: chi(:,:)
! G+q vectors in Cartesian coordinates
real(8), allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: gq0(:)
! Coulomb potential 
real(8), allocatable :: vc(:)

complex(8), allocatable :: epsilon_GqGq(:)

integer nspin_chi0

! allocatable arrays
real(8), allocatable :: func(:,:)
complex(8), allocatable :: mtrx1(:,:)
integer, allocatable :: ipiv(:)

integer ie,ig,ngsh_me_,info,i,j
complex(8) chi_scalar
character*100 fname
character*4 name1,name2,name3

if (iproc.eq.0) then

  write(150,*)
  write(150,'("Calculation of charge polarizability chi")')  

  write(fname,'("CHI0[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
    ivq0m(1),ivq0m(2),ivq0m(3)
  write(150,'("Reading file ",A40)')trim(fname)
  open(160,file=trim(fname),form='unformatted',status='old')
  read(160)ngsh_me,ngvec_me,nepts,igq0
  allocate(w(nepts))
  read(160)w(1:nepts)
  read(160)vq0l(1:3)
  read(160)vq0rl(1:3)
  read(160)vq0c(1:3)
  read(160)vq0rc(1:3)
  read(160)spin_me,nspin_chi0
  allocate(chi0(ngvec_me,ngvec_me,nepts,nspin_chi0))
  do ie=1,nepts
    read(160)chi0(1:ngvec_me,1:ngvec_me,ie,1:nspin_chi0)
  enddo
  if (task.eq.402) close(160)
  if (task.eq.403) close(160,status='delete')
  
  write(150,'("chi0 was calculated for ",I4," G-vector(s) (",I4,&
    & " G-shell(s))")')ngvec_me,ngsh_me
  if (spinpol) then 
    if (spin_me.eq.1) write(150,'("chi0 was calculated for spin up")')
    if (spin_me.eq.2) write(150,'("chi0 was calculated for spin dn")')
    if (spin_me.eq.3) write(150,'("chi0 was calculated for both spins")')
  endif
  
  if (igq0.gt.ngvec_chi) then
    write(*,*)
    write(*,'("Error(response_chi): not enough G-vectors for calculation of &
      &chi")')
    write(*,'("  Increase number of G-shells for chi")')
    write(*,*)
    call pstop
  endif
  if (ngsh_chi.gt.ngsh_me) then
    write(*,*)
    write(*,'("Error(response_chi): ngsh_chi > ngsh_me")')
    write(*,'("  ngsh_chi = ",I4)'),ngsh_chi
    write(*,'("  ngsh_me = ",I4)'),ngsh_me
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
  
  allocate(chi(ngvec_chi,nepts))
  
  allocate(vgq0c(3,ngvec_chi))
  allocate(vc(ngvec_chi))
  allocate(gq0(ngvec_chi))
  allocate(mtrx1(ngvec_chi,ngvec_chi))
  allocate(epsilon_GqGq(nepts))
  
  write(150,*)
  write(150,'("Coulomb potential matrix elements:")')
  write(150,'("   ig        |G+q|        Vc    ")')
  write(150,'(" ------------------------------ ")')
  do ig=1,ngvec_chi
! generate G+q vectors  
    vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
    gq0(ig)=sqrt(vgq0c(1,ig)**2+vgq0c(2,ig)**2+vgq0c(3,ig)**2)
    vc(ig)=fourpi/gq0(ig)**2 
    write(150,'(1X,I4,2X,2F12.6)')ig,gq0(ig),vc(ig)
  enddo
  
  allocate(ipiv(ngvec_chi))
  do ie=1,nepts
! compute 1-chi0*V
    do i=1,ngvec_chi
      do j=1,ngvec_chi
        mtrx1(i,j)=-chi0(i,j,ie,1)*vc(j)
      enddo
      mtrx1(i,i)=dcmplx(1.d0,0.d0)+mtrx1(i,i)
    enddo
    epsilon_GqGq(ie)=mtrx1(igq0,igq0)
! solve [1-chi0*V]^{-1}*chi=chi0
    chi(1:ngvec_chi,ie)=chi0(1:ngvec_chi,igq0,ie,1)
    call zgesv(ngvec_chi,1,mtrx1,ngvec_chi,ipiv,chi(1,ie),ngvec_chi,info)
    if (info.ne.0) then
      write(*,*)'Error solving linear equations'
      write(*,*)'info=',info
      stop
    endif
  enddo !ie
  deallocate(ipiv)
  
  chi0=chi0/ha2ev/(au2ang)**3
  chi=chi/ha2ev/(au2ang)**3
!  epsilon_GqGq=epsilon_GqGq/(au2ang)**3
  
!  write(name1,'(I4.3)')ivq0l(1)
!  write(name2,'(I4.3)')ivq0l(2)
!  write(name3,'(I4.3)')ivq0l(3)
!  fname='response['//trim(adjustl(name1))//','//trim(adjustl(name2))//','//trim(adjustl(name3))//'].dat'

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
  write(160,'("#   number of G-shells               : ",I4)')ngsh_chi
  write(160,'("#   number of G-vectors              : ",I4)')ngvec_chi
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
  write(160,'("#   9:  Re epsilon_eff   [eV*A^3]  ")')
  write(160,'("#  10:  Im epsilon_eff   [eV*A^3]  ")')
  write(160,'("#  11:  Re epsilon_GqGq            ")')
  write(160,'("#  12:  Im epsilon_GqGq            ")')
  write(160,'("#")')
  allocate(func(12,nepts))
  do ie=1,nepts
    chi_scalar=chi0(igq0,igq0,ie,1)/(1.d0-vc(igq0)*chi0(igq0,igq0,ie,1))
    !epsilon_GqGq=1.d0-vc(igq0)*chi0(igq0,igq0,ie,1)
    func(1,ie)=dreal(w(ie))*ha2ev
    func(2,ie)=-dreal(chi0(igq0,igq0,ie,1))
    func(3,ie)=-dimag(chi0(igq0,igq0,ie,1))
    func(4,ie)=-dreal(chi(igq0,ie))
    func(5,ie)=-dimag(chi(igq0,ie))
    func(6,ie)=-2.d0*dimag(chi(igq0,ie))
    func(7,ie)=-dreal(chi_scalar)
    func(8,ie)=-dimag(chi_scalar)
    func(9,ie)=dreal(0.5d0/chi(igq0,ie))
    func(10,ie)=dimag(0.5d0/chi(igq0,ie))
    func(11,ie)=dreal(epsilon_GqGq(ie))
    func(12,ie)=dimag(epsilon_GqGq(ie))
    write(160,'(16F16.8)')func(1:12,ie)
  enddo
  deallocate(func)
  close(160)
 
  deallocate(w)
  deallocate(chi0)
  deallocate(chi)
  deallocate(vgq0c)
  deallocate(vc)
  deallocate(gq0)
  deallocate(mtrx1)
  deallocate(epsilon_GqGq)
  
  write(150,*)
  write(150,'("Done.")')

endif

return
end  
