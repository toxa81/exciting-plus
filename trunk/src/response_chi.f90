subroutine response_chi(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)

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
complex(8), allocatable :: chi0(:,:,:)
! true polarisability
complex(8), allocatable :: chi(:,:)
! number of G-vectors for chi
integer ngvec_chi
! G+q vectors in Cartesian coordinates
real(8), allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: gq0(:)
! Coulomb potential 
real(8), allocatable :: vc(:)

! allocatable arrays
complex(8), allocatable :: mtrx1(:,:)
integer, allocatable :: ipiv(:)

integer ie,ig,ngsh_me_,info,i,j
character*100 fname
character*4 name1,name2,name3

if (iproc.eq.0) then
  write(fname,'("CHI0[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
    ivq0m(1),ivq0m(2),ivq0m(3)
  open(160,file=trim(fname),form='unformatted',status='old')
  read(160)ngsh_me_,ngvec_me,nepts,igq0
  if (ngsh_me_.ne.ngsh_me) then
    write(*,*)
    write(*,'("Error(response_chi): number of G-shells for matrix elements &
      &was changed")')
    write(*,*)
    call pstop
  endif
  allocate(w(nepts))
  allocate(chi0(ngvec_me,ngvec_me,nepts))
  read(160)w(1:nepts)
  read(160)vq0l(1:3)
  read(160)vq0rl(1:3)
  read(160)vq0c(1:3)
  read(160)vq0rc(1:3)
  do ie=1,nepts
    read(160)chi0(1:ngvec_me,1:ngvec_me,ie)
  enddo
  close(160)
  
! find number of G-vectors by given number of G-shells
  call getngvec(ngsh_chi,ngvec_chi)

  if (igq0.gt.ngvec_chi) then
    write(*,*)
    write(*,'("Error(response_chi): not enough G-vectors for calculation of &
      &chi")')
    write(*,'("  Increase number of G-shells for chi")')
    write(*,*)
    call pstop
  endif
  
  allocate(chi(ngvec_chi,nepts))
  
  allocate(vgq0c(3,ngvec_chi))
  allocate(vc(ngvec_chi))
  allocate(gq0(ngvec_chi))
  allocate(mtrx1(ngvec_chi,ngvec_chi))
  
  write(150,*)
  write(150,'("Coulomb potential matrix elements:")')
  write(150,'("ig    |G+q|    V")')
  write(150,'("----------------")')
  do ig=1,ngvec_chi
! generate G+q vectors  
    vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
    gq0(ig)=sqrt(vgq0c(1,ig)**2+vgq0c(2,ig)**2+vgq0c(3,ig)**2)
    vc(ig)=fourpi/gq0(ig)**2 
    write(150,'(I4,2x,2G18.10)')ig,gq0(ig),vc(ig)
  enddo

  allocate(ipiv(ngvec_chi))
  do ie=1,nepts
! compute 1-chi0*V
    do i=1,ngvec_chi
      do j=1,ngvec_chi
        mtrx1(i,j)=-chi0(i,j,ie)*vc(j)
      enddo
      mtrx1(i,i)=dcmplx(1.d0,0.d0)+mtrx1(i,i)
    enddo
! solve [1-chi0*V]^{-1}*chi=chi0
    chi(1:ngvec_chi,ie)=chi0(1:ngvec_chi,igq0,ie)
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
  write(160,'("#")')
  write(160,'("# Definition of columns")')
  write(160,'("#   1: energy            [eV]")')
  write(160,'("#   2: -Re chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   3: -Im chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   4: -Re chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   5: -Im chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#")')
  do ie=1,nepts
    write(160,'(5F18.10)')dreal(w(ie))*ha2ev,-dreal(chi0(igq0,igq0,ie)),&
      -dimag(chi0(igq0,igq0,ie)),-dreal(chi(igq0,ie)),-dimag(chi(igq0,ie))
  enddo
  close(160)
 
  deallocate(w)
  deallocate(chi0)
  deallocate(chi)
  deallocate(vgq0c)
  deallocate(vc)
  deallocate(gq0)
  deallocate(mtrx1)

endif

return
end  
