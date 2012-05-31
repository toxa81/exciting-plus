module mod_addons_q
!-------------------------!
!     q and G+q vectors   !
!-------------------------!
! number of q-vectors
integer nvq
! q-vectors in k-mesh coordinates
integer, allocatable :: vqm(:,:)
! non-reduced (to first BZ) q-vectors in lattice coordinates
real(8), allocatable :: vqlnr(:,:)
! non-reduced (to first BZ) q-vectors in Cartesian coordinates
real(8), allocatable :: vqcnr(:,:)
! reduced to first BZ q-vectors in lattice coordinates
real(8), allocatable :: vql(:,:)
! reduced to first BZ q-vectors in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! cutoff for |G+q|
real(8) gqmax
data gqmax/2.d0/
integer gqsh
data gqsh/2/
logical tgqsh
data tgqsh/.false./
! maximum number of G+q vectors
integer ngqmax
! global index of Gq-vector, which brigs q-vector to first BZ
integer, allocatable :: ig0q(:)
! index of Gq vector in the range[1,ngq(iq)]
integer iig0q
! number of G+q vectors
integer, allocatable :: ngq(:)
! G+q vectors in Cartesian coordinates
real(8), allocatable ::  vgqc(:,:,:)
! length of |G+q| vectors
real(8), allocatable :: gq(:,:)
! global index of G of G+q vector  
integer, allocatable :: igqig(:,:)
! 4*Pi/|G+q|^2 (Fourier transform of Hartree potential)
real(8), allocatable :: vhgq(:,:)
! weighted G+q components of Hartee potential
real(8), allocatable :: wtvhgq(:,:)

real(8), allocatable :: tpgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: ylmgq(:,:)

! true if Gamma point included in the integration over Brillouin zone
logical tq0bz
data tq0bz/.true./
integer nvq0
data nvq0/0/
real(8) vq0c(3)
real(8) q0wt

contains 

logical function vq_gamma(iq)
implicit none
integer, intent(in) :: iq
!
vq_gamma=.false.
if (sum(vql(:,iq)**2).lt.1d-10) vq_gamma=.true.
return
end function

integer function getngvecme()
use modmain
implicit none
integer ngsh
integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))
call getgshells(ngsh,igishell,ishellng)
getngvecme=ishellng(gqsh,2)
deallocate(igishell)
deallocate(ishellng)
return
end function

subroutine init_q_mesh(nvq0_)
use modmain
implicit none
integer, intent(in) :: nvq0_
integer j,i1,i2,i3
!
nvq0=nvq0_
if (allocated(vqm)) deallocate(vqm)
nvq=nkptnr-1+nvq0
allocate(vqm(3,nvq))
vqm=0
j=nvq0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        vqm(:,j)=(/i1,i2,i3/)
      endif
    enddo
  enddo
enddo
return
end subroutine

! generate q-vector related data
subroutine genvq
use modmain
implicit none
integer iq,ig,v1(3),i,i1,i2,i3
logical f
!
if (allocated(vqlnr)) deallocate(vqlnr)
allocate(vqlnr(3,nvq))
if (allocated(vqcnr)) deallocate(vqcnr)
allocate(vqcnr(3,nvq))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nvq))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nvq))
if (allocated(ig0q)) deallocate(ig0q)
allocate(ig0q(nvq))

do iq=1,nvq
! non-reduced q-vector in lattice coordinates
  vqlnr(:,iq)=dble(vqm(:,iq))/dble(ngridk(:))
! non-reduced q-vector in Cartesian coordinates 
  call r3mv(bvec,vqlnr(:,iq),vqcnr(:,iq))
! find Gq vector which reduces q to first BZ
  f=.false.
  do ig=1,ngvec
! q-G vector in k-mesh coordinates 
    v1(:)=vqm(:,iq)-ngridk(:)*ivg(:,ig)
    if (v1(1).ge.0.and.v1(1).lt.ngridk(1).and.&
        v1(2).ge.0.and.v1(2).lt.ngridk(2).and.&
        v1(3).ge.0.and.v1(3).lt.ngridk(3).and..not.f) then
      ig0q(iq)=ig
      f=.true.
    endif
  enddo !ig
  if (.not.f) then
    write(*,*)
    write(*,'("Error(genvq): Gq-vector is not found because &
      &q-vector is too large")')
    call pstop
  endif
! reduced q-vector in lattice coordinates
  vql(:,iq)=vqlnr(:,iq)-dble(ivg(:,ig0q(iq)))
! reduced q-vector in Cartesian coordinates
  call r3mv(bvec,vql(:,iq),vqc(:,iq))
enddo
! init q=0 points
if (nvq0.eq.1) then
  if (.not.vq_gamma(1)) then
    write(*,'("[genvq] : firt q-pont in the list is not a Gamma")')
    call pstop
  endif
  vqc(:,1)=vq0c(:)
else if (nvq0.eq.8) then
  i=0
  do i1=0,1
    do i2=0,1
      do i3=0,1
        i=i+1
! take small vector on the diagonal of the eight nearest corner-sharing micro-cells
        vqc(:,i)=0.01d0*((i1-0.5d0)*bvec(:,1)/ngridk(1)+(i2-0.5d0)*bvec(:,2)/ngridk(2)+&
          &(i3-0.5d0)*bvec(:,3)/ngridk(3))
      enddo
    enddo
  enddo
else if (nvq0.ne.0) then
  write(*,'("[genvq] : nvq0= ",I4," is not implemented")')nvq0
  call pstop
endif
return
end subroutine

! generate G+q-vector related data
subroutine genvgq
use modmain
!use mod_expigqr
implicit none
real(8) t2,gqmax2
integer iq,ig,i
real(8) v2(3)
logical tautogqmax
!
! find |G+q| cutoff
tautogqmax=.true.
if (tautogqmax) then
  do iq=1,nvq
    t2=sqrt(sum(vqcnr(:,iq)**2))
    gqmax=max(gqmax,t2*1.01)
  enddo
endif
! find maximum number of G+q vectors
ngqmax=0
gqmax2=gqmax**2
do iq=1,nvq
  i=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if (t2.le.gqmax2) i=i+1
  enddo
  ngqmax=max(ngqmax,i)
enddo
if (tgqsh) then
  ngqmax=getngvecme()
endif
! generate G+q vectors
if (allocated(ngq)) deallocate(ngq)
allocate(ngq(nvq))
if (allocated(gq)) deallocate(gq)
allocate(gq(ngqmax,nvq))
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngqmax,nvq))
if (allocated(igqig)) deallocate(igqig)
allocate(igqig(ngqmax,nvq))
if (allocated(vhgq)) deallocate(vhgq)
allocate(vhgq(ngqmax,nvq))
vhgq(:,:)=0.d0
do iq=1,nvq
  ngq(iq)=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if ((.not.tgqsh.and.t2.le.gqmax2).or.(tgqsh.and.ig.le.ngqmax)) then
      ngq(iq)=ngq(iq)+1
      gq(ngq(iq),iq)=sqrt(t2)
      vgqc(:,ngq(iq),iq)=v2(:)
      igqig(ngq(iq),iq)=ig
      vhgq(ngq(iq),iq)=fourpi/t2
    endif
  enddo !ig
enddo
! compute weights for 1/|q|^2 integral
if (nvq0.ne.0) call genwtvhgq
return
end subroutine

subroutine genwtvhgq
use modmain
implicit none
integer i,n,i1,i2,i3,ig,iq
real(8) x(2),y(2)
real(8) vq(3),p1,p2,p3,vol
!
if (allocated(wtvhgq)) deallocate(wtvhgq)
allocate(wtvhgq(ngqmax,nvq))
!
n=100
y(:)=0.d0
do i=1,2
  do i1=-n,n-1
    p1=(0.5d0+i1)/dble(2*n)
    do i2=-n,n-1
      p2=(0.5d0+i2)/dble(2*n)
      do i3=-n,n-1
        p3=(0.5d0+i3)/dble(2*n)
        vq(:)=(p1*bvec(:,1)/ngridk(1)+p2*bvec(:,2)/ngridk(2)+p3*bvec(:,3)/ngridk(3))
        y(i)=y(i)+1.d0/dot_product(vq,vq)
      enddo
    enddo
  enddo
  vol=(twopi**3)/omega/nkptnr/((2*n)**3)
  y(i)=y(i)*vol
  x(i)=vol**(1/3.d0)
  n=n+40
enddo
! numerical derivative to get y(x=0), where x is the effective size of the micro-cell
q0wt=y(2)-x(2)*(y(1)-y(2))/(x(1)-x(2))

wtvhgq(:,:)=vhgq(:,:)
do iq=1,nvq
  if (vq_gamma(iq)) then
    do ig=1,ngq(iq)
      if (igqig(ig,iq).eq.1) then
        wtvhgq(ig,iq)=fourpi*q0wt*nkptnr*omega/(twopi**3)
      endif
    enddo
    wtvhgq(:,iq)=wtvhgq(:,iq)/dble(nvq0)
  endif
enddo
return
end subroutine

end module
