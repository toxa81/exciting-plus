subroutine init_q_gq
use modmain
use mod_addons_q
implicit none
real(8) t0,t2
integer iq,ig,i
integer v1(3)
real(8) v2(3)
logical tgsh

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

tgsh=.false.

t0=gqmax**2
ngqmax=0
do iq=1,nvq
  vqlnr(:,iq)=dble(vqm(:,iq))/ngridk(:)
  do ig=1,ngvec
    v1(:)=vqm(:,iq)-ngridk(:)*ivg(:,ig)
    if (v1(1).ge.0.and.v1(1).lt.ngridk(1).and.&
        v1(2).ge.0.and.v1(2).lt.ngridk(2).and.&
        v1(3).ge.0.and.v1(3).lt.ngridk(3)) then
      ig0q(iq)=ig
      vql(:,iq)=dble(v1(:))/ngridk(:)
      goto 10
    endif
  enddo !ig
  write(*,*)
  write(*,'("Error(init_q_gq): G0-vector is not found because &
    &q-vector is too large")')
  call pstop
10 continue
  call r3mv(bvec,vqlnr(:,iq),vqcnr(:,iq))
  call r3mv(bvec,vql(:,iq),vqc(:,iq))
  i=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if (t2.le.t0) i=i+1
  enddo
  ngqmax=max(ngqmax,i)
enddo
if (tgsh) then
  call getngvecme
  ngqmax=ngvecme
endif
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
do iq=1,nvq
  ngq(iq)=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if ((.not.tgsh.and.t2.le.t0).or.(tgsh.and.ig.le.ngvecme)) then
      ngq(iq)=ngq(iq)+1
      gq(ngq(iq),iq)=sqrt(t2)
      vgqc(:,ngq(iq),iq)=v2(:)
      igqig(ngq(iq),iq)=ig
      if (all(vqm(:,iq).eq.0).and.ig.eq.1) then
        if (nvq0.eq.1) then
          vhgq(ngq(iq),iq)=fourpi*aq0(1)
        endif
      else
        vhgq(ngq(iq),iq)=fourpi/t2
      endif
    endif
  enddo !ig
enddo    

return
end