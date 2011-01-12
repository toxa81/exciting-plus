subroutine sic_genfvprj_ir
use modmain
use mod_sic
implicit none
! local varibales
integer ik,ikloc,h,j,ispn,it,ig,ir
real(8) d1
integer ngk1
integer, allocatable :: igkig1(:)
real(8), allocatable :: gkc1(:)
real(8), allocatable :: tpgkc1(:,:)
real(8), allocatable :: vgkl1(:,:)
real(8), allocatable :: vgkc1(:,:)
complex(8), allocatable :: wanirk(:)
complex(8), allocatable :: wvirk(:)
complex(8), allocatable :: wgk(:)
complex(8), allocatable :: wvgk(:)
!
allocate(igkig1(ngkmax))
allocate(gkc1(ngkmax))
allocate(tpgkc1(2,ngkmax))
allocate(vgkl1(3,ngkmax))
allocate(vgkc1(3,ngkmax))
allocate(wanirk(ngrloc))
allocate(wvirk(ngrloc))
allocate(wgk(ngkmax))
allocate(wvgk(ngkmax))
sic_wgk=zzero
sic_wvgk=zzero
do ik=1,nkpt
  ikloc=mpi_grid_map(nkpt,dim_k,x=h,glob=ik)
  call gengpvec(vkl(1,ik),vkc(1,ik),ngk1,igkig1,vgkl1,vgkc1,gkc1,tpgkc1)
  do j=1,sic_wantran%nwan
    do ispn=1,nspinor
      wanirk=zzero
      wvirk=zzero
      do it=1,sic_orbitals%ntr
        d1=dot_product(vkc(:,ik),sic_orbitals%vtc(:,it))
        wanirk(:)=wanirk(:)+dconjg(sic_orbitals%wanir(:,it,ispn,j))*exp(zi*d1)
        wvirk(:)=wvirk(:)+dconjg(sic_orbitals%wvir(:,it,ispn,j))*exp(zi*d1)
      enddo !it
      wgk=zzero
      wvgk=zzero
      do ig=1,ngk1
        do ir=1,ngrloc
          d1=dot_product(vgkc1(:,ig),vgrc(:,ir+groffs))
          wgk(ig)=wgk(ig)+wanirk(ir)*cfunir(ir+groffs)*exp(zi*d1)/sqrt(omega)
          wvgk(ig)=wvgk(ig)+wvirk(ir)*cfunir(ir+groffs)*exp(zi*d1)/sqrt(omega)
        enddo !ir
        wgk(ig)=wgk(ig)*omega/ngrtot
        wvgk(ig)=wvgk(ig)*omega/ngrtot
      enddo !ig
      call mpi_grid_reduce(wgk(1),ngkmax,outb=sic_wgk(1,j,ispn,ikloc),&
        root=(/h/))
      call mpi_grid_reduce(wvgk(1),ngkmax,outb=sic_wvgk(1,j,ispn,ikloc),&
        root=(/h/))
    enddo !ispn
  enddo !j
enddo !ik
deallocate(igkig1,gkc1,tpgkc1,vgkl1,vgkc1,wanirk,wvirk,wgk,wvgk)
return
end
