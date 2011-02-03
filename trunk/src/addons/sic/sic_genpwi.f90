subroutine sic_genpwi
use modmain
use mod_sic
implicit none
integer n,j,ik,ikloc,ir,itp,ispn,ig,ngvec1,lm
real(8) x(3),x1(3)
integer, allocatable :: stepf(:)
complex(8) zt1(nspinor),pwtp,expikr
complex(8), allocatable :: wtp(:,:)
complex(8), allocatable :: wvtp(:,:)
complex(8), allocatable :: expigr(:)
complex(8), allocatable :: expigkx0(:,:,:)
!real(8), allocatable :: stepf(:,:)
complex(8), external :: zdotu
!
allocate(wtp(sic_wantran%nwan,nspinor))
allocate(wvtp(sic_wantran%nwan,nspinor))
allocate(stepf(sic_wantran%nwan))
!allocate(stepf(s_ntp,s_nr))
!call s_gen_stepf(stepf)
ngvec1=0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do ig=1,ngk(1,ik)
    ngvec1=max(ngvec1,igkig(ig,1,ikloc))
  enddo
enddo
allocate(expigr(ngvec1))
! note: \int W_n^{*}(x-x0) e^{i(G+k)x} dx = [x-x0=x';x=x0+x'] 
!  = \int W_n^{*}(x') e^{i(G+k)(x0+x')} dx'
allocate(expigkx0(ngkmax,sic_wantran%nwan,nkptloc))
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    do ig=1,ngk(1,ik)
      expigkx0(ig,j,ikloc)=exp(zi*dot_product(vgkc(:,ig,1,ikloc),wanpos(:,n)))
    enddo
  enddo
enddo
sic_wgk=zzero
sic_wvgk=zzero
do ir=1,s_nr
  do itp=1,s_ntp
    x(:)=s_spx(:,itp)*s_r(ir) 
    do ig=1,ngvec1
      expigr(ig)=exp(zi*dot_product(vgc(:,ig),x(:)))*s_tpw(itp)*&
        s_rw(ir)/sqrt(omega)
    enddo
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      do ispn=1,nspinor
        wtp(j,ispn)=zdotu(lmmaxwan,s_wanlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)
        wvtp(j,ispn)=zdotu(lmmaxwan,s_wvlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)
      enddo
      x1(:)=x(:)+wanpos(:,n)
      stepf(j)=s_gen_stepf(x1)
    enddo
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      expikr=exp(zi*dot_product(vkc(:,ik),x(:)))
      do ispn=1,nspinor
        do j=1,sic_wantran%nwan
          if (stepf(j).ne.0) then
            do ig=1,ngk(1,ik)
              sic_wgk(ig,j,ispn,ikloc)=sic_wgk(ig,j,ispn,ikloc)+&
                expigr(igkig(ig,1,ikloc))*expikr*dconjg(wtp(j,ispn))*&
                expigkx0(ig,j,ikloc)
              sic_wvgk(ig,j,ispn,ikloc)=sic_wvgk(ig,j,ispn,ikloc)+&
                expigr(igkig(ig,1,ikloc))*expikr*dconjg(wvtp(j,ispn))*&
                expigkx0(ig,j,ikloc)
            enddo !ig
          endif
        enddo
      enddo
    enddo !ikloc
  enddo !itp
enddo !ir
deallocate(wtp,wvtp)
deallocate(stepf)
deallocate(expigr)
deallocate(expigkx0)
return
end
