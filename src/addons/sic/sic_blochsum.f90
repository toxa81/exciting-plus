subroutine sic_blochsum
use modmain
use mod_sic
implicit none
integer ik,ikloc,j,n,jas,it,ias,is,ir,itp,ispn
real(8) x(3)
real(8), allocatable :: tp(:,:)
complex(8) zt1,zt2,expikt
complex(8), external :: zfinp_
!
s_wankmt=zzero
s_wankir=zzero
s_wvkmt=zzero
s_wvkir=zzero

if (.not.tsic_wv) return

allocate(tp(2,lmmaxvr))
call sphcover(lmmaxvr,tp)

do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! make Bloch sums
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    jas=wan_info(1,n)
    do it=1,sic_orbitals%ntr
      expikt=exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
! muffin-tins
      do ias=1,natmtot
        is=ias2is(ias)
        do ir=1,nrmt(is)
          do itp=1,lmmaxvr
            x(:)=(/sin(tp(1,itp))*cos(tp(2,itp)), &
                   sin(tp(1,itp))*sin(tp(2,itp)), &
                   cos(tp(1,itp))/)*spr(ir,is)+atposc(:,ias2ia(ias),ias2is(ias))+&
                   sic_orbitals%vtc(:,it)-atposc(:,ias2ia(jas),ias2is(jas))
            do ispn=1,nspinor
              s_wankmt(itp,ir,ias,ispn,j,ikloc)=s_wankmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt*s_func_val(x,s_wanlm(1,1,ispn,j))
              s_wvkmt(itp,ir,ias,ispn,j,ikloc)=s_wvkmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt*s_func_val(x,s_wvlm(1,1,ispn,j))
            enddo
          enddo !itp
        enddo !ir
      enddo !ias
! interstitial
      do ir=1,ngrtot
        x(:)=vgrc(:,ir)+sic_orbitals%vtc(:,it)-atposc(:,ias2ia(jas),ias2is(jas))
        do ispn=1,nspinor
          s_wankir(ir,ispn,j,ikloc)=s_wankir(ir,ispn,j,ikloc)+&
            s_func_val(x,s_wanlm(1,1,ispn,j))*expikt 
          s_wvkir(ir,ispn,j,ikloc)=s_wvkir(ir,ispn,j,ikloc)+&
            s_func_val(x,s_wvlm(1,1,ispn,j))*expikt 
        enddo
      enddo
    enddo !it
    zt1=zzero
    zt2=zzero
    do ispn=1,nspinor
      zt1=zt1+zfinp_(s_wankmt(1,1,1,1,j,ikloc),s_wankmt(1,1,1,1,j,ikloc),&
        s_wankir(1,1,j,ikloc),s_wankir(1,1,j,ikloc))
      zt2=zt2+zfinp_(s_wvkmt(1,1,1,1,j,ikloc),s_wankmt(1,1,1,1,j,ikloc),&
        s_wvkir(1,1,j,ikloc),s_wankir(1,1,j,ikloc))
    enddo
    write(*,*)"n =",n," <W_nk|W_nk> =",zt1," <(W*V)_nk|W_nk> = ",zt2
  enddo !j
enddo !ikloc 
deallocate(tp)
return
end
