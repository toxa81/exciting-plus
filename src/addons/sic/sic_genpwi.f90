subroutine sic_genpwi
use modmain
use mod_sic
implicit none
integer n,j,ik,ikloc,ir,itp,ispn,ig,ngvec1,lm
real(8) x(3)
complex(8) zt1(nspinor),pwtp,expikr
complex(8), allocatable :: wantp(:,:,:)
complex(8), allocatable :: wantp1(:,:)
complex(8), allocatable :: wvtp1(:,:)
complex(8), allocatable :: expigr(:)
real(8), allocatable :: stepf(:,:)
complex(8), external :: zdotu
!
!call timer_start(60,reset=.true.)
allocate(wantp1(sic_wantran%nwan,nspinor))
allocate(wvtp1(sic_wantran%nwan,nspinor))
allocate(stepf(s_ntp,s_nr))
call s_gen_stepf(stepf)
ngvec1=0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do ig=1,ngk(1,ik)
    ngvec1=max(ngvec1,igkig(ig,1,ikloc))
  enddo
enddo
allocate(expigr(ngvec1))
sic_wgk=zzero
sic_wvgk=zzero
do ir=1,s_nr
  do itp=1,s_ntp
    x(:)=s_spx(:,itp)*s_r(ir) !+pos of W_n <-- extra phase
    do ig=1,ngvec1
      expigr(ig)=exp(zi*dot_product(vgc(:,ig),x(:)))*s_tpw(itp)*s_rw(ir)*&
        stepf(itp,ir)/sqrt(omega)
    enddo
    wantp1=zzero
    wvtp1=zzero
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
        wantp1(j,ispn)=zdotu(lmmaxwan,s_wanlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)
        wvtp1(j,ispn)=zdotu(lmmaxwan,s_wvlm(1,ir,ispn,j),1,s_ylmf(1,itp),1)
      enddo
    enddo
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      expikr=exp(zi*dot_product(vkc(:,ik),x(:)))
      do ispn=1,nspinor
        do j=1,sic_wantran%nwan
          !n=sic_wantran%iwan(j)
          do ig=1,ngk(1,ik)
            sic_wgk(ig,j,ispn,ikloc)=sic_wgk(ig,j,ispn,ikloc)+&
              expigr(igkig(ig,1,ikloc))*expikr*dconjg(wantp1(j,ispn))
            sic_wvgk(ig,j,ispn,ikloc)=sic_wvgk(ig,j,ispn,ikloc)+&
              expigr(igkig(ig,1,ikloc))*expikr*dconjg(wvtp1(j,ispn))
          enddo
        enddo
      enddo
    enddo !ikloc
  enddo !itp
enddo !ir
deallocate(wantp1,wvtp1)
deallocate(stepf)
deallocate(expigr)
      
 
!allocate(wantp(s_ntp,s_nr,nspinor))
!do j=1,sic_wantran%nwan
!  n=sic_wantran%iwan(j)
!! convert functions to spherical coordinates
!  do ispn=1,nspinor
!    call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,&
!      s_wanlm(1,1,ispn,j),lmmaxwan,zzero,wantp(1,1,ispn),s_ntp)
!  enddo
!! not very efficient but saves memory
!  do ikloc=1,nkptloc
!    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!    do ig=1,ngk(1,ik)
!      zt1=zzero
!      do ir=1,s_nr
!        do itp=1,s_ntp
!          x(:)=s_spx(:,itp)*s_r(ir) !+pos of W_n
!          call s_get_pwval(ikloc,ig,x,pwtp)    
!          zt1(:)=zt1(:)+dconjg(wantp(itp,ir,:))*pwtp*s_tpw(itp)*s_rw(ir)
!        enddo
!      enddo !ir
!      sic_wgk(ig,j,:,ikloc)=zt1(:)
!    enddo !ig
!  enddo !ikloc
!enddo !j
!call timer_stop(60)
!write(*,*)"time : ",timer_get_value(60)
!deallocate(wantp)
return
end