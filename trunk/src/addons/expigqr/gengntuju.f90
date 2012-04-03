subroutine gengntuju(iq,lmaxexp)
use modmain
use mod_addons_q
use mod_expigqr
use mod_util
implicit none
! arguments
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp

integer ig,is,ir,n,ias,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,ic
! for parallel
integer igloc,ngvecmeloc
integer lmmaxexp
complex(8) zt1
real(8) fr(nrmtmax)
real(8), allocatable :: jl(:,:)
real(8), allocatable :: uju(:,:,:,:,:)
real(8), allocatable :: gnt(:,:,:)
real(8), external :: gaunt
real(8), external :: rfinteg
complex(8), allocatable :: zm(:,:,:)
integer, parameter :: ngvb=2
integer i

lmmaxexp=(lmaxexp+1)**2
allocate(gnt(lmmaxexp,lmmaxvr,lmmaxvr))
do l1=0,lmaxvr
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    do l2=0,lmaxvr
      do m2=-l2,l2
        lm2=idxlm(l2,m2)
        do l3=0,lmaxexp
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            gnt(lm3,lm2,lm1)=gaunt(l2,l1,l3,m2,m1,m3)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
allocate(jl(nrmtmax,0:lmaxexp))
allocate(uju(0:lmaxexp,0:lmaxvr,0:lmaxvr,nufrmax,nufrmax))
allocate(zm(0:lmaxexp,lmmaxvr,lmmaxvr))
igntuju=0
ngntuju=0
gntuju=zzero
! loop over G-vectors
ngvecmeloc=mpi_grid_map(ngvecme,dim_k)
do igloc=1,ngvecmeloc
  ig=mpi_grid_map(ngvecme,dim_k,loc=igloc)
! precompute atom-independent array
  do l1=0,lmaxvr
    do m1=-l1,l1 
      lm1=idxlm(l1,m1)
      do l2=0,lmaxvr
        do m2=-l2,l2
          lm2=idxlm(l2,m2)
          do l3=0,lmaxexp
            zt1=zzero
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              zt1=zt1+gnt(lm3,lm2,lm1)*ylmgq(lm3,ig)*dconjg(zi**l3)
            enddo !m3
            zm(l3,lm2,lm1)=zt1*fourpi
          enddo !l3
        enddo
      enddo
    enddo
  enddo
! loop over atom classes
  do ic=1,natmcls
    ias=ic2ias(ic)
    is=ias2is(ias)
! generate Bessel functions j_l(|G+q'|x)
    do ir=1,nrmt(is)
      call sbessel(lmaxexp,gq(ig,iq)*spr(ir,is),jl(ir,:))
    enddo
! compute radial integrals <u_{l1,io1} | j_{l3}(|G+q|x) | u_{l2,io2}>
    do l3=0,lmaxexp
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do io1=1,nufr(l1,is)
            do io2=1,nufr(l2,is)
              do ir=1,nrmt(is)
                fr(ir)=ufr(ir,l1,io1,ic)*ufr(ir,l2,io2,ic)*jl(ir,l3)
              enddo !ir
              uju(l3,l1,l2,io1,io2)=rintegrate(nrmt(is),spr(1,is),fr)
            enddo !io2
          enddo !io1
        enddo !l2
      enddo !l1
    enddo !l3
! compute muffin-tin integrals
!  1) sfacgq and ylmgq are generated for exp^{+i(G+q)x}
!     expansion of a plane-wave: 
!       exp^{+igx}=4\pi \sum_{l_3 m_3} i^{l_3} j_{l_3}(gr)Y_{l_3 m_3}^{*}(\hat g)Y_{l_3 m_3}(\hat r)
!     but we need exp^{-i(G+q)x}, so expansion terms will be conjugated
!  2) angular part of integral:
!     <Y_{l_1 m_1} | e^{-i{G+x}x} | Y_{l_2 m_2}> =
!       = \int d \Omega Y_{l_1 m_1}^{*}Y_{l_3 m_3}^{*} Y_{l_2 m_2} = gaunt coeff, which is real
!     so we can conjugate the integral:
!     \int d \Omega Y_{l_1 m_1} Y_{l_3 m_3} Y_{l_2 m_2}^{*} = gaunt(lm2,lm1,lm3)
!  2*) gnt array has different order of indices: gnt(lm3,lm2,lm1)=gaunt(lm2,lm1,lm3)
!  3) we can sum over lm3 index of a plane-wave expansion
!  3*) we can sum over m3 index and get atom-idependent array
!  4) structure factor sfacgq of a (G+q) plane wave is taken into 
!     account in genmegqblh subroutine; this allows to keep radial integrals
!     for atom classes only (not for all atoms)
    do l1=0,lmaxvr
      do m1=-l1,l1 
        lm1=idxlm(l1,m1)
        do l2=0,lmaxvr
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            do io1=1,nufr(l1,is)
              do io2=1,nufr(l2,is)
                zt1=zzero
                do l3=0,lmaxexp
                  !do m3=-l3,l3
                  !  lm3=idxlm(l3,m3)
                  !  zt1=zt1+gnt(lm3,lm2,lm1)*uju(l3,l1,l2,io1,io2)*&
                  !    ylmgq(lm3,ig)*dconjg(zi**l3)
                  !enddo !m3
                  zt1=zt1+zm(l3,lm2,lm1)*uju(l3,l1,l2,io1,io2)
                enddo !l3
                !zt1=zt1*fourpi
                if (abs(zt1).gt.1d-12) then
                  ngntuju(ic,ig)=ngntuju(ic,ig)+1
                  n=ngntuju(ic,ig)
                  gntuju(n,ic,ig)=zt1
                  igntuju(1,n,ic,ig)=lm1+(io1-1)*lmmaxvr
                  igntuju(2,n,ic,ig)=lm2+(io2-1)*lmmaxvr
                endif
              enddo !io2
            enddo !io1
          enddo !m2
        enddo !l2
      enddo !m1
    enddo !l1
  enddo !ic
enddo !ig
! syncronize all values along auxiliary k-direction
!call mpi_grid_reduce(gntuju(1,1,1),ngntujumax*natmcls*ngvecme,dims=(/dim_k/),all=.true.)
!call mpi_grid_barrier(dims=(/dim_k/))
!call mpi_grid_reduce(igntuju(1,1,1,1),2*ngntujumax*natmcls*ngvecme,dims=(/dim_k/),all=.true.)
!call mpi_grid_barrier(dims=(/dim_k/))
!call mpi_grid_reduce(ngntuju(1,1),natmcls*ngvecme,dims=(/dim_k/),all=.true.)    
!call mpi_grid_barrier(dims=(/dim_k/))
!do ig=1,ngvecme
!  call mpi_grid_reduce(gntuju(1,1,ig),ngntujumax*natmcls,dims=(/dim_k/),&
!    all=.true.)
!  call mpi_grid_reduce(igntuju(1,1,1,ig),2*ngntujumax*natmcls,dims=(/dim_k/),&
!    all=.true.)
!  call mpi_grid_reduce(ngntuju(1,ig),natmcls,dims=(/dim_k/),all=.true.)    
!  call mpi_grid_barrier(dims=(/dim_k/))
!enddo

! synchronize blocks of G-vectors arrays
i=ngvecme/ngvb
do ig=1,i
  call mpi_grid_reduce(gntuju(1,1,(ig-1)*ngvb+1),ngvb*ngntujumax*natmcls,&
    dims=(/dim_k/),all=.true.)
  call mpi_grid_reduce(igntuju(1,1,1,(ig-1)*ngvb+1),ngvb*2*ngntujumax*natmcls,&
    dims=(/dim_k/),all=.true.)
  call mpi_grid_reduce(ngntuju(1,(ig-1)*ngvb+1),ngvb*natmcls,dims=(/dim_k/),&
    all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo
do ig=i*ngvb+1,ngvecme
  call mpi_grid_reduce(gntuju(1,1,ig),ngntujumax*natmcls,dims=(/dim_k/),&
    all=.true.)
  call mpi_grid_reduce(igntuju(1,1,1,ig),2*ngntujumax*natmcls,dims=(/dim_k/),&
    all=.true.)
  call mpi_grid_reduce(ngntuju(1,ig),natmcls,dims=(/dim_k/),all=.true.)    
  call mpi_grid_barrier(dims=(/dim_k/))
enddo

deallocate(jl)
deallocate(uju)
deallocate(gnt)
deallocate(zm)
return
end
