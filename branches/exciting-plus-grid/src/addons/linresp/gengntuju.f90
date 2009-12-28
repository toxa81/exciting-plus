subroutine gengntuju(lmaxexp,ngntujumax,ngntuju,igntuju,gntuju)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
integer, intent(in) :: ngntujumax
integer, intent(out) :: ngntuju(natmtot,ngvecmeloc)
integer, intent(out) :: igntuju(4,ngntujumax,natmtot,ngvecmeloc)
complex(8), intent(out) :: gntuju(ngntujumax,natmtot,ngvecmeloc)


integer ig,is,ir,n,ias,i,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3
integer igloc,igloc1,ngvecmeloc1
! for parallel
integer cart_size,cart_rank,cart_group,ierr
integer tmp_comm,tmp_group,tmp_size
integer, allocatable :: ranks(:)
integer idx0,bs
complex(8) zt1
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8), allocatable :: jl(:,:)
real(8), allocatable :: uju(:,:,:,:,:)
real(8), external :: gaunt

! use first direction in the grid (direction of k-points) to compute
!   MT integrals in parallel
ngvecmeloc1=mpi_grid_map(ngvecmeloc,dim_k)

allocate(jl(nrmtmax,0:lmaxexp))
allocate(uju(0:lmaxexp,0:lmaxvr,0:lmaxvr,nrfmax,nrfmax))

igntuju=0
ngntuju=0
gntuju=dcmplx(0.d0,0.d0)
! loop over G-vectors
do igloc1=1,ngvecmeloc1
  igloc=mpi_grid_map(ngvecmeloc,dim_k,loc=igloc1)
  ig=mpi_grid_map(ngvecme,dim_g,loc=igloc)
! loop over atoms
  do ias=1,natmtot
    is=ias2is(ias)
! generate Bessel functions j_l(|G+q'|x)
    do ir=1,nrmt(is)
      call sbessel(lmaxexp,lr_gq0(ig)*spr(ir,is),jl(ir,:))
    enddo
! compute radial integrals <u_{l1,io1} | j_{l3}(|G+q|x) | u_{l2,io2}>
    do l3=0,lmaxexp
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do io1=1,nrfmax
            do io2=1,nrfmax
              do ir=1,nrmt(is)
                fr(ir)=urf(ir,l1,io1,ias)*urf(ir,l2,io2,ias)*&
                  jl(ir,l3)*(spr(ir,is)**2)
              enddo !ir
              call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
              uju(l3,l1,l2,io1,io2)=gr(nrmt(is))
            enddo !io2
          enddo !io1
        enddo !l2
      enddo !l1
    enddo !l3
! compute muffin-tin integrals
    do io1=1,nrfmax
      do io2=1,nrfmax
        do l1=0,lmaxvr
        do m1=-l1,l1 
          lm1=idxlm(l1,m1)
          do l2=0,lmaxvr
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
!  1) sfacgq0 and ylmgq0 are generated for exp^{+i(G+q)x}
!     expansion of a plane-wave: 
!       exp^{+igx}=4\pi \sum_{l_3 m_3} i^{l_3} j_{l_3}(gr)Y_{l_3 m_3}^{*}(\hat g)Y_{l_3 m_3}(\hat r)
!     but we need exp^{-i(G+q)x}, so expansion terms will be conjugated
!  2) angular part of integral:
!     <Y_{l_1 m_1} | e^{-i{G+x}x} | Y_{l_2 m_2}> =
!       = \int d \Omega Y_{l_1 m_1}^{*}Y_{l_3 m_3}^{*} Y_{l_2 m_2} = gaunt coeff, which is real
!     so we can conjugate the integral:
!     \int d \Omega Y_{l_1 m_1} Y_{l_3 m_3} Y_{l_2 m_2}^{*} = gaunt(lm2,lm1,lm3)
!  3) we can sum over lm3 index of a plane-wave expansion           
            zt1=zzero
            do l3=0,lmaxexp
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              zt1=zt1+gaunt(l2,l1,l3,m2,m1,m3)*uju(l3,l1,l2,io1,io2)*&
                lr_ylmgq0(lm3,ig)*dconjg(zi**l3)*fourpi*dconjg(lr_sfacgq0(ig,ias))
            enddo
            enddo
            if (abs(zt1).gt.1d-16) then
              ngntuju(ias,igloc)=ngntuju(ias,igloc)+1
              n=ngntuju(ias,igloc)
              gntuju(n,ias,igloc)=zt1
              igntuju(1,n,ias,igloc)=lm1
              igntuju(2,n,ias,igloc)=lm2
              igntuju(3,n,ias,igloc)=io1
              igntuju(4,n,ias,igloc)=io2
            endif
          enddo !m2
          enddo !l2
        enddo !m1
        enddo !l1
      enddo !io2
    enddo !io1
  enddo !ias
enddo !igloc1
! syncronize all values along auxiliary k-direction
do igloc=1,ngvecmeloc
  call mpi_grid_reduce(gntuju(1,1,igloc),ngntujumax*natmtot,dims=(/dim_k/),&
    all=.true.)
  call mpi_grid_barrier(dims=(/dim_k/))
enddo
deallocate(jl)
deallocate(uju)
return
end
