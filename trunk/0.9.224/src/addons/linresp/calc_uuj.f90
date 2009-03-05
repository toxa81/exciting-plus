subroutine calc_uuj(uuj,lmaxexp,gq0)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
real(8), intent(in) :: gq0(ngvecme)
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
! local variables
integer ia,is,ias,ig,l1,l2,l3,io1,io2,ir
real(8), allocatable :: jlgq0r(:,:,:,:)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) t1,jl(0:lmaxexp)

allocate(jlgq0r(nrmtmax,0:lmaxexp,nspecies,ngvecme))

! generate Bessel functions j_l(|G+q'|x)
do ig=1,ngvecme
  do is=1,nspecies
    do ir=1,nrmt(is)
      t1=gq0(ig)*spr(ir,is)
      call sbessel(lmaxexp,t1,jl)
      jlgq0r(ir,:,is,ig)=jl(:)
    enddo
  enddo
enddo

uuj=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvecme
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do l3=0,lmaxexp
            do io1=1,nrfmax
              do io2=1,nrfmax
                do ir=1,nrmt(is)
                  fr(ir)=urf(ir,l1,io1,ias)*urf(ir,l2,io2,ias)*jlgq0r(ir,l3,is,ig)*(spr(ir,is)**2)
                enddo
                call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
                uuj(l1,l2,l3,io1,io2,ias,ig)=gr(nrmt(is))
              enddo !io2
            enddo !io1
          enddo !l3
        enddo !l2
      enddo !l1   
    enddo !ig
  enddo !ia
enddo !is

deallocate(jlgq0r)
return
end   
