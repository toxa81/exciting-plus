subroutine genufrp
use modmain
implicit none
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
integer is,ias,l,io1,io2,ir,ic
ufrp=0.d0
do ic=1,natmcls
  ias=ic2ias(ic)
  is=ias2is(ias)
  do l=0,lmaxapw
    do io1=1,nufr(l,is)
      do io2=1,nufr(l,is)
        do ir=1,nrmt(is)
          fr(ir)=ufr(ir,l,io1,ic)*ufr(ir,l,io2,ic)*(spr(ir,is)**2)                                                        
        enddo
        call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
        ufrp(l,io1,io2,ic)=gr(nrmt(is))
      enddo
    enddo
  enddo 
enddo
return
end

