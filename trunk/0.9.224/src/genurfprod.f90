subroutine genurfprod
use modmain
implicit none

real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
integer is,ia,ias,l,io1,io2,ir

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do l=0,lmaxvr
      do io1=1,nrfmax
        do io2=1,nrfmax
          do ir=1,nrmt(is)
            fr(ir)=urf(ir,l,io1,ias)*urf(ir,l,io2,ias)*(spr(ir,is)**2)                                                        
          enddo
          call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
          urfprod(l,io1,io2,ias)=gr(nrmt(is))
        enddo
      enddo
    enddo 
  enddo
enddo

return
end

