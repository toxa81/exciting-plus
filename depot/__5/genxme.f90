subroutine genxme
use modmain
use mod_sic
implicit none
integer lm1,lm2,j1,j2,n1,n2,ispn,ir
complex(8) zt1,zt2
complex(8), external :: gauntyry

do j1=1,sic_wantran%nwan
  n1=sic_wantran%iwan(j1)
  do j2=1,sic_wantran%nwan
    n2=sic_wantran%iwan(j2)
    zt2=zzero
    do lm1=1,lmmaxwan
      do lm2=1,lmmaxwan
        zt1=gauntyry(lm2l(lm1),1,lm2l(lm2),lm2m(lm1),-1,lm2m(lm2))
        if (abs(zt1).gt.1d-10) then
          do ispn=1,nspinor
            do ir=1,s_nr
              zt2=zt2+dconjg(s_wlm(lm1,ir,ispn,j1))*(-2.d0*sqrt(pi/3)*s_r(ir))*&
                s_wlm(lm2,ir,ispn,j2)*zt1*s_rw(ir)
            enddo
          enddo !ispn
        endif
      enddo
    enddo
    write(*,*)"n1,n2=",n1,n2,"  y=",zt2
  enddo
enddo

return
end


