subroutine genwffvwann(ik,ngp,igpig,wffvmt,evecfv,wffvwann)
use modmain
implicit none
integer, intent(in) :: ik
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: wffvmt(nstfv,nrfmax,lmmaxvr,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(out) :: wffvwann(nstfv,wann_nmax,wann_nspin)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:) 
integer is,ia,ias,ig1,ig2,istfv,io1,io2,l,m,lm,i,ispn,n
integer iv2g(3)
real(8) v1(3),v2(3),tp2g(2),len2g
complex(8) sfac2g(natmtot),z1




wffvwann=dcmplx(0.d0,0.d0)

allocate(a(ngp,nstfv))
allocate(mit(ngp,ngp))

mit=dcmplx(0.d0,0.d0)
do ig1=1,ngp
  do ig2=1,ngp
! G1-G2
    iv2g(:)=ivg(:,igpig(ig1))-ivg(:,igpig(ig2))
    if (sum(abs(iv2g)).eq.0) mit(ig1,ig2)=dcmplx(1.d0,0.d0)
    v2(:)=1.d0*iv2g(:)
    call r3mv(bvec,v2,v1)
    call sphcrd(v1,len2g,tp2g)
    call gensfacgp(1,v1,1,sfac2g)
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        if (len2g.lt.1d-8) then
          mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac2g(ias))*(rmt(is)**3)/3.d0
        else
          mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac2g(ias)) * &
            (-(rmt(is)/len2g**2)*cos(len2g*rmt(is))+(1/len2g**3)*sin(len2g*rmt(is)))
        endif
      enddo !ia
    enddo !is
  enddo
enddo

a=dcmplx(0.d0,0.d0)
do i=1,nstfv
  do ig2=1,ngp
    do ig1=1,ngp
      a(ig2,i)=a(ig2,i)+dconjg(evecfv(ig1,i))*mit(ig1,ig2)
    enddo
  enddo
enddo

do ispn=1,wann_nspin
  do n=1,nwann(ispn)
    do istfv=1,nstfv
      do ig2=1,ngp
! interstitial contribution       
        wffvwann(istfv,n,ispn)=wffvwann(istfv,n,ispn) + &
          wann_unkit(ig2,n,ispn,ik)*a(ig2,istfv)
      enddo
! muffin-tin contribution
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          do l=0,lmaxvr
            do io1=1,nrfmax
              do io2=1,nrfmax
                do m=-l,l
                  lm=idxlm(l,m)
                  wffvwann(istfv,n,ispn)=wffvwann(istfv,n,ispn) + &
                    dconjg(wffvmt(istfv,io1,lm,ias)) * &
                    wann_unkmt(lm,io2,ias,n,ispn,ik)*urfprod(l,io1,io2,ias)
                enddo !m
              enddo
            enddo
          enddo !l
        enddo !ia
      enddo !is
    enddo !istfv
  enddo !n
enddo !ispn

deallocate(mit,a)

return
end

