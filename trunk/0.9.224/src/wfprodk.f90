subroutine wfsvprodk(ngp,igpig,wfsvmt,wfsvit,wfnrmdev)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nstsv)
complex(8), intent(in) :: wfsvit(nmatmax,nstsv)
real(8), intent(out) :: wfnrmdev(nstsv*(nstsv+1)/2)

complex(8) norm
real(8) t1

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:) 

integer is,ia,ias,ig1,ig2,ist1,ist2,io1,io2,l,m,lm,i,j
integer iv2g(3)
real(8) v1(3),v2(3),tp2g(2),len2g
complex(8) sfac2g(natmtot)

allocate(a(ngp,nstsv))
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
do i=1,nstsv
  do ig2=1,ngp
    do ig1=1,ngp
      a(ig2,i)=a(ig2,i) + dconjg(wfsvit(ig1,i))*mit(ig1,ig2)
    enddo
  enddo
enddo

j=0
do ist1=1,nstsv
  do ist2=ist1,nstsv
    j=j+1
    norm=dcmplx(0.0,0.d0)
! interstitial contribution
    do ig2=1,ngp
      norm=norm+wfsvit(ig2,ist2)*a(ig2,ist1)
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
	        norm=norm+dconjg(wfsvmt(lm,io1,ias,ist1))* &
		  wfsvmt(lm,io2,ias,ist2)*urfprod(l,io1,io2,ias)
	      enddo !m
	    enddo
	  enddo
	enddo !l
      enddo !ia
    enddo !is
    t1=0.d0
    if (ist1.eq.ist2) t1=1.d0
    wfnrmdev(j)=abs(norm-t1)
  enddo !ist1 
enddo !ist2

deallocate(mit,a)

return
end

