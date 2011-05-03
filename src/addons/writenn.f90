subroutine writenn
use modmain
implicit none
integer i,ia,is,ias,ja,js,jas
real(8) dist,v1(3),v2(3),v3(3)
call getnghbr(-0.d0,nn_maxdist)
open(50,file='NGHBR.OUT',status='replace',form='formatted')
write(50,*)
write(50,'("Radius of cluster of nearest neighbours : ",G18.10)')nn_maxdist
write(50,*)
write(50,'("Atom list : ")')
write(50,'(" iatom (is ia ias)",15(" "),"pos(lat)",24(" "),"pos(Cart)")')
write(50,'(85("-"))')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(3X,A,T6," (",I2,1X,I2,1X,I3,")",T22,3F10.5,2X,3F10.5)') &
      trim(spsymb(is)),is,ia,ias,atposl(:,ia,is),atposc(:,ia,is)
  enddo
enddo
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  write(50,*)
  write(50,'("Cluster around ",A," (is ia ias : ",I2,1X,I2,1X,I3")")')trim(spsymb(is)),is,ia,ias
  write(50,'(" jatom (js ja jas)   D(a.u.)     D(A)",8(" "),"T",18(" "),"R(Cart)",12(" "),"shell")')
  write(50,'(88("-"))')
  do i=1,nnghbr(ias)
    jas=inghbr(1,i,ias)
    js=ias2is(jas)
    ja=ias2ia(jas)
    v1=atposl(:,ja,js)+inghbr(3:5,i,ias)
    v2=atposc(:,ja,js)+inghbr(3,i,ias)*avec(:,1)+&
      inghbr(4,i,ias)*avec(:,2)+inghbr(5,i,ias)*avec(:,3)
    v3=v2(:)-atposc(:,ia,is)
    dist=inghbr(2,i,ias)/1000000.d0
    if (inghbr(6,i,ias).le.nn_maxsh) then
      write(50,'(3X,A,T6," (",I2,1X,I2,1X,I3") ",2F10.5,2X,3I3,2X,3F10.5,2X,I4)') &
        trim(spsymb(js)),js,ja,jas,dist,dist*au2ang,inghbr(3:5,i,ias),v3,inghbr(6,i,ias)
    endif
  enddo
enddo    
close(50)
end
