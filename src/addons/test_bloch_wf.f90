subroutine test_bloch_wf
use modmain
use mod_nrkp
implicit none
integer ikloc,ik,ispn,ig,ist1,ist2,ias,is,ic,io1,io2,l,m,lm,ir
real(8), allocatable :: maxdev(:)
complex(8) zprod
complex(8), allocatable :: wfir1(:,:)
complex(8), allocatable :: wfir2(:,:)
!
call init0
call init1
call readstate
call gencore
call linengy
call genapwfr
call genlofr
call getufr
call genufrp
wproc=mpi_grid_root()
call genwfnr(6,.false.)
allocate(wfir1(ngrtot,nspinor))
allocate(wfir2(ngrtot,nspinor))
allocate(maxdev(nkptnr))
maxdev=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do ist1=1,nstsv
    wfir1=zzero
    do ispn=1,nspinor
      do ig=1,ngknr(ikloc)
        wfir1(igfft(igkignr(ig,ikloc)),ispn)=wfsvitnrloc(ig,ispn,ist1,ikloc)/sqrt(omega)
      enddo
      call zfftifc(3,ngrid,1,wfir1(1,ispn))
    enddo
    do ist2=1,nstsv
      wfir2=zzero
      do ispn=1,nspinor
        do ig=1,ngknr(ikloc)
          wfir2(igfft(igkignr(ig,ikloc)),ispn)=wfsvitnrloc(ig,ispn,ist2,ikloc)/sqrt(omega)
        enddo
        call zfftifc(3,ngrid,1,wfir2(1,ispn))
      enddo
      zprod=zzero
      do ispn=1,nspinor
! muffin-tin part
        do ias=1,natmtot
          is=ias2is(ias)
          ic=ias2ic(ias)
          do l=0,lmaxapw
            do io1=1,nufr(l,is); do io2=1,nufr(l,is)
              do m=-l,l
                lm=idxlm(l,m)
                zprod=zprod+dconjg(wfsvmtnrloc(lm,io1,ias,ispn,ist1,ikloc))*&
                  wfsvmtnrloc(lm,io2,ias,ispn,ist2,ikloc)*&
                  ufrp(l,io1,io2,ic)
              enddo
            enddo; enddo
          enddo !l
        enddo !ias
! interstitial contribution
        do ir=1,ngrtot
          zprod=zprod+cfunir(ir)*dconjg(wfir1(ir,ispn))*wfir2(ir,ispn)*omega/dble(ngrtot)
        enddo
      enddo !ispn
      if (ist1.eq.ist2) zprod=zprod-zone
      maxdev(ik)=max(maxdev(ik),abs(zprod))
    enddo !ist2
  enddo !ist1
enddo !ikloc
call mpi_grid_reduce(maxdev(1),nkptnr)
if (mpi_grid_root()) then
  do ik=1,nkptnr
    write(*,'("k-point : ",I4,"    maximum deviation from norm : ",G18.10)')&
      ik,maxdev(ik)
  enddo
endif
deallocate(wfir1,wfir2)
return
end