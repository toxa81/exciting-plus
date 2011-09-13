subroutine printmegqblh2(iq)
use modmain
implicit none
integer, intent(in) :: iq
integer i,ig,ikloc,ik
integer ist1,ist2
character*100 fname
complex(8) t1me
write(fname,'("wst_MEGQBLH_iq_",I4.4,".OUT")')iq
!if (iproc.eq.0) then
!  open(200,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
!  do ikloc=1,nkptnrloc
!    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!    write(200,'("ikloc : ",I4,"  idxkq : ",2I4)')ikloc,idxkq(:,ik)
!    do ist2=1,nstsv
!      write(200,'("  ist2 : ",I4)')ist2
!      do ist1=1,nstsv
!        do i=1,nmegqblh(ikloc)
!          if (bmegqblh(1,i,ikloc).eq.ist1.and.bmegqblh(2,i,ikloc).eq.ist2) then
!            write(200,'("    ist1 : ",I4)')ist1
!            do ig=1,ngvecme
!              write(200,'("      ig : ",I4,"   ",3G18.10)')&
!                ig,dreal(megqblh(i,ig,ikloc)),-dimag(megqblh(i,ig,ikloc)),&
!                abs(megqblh(i,ig,ikloc))**2
!            enddo
!          endif
!        enddo
!      enddo
!    enddo
!  enddo !ikloc
!  close(200)
!endif

if (iproc.eq.0) then
  open(200,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do ist2=1,nstsv
      do ist1=1,nstsv
        do i=1,nmegqblh(ikloc)
          if (bmegqblh(1,i,ikloc).eq.ist1.and.bmegqblh(2,i,ikloc).eq.ist2) then
            do ig=1,ngvecme
              t1me = t1me + megqblh(i,ig,ikloc);
            enddo
            write(200,'4I6,"    ",3G18.10')ikloc,idxkq(:,ik),&
              ist2,ist1,dreal(t1me),-dimag(t1me),abs(t1me)**2
          endif
        enddo
      enddo
    enddo
  enddo !ikloc
  close(200)
endif
return
end
