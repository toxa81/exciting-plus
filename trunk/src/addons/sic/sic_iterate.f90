subroutine sic_iterate(fout)
use modmain
use mod_sic
use mod_nrkp
implicit none
integer, intent(in) :: fout
integer j,n,iter

do iter=1,sic_niter_umtrx
  if (iter.gt.1) then
    call sic_update_umtrx
  endif
  call wancnr_transform(sic_wan_umtrx)
  if (wproc) then
    write(fout,*)
    write(fout,'(80("="))')
    write(fout,'("SIC minimization")')
    write(fout,'(80("="))')
    write(fout,'("energies of Wannier functions")')
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      write(151,'("  n : ",I4,"    wann_ene : ",F12.6)')n,wann_ene(n)
    enddo
  endif
! generate Wannier functions and corresponding potential
  call sic_wan(fout)
! matrix elements
  call sic_genvme(fout)
enddo

return
end subroutine
