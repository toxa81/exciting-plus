
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepvnl(vnlcv,vnlvv)
use modmain
implicit none
! arguments
complex(8), intent(out) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(out) :: vnlvv(nstsv,nstsv,nkpt)
! local variables
integer ik
do ik=1,nkpt
  write(*,'("Info(oepvnl): ",I6," of ",I6," k-points")') ik,nkpt
  call oepvnlk(ik,vnlcv(:,:,:,ik),vnlvv(:,:,ik))
end do
return
end subroutine

