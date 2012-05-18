subroutine write_chi_header(fout,iq,fxca)
use modmain
use mod_addons_q
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: fout
integer, intent(in) :: iq
real(8), intent(in) :: fxca
real(8) t1

if (lrtype.eq.0) write(fout,'("# charge density response")')
if (lrtype.eq.1) write(fout,'("# magnetization density response")')
write(fout,'("#")')
write(fout,'("# Included band interval (Ha)        : ",2F8.2)')&
  chi0_include_bands(1),chi0_include_bands(2)
write(fout,'("# Excluded band interval (Ha)        : ",2F8.2)')&
  chi0_exclude_bands(1),chi0_exclude_bands(2)
write(fout,'("# Approximate number of transitions  : ",I8)')nmegqblh(1)
write(fout,'("#")')
write(fout,'("# k-mesh division                    : ",3I4)') &
  ngridk(1),ngridk(2),ngridk(3)
write(fout,'("# Energy mesh parameters             : ")')
write(fout,'("#   energy interval [eV]             : ", 2F9.4)')lr_w0,lr_w1
write(fout,'("#   energy step     [eV]             : ", F9.4)')lr_dw
write(fout,'("#   eta             [eV]             : ", F9.4)')lr_eta
write(fout,'("# q-vector information               : ")')
write(fout,'("#   q-vector (mesh coord.)           : ",3I4)')vqm(:,iq)
write(fout,'("#   q-vector (lat. coord.)           : ",3F18.10)')vqlnr(:,iq)
write(fout,'("#   q-vector (Cart. coord.) [1/a.u.] : ",3F18.10)')vqcnr(:,iq)
t1=sqrt(vqcnr(1,iq)**2+vqcnr(2,iq)**2+vqcnr(3,iq)**2)
write(fout,'("#   q-vector length         [1/a.u.] : ",3F18.10)') t1
write(fout,'("#   q-vector (Cart. coord.) [1/A]    : ",3F18.10)') &
  vqcnr(:,iq)/au2ang
write(fout,'("#   q-vector length         [1/A]    : ",3F18.10)') t1/au2ang
write(fout,'("# G-vector information               : ")')
write(fout,'("#   number of G-vectors              : ",2I4)')ngvecme
write(fout,'("#   index of Gq vector               : ",I4)')ig0q(iq)
if (lrtype.eq.0) then
  write(fout,'("#")')
  write(fout,'("# fxc type : ",I1)')fxctype
  write(fout,'("# fxc A : ",F8.4)')fxca
endif

return
end