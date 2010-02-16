subroutine write_chi_header(fout,igq0,ivq0m,fxca)
use modmain
implicit none
integer, intent(in) :: fout
integer, intent(in) :: igq0
integer, intent(in) :: ivq0m(3)
real(8), intent(in) :: fxca

if (lrtype.eq.0) write(fout,'("# charge density response")')
if (lrtype.eq.1) write(fout,'("# magnetization density response")')
write(fout,'("#")')
write(fout,'("# Band interval (Ha) : ",2F8.2)')lr_e1,lr_e2
write(fout,'("#")')
write(fout,'("# k-mesh division                    : ",3I4)') &
  ngridk(1),ngridk(2),ngridk(3)
write(fout,'("# Energy mesh parameters             : ")')
write(fout,'("#   maximum energy [eV]              : ", F9.4)')maxomega
write(fout,'("#   energy step    [eV]              : ", F9.4)')domega
write(fout,'("#   eta            [eV]              : ", F9.4)')lr_eta
write(fout,'("# q-vector information               : ")')
write(fout,'("#   q-vector (mesh coord.)           : ",3I4)')ivq0m
write(fout,'("#   q-vector (lat. coord.)           : ",3F18.10)')vq0l
write(fout,'("#   q-vector (Cart. coord.) [a.u.]   : ",3F18.10)')vq0c
write(fout,'("#   q-vector length         [a.u.]   : ",3F18.10)') &
  sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
write(fout,'("#   q-vector (Cart. coord.) [1/A]    : ",3F18.10)')vq0c/au2ang
write(fout,'("#   q-vector length         [1/A]    : ",3F18.10)') &
  sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
write(fout,'("# G-vector information               : ")')
write(fout,'("#   G-shells                         : ",2I4)')gshme1,gshme2
write(fout,'("#   G-vectors                        : ",2I4)')gvecme1,gvecme2
write(fout,'("#   index of Gq vector               : ",I4)')igq0
if (lrtype.eq.0) then
  write(fout,'("#")')
  write(fout,'("# fxc type : ",I1)')fxctype
  write(fout,'("# fxc A : ",F8.4)')fxca
endif

return
end