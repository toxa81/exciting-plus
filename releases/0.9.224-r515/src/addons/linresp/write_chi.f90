subroutine write_chi(igq0,ivq0m,chi_,epsilon_,fxca)
use modmain
implicit none
integer, intent(in) :: igq0
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: chi_(7,nepts)
complex(8), intent(in) :: epsilon_(6,nepts)
real(8), intent(in) :: fxca

real(8), allocatable :: func(:,:)
character*100 fname,qnm
character*10 c1,c2,c3,c4,c5
integer ie,i
complex(8) z1

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
write(c2,'(F6.3)')fxca
write(c3,'(I8)')ngvecchi
write(c4,'(F6.3)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
write(c5,'(F5.3)')lr_eta

if (lrtype.eq.0) then
  fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
    "__A_"//trim(adjustl(c2))//"__.dat"
else
  fname=trim(qnm)//"_A"//c2//"_s"//c1//".dat"
endif
open(160,file=trim(fname),form='formatted',status='replace')
if (lrtype.eq.0) write(160,'("# charge density response")')
if (lrtype.eq.1) write(160,'("# magnetization density response")')
write(160,'("#")')
!write(160,'("# Band interval (Ha) : ",2F8.2)')lr_e1,lr_e2
!write(160,'("#")')
write(160,'("# k-mesh division                    : ",3I4)')ngridk(1),ngridk(2),ngridk(3)
write(160,'("# Energy mesh parameters             : ")')
write(160,'("#   maximum energy [eV]              : ", F9.4)')maxomega
write(160,'("#   energy step    [eV]              : ", F9.4)')domega
write(160,'("#   eta            [eV]              : ", F9.4)')lr_eta
write(160,'("# q-vector information               : ")')
write(160,'("#   q-vector (mesh coord.)           : ",3I4)')ivq0m
write(160,'("#   q-vector (lat. coord.)           : ",3F18.10)')vq0l
write(160,'("#   q-vector (Cart. coord.) [a.u.]   : ",3F18.10)')vq0c
write(160,'("#   q-vector length         [a.u.]   : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
write(160,'("#   q-vector (Cart. coord.) [1/A]    : ",3F18.10)')vq0c/au2ang
write(160,'("#   q-vector length         [1/A]    : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
write(160,'("# G-vector information               : ")')
write(160,'("#   G-shells                         : ",2I4)')gshchi1,gshchi2
write(160,'("#   G-vectors                        : ",2I4)')gvecchi1,gvecchi2
write(160,'("#   index of Gq vector               : ",I4)')igq0
if (lrtype.eq.0) then
  write(160,'("#")')
  write(160,'("# fxc type : ",I1)')fxctype
  write(160,'("# fxc A : ",F8.4)')fxca
endif
write(160,'("#")')
write(160,'("# Definition of columns")')
write(160,'("#   1: energy            [eV]")')
write(160,'("#   2: -Re chi0(Gq,Gq)   [1/eV/A^3]    ")')
write(160,'("#   3: -Im chi0(Gq,Gq)   [1/eV/A^3]    ")')
write(160,'("#   4: -Re chi(Gq,Gq)    [1/eV/A^3]    ")')
write(160,'("#   5: -Im chi(Gq,Gq)    [1/eV/A^3]    ")')
write(160,'("#   6:  S(q,w)           [1/eV/A^3]    ")')
write(160,'("#   7: -Re chi_scalar    [1/eV/A^3]    ")')
write(160,'("#   8: -Im chi_scalar    [1/eV/A^3]    ")')
write(160,'("#   9:  Re epsilon_eff                 ")')
write(160,'("#  10:  Im epsilon_eff                 ")')
write(160,'("#  11:  Re epsilon_eff_scalar          ")')
write(160,'("#  12:  Im epsilon_eff_scalar          ")')
write(160,'("#  13:  Re sigma [eV]                  ")')
write(160,'("#  14:  Im sigma [eV]                  ")')
write(160,'("#  15:  Re sigma_scalar [eV]           ")')
write(160,'("#  16:  Im sigma_scalar [eV]           ")')
write(160,'("#  17: -Re chi0_wf_full(Gq,Gq)   [1/eV/A^3]    ")')
write(160,'("#  18: -Im chi0_wf_full(Gq,Gq)   [1/eV/A^3]    ")')
write(160,'("#  19: -Re chi0_wf(Gq,Gq)        [1/eV/A^3]    ")')
write(160,'("#  20: -Im chi0_wf(Gq,Gq)        [1/eV/A^3]    ")')
write(160,'("#  21: -Re chi_wf(Gq,Gq)         [1/eV/A^3]    ")')
write(160,'("#  22: -Im chi_wf(Gq,Gq)         [1/eV/A^3]    ")')
write(160,'("#  23:  Re epsilon_eff_wf              ")')
write(160,'("#  24:  Im epsilon_eff_wf              ")')
write(160,'("#  25:  loss_function                  ")')
write(160,'("#  26:  loss_function_wf               ")')
write(160,'("#")')
allocate(func(26,nepts))
do ie=1,nepts
  func(1,ie)=dreal(lr_w(ie))*ha2ev
  func(2,ie)=-dreal(chi_(1,ie))/ha2ev/(au2ang)**3
  func(3,ie)=-dimag(chi_(1,ie))/ha2ev/(au2ang)**3
  func(4,ie)=-dreal(chi_(4,ie))/ha2ev/(au2ang)**3
  func(5,ie)=-dimag(chi_(4,ie))/ha2ev/(au2ang)**3
  func(6,ie)=-2.d0*dimag(chi_(4,ie))/ha2ev/(au2ang)**3
  func(7,ie)=-dreal(chi_(2,ie))/ha2ev/(au2ang)**3
  func(8,ie)=-dimag(chi_(2,ie))/ha2ev/(au2ang)**3
  func(9,ie)=dreal(epsilon_(4,ie))
  func(10,ie)=dimag(epsilon_(4,ie))
  func(11,ie)=dreal(epsilon_(5,ie))
  func(12,ie)=dimag(epsilon_(5,ie))
  z1=zi*dreal(lr_w(ie))*(zone-epsilon_(4,ie))/fourpi
  func(13,ie)=dreal(z1)*ha2ev
  func(14,ie)=dimag(z1)*ha2ev
  z1=zi*dreal(lr_w(ie))*(zone-epsilon_(5,ie))/fourpi
  func(15,ie)=dreal(z1)*ha2ev
  func(16,ie)=dimag(z1)*ha2ev
  func(17,ie)=-dreal(chi_(5,ie))/ha2ev/(au2ang)**3
  func(18,ie)=-dimag(chi_(5,ie))/ha2ev/(au2ang)**3
  func(19,ie)=-dreal(chi_(6,ie))/ha2ev/(au2ang)**3
  func(20,ie)=-dimag(chi_(6,ie))/ha2ev/(au2ang)**3
  func(21,ie)=-dreal(chi_(7,ie))/ha2ev/(au2ang)**3
  func(22,ie)=-dimag(chi_(7,ie))/ha2ev/(au2ang)**3
  func(23,ie)=dreal(epsilon_(6,ie))
  func(24,ie)=dimag(epsilon_(6,ie))
  func(25,ie)=-imag(epsilon_(7,ie))
  func(26,ie)=-imag(epsilon_(8,ie))
  write(160,'(26G14.6)')func(1:26,ie)
enddo
deallocate(func)
close(160)
! eigen-values of denominator matrix
!fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
!  "__A_"//trim(adjustl(c2))//"_epseval__.dat"
!open(160,file=trim(fname),form='formatted',status='replace')
!do i=1,ngvecchi
!  do ie=1,nepts
!    write(160,'(4G18.10)')dreal(lr_w(ie))*ha2ev,abs(1.d0/lmbd(i,ie)),dreal(lmbd(i,ie)),dimag(lmbd(i,ie))  
!  enddo
!  write(160,'(" ")')
!enddo
!close(160)
return
end
