subroutine write_response_f(iq,ivq0m,ifxc)
use modmain
use mod_addons_q
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: iq
integer, intent(in) :: ivq0m(3)
integer, intent(in) :: ifxc

real(8) fxca
real(8), allocatable :: func(:,:)
real(8) fsum
character*100 fname,qnm,qdir
character*10 c1,c2,c3,c4,c5
integer ie
real(8) t1,t2

fsum=0.d0
do ie=1,lr_nw-1
  t1=dimag(-1.d0/f_response(f_epsilon_eff,ie,ifxc))*dreal(lr_w(ie))
  t2=dimag(-1.d0/f_response(f_epsilon_eff,ie+1,ifxc))*dreal(lr_w(ie+1))
  fsum=fsum+0.5d0*dreal(lr_w(ie+1)-lr_w(ie))*(t1+t2)
enddo

fxca=fxca0+(ifxc-1)*fxca1

call getqdir(iq,ivq0m,qdir)
call getqname(ivq0m,qnm)
qnm=trim(qdir)//"/"//trim(qnm)
t1=sqrt(vqcnr(1,iq)**2+vqcnr(2,iq)**2+vqcnr(3,iq)**2)/au2ang
write(c2,'(F7.3)')fxca
write(c3,'(I8)')ngq(iq)
write(c4,'(F6.3)')t1
write(c5,'(F5.3)')lr_eta

if (lrtype.eq.0) then
  fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
    &"__A_"//trim(adjustl(c2))//"__.dat"
else
  fname=trim(qnm)//"_A"//c2//"_s"//c1//".dat"
endif
open(160,file=trim(fname),form='formatted',status='replace')
call write_response_f_header(160,iq,fxca)
write(160,'("#")')
write(160,'("# f-sum rule : ",G18.10)')fsum
write(160,'("# 2*Pi^2*rho_avg : ",G18.10)')2*pi*pi*(chgval/omega)
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
write(160,'("#  17:  loss_function                  ")')
write(160,'("#  18:  loss_function_scalar           ")')
write(160,'("#")')
allocate(func(20,lr_nw))
do ie=1,lr_nw
  func(1,ie)=dreal(lr_w(ie))*ha2ev
  func(2,ie)=-dreal(f_response(f_chi0,ie,ifxc))/ha2ev/(au2ang)**3
  func(3,ie)=-dimag(f_response(f_chi0,ie,ifxc))/ha2ev/(au2ang)**3
  func(4,ie)=-dreal(f_response(f_chi,ie,ifxc))/ha2ev/(au2ang)**3
  func(5,ie)=-dimag(f_response(f_chi,ie,ifxc))/ha2ev/(au2ang)**3
  func(6,ie)=-2.d0*dimag(f_response(f_chi,ie,ifxc))/ha2ev/(au2ang)**3
  func(7,ie)=-dreal(f_response(f_chi_scalar,ie,ifxc))/ha2ev/(au2ang)**3
  func(8,ie)=-dimag(f_response(f_chi_scalar,ie,ifxc))/ha2ev/(au2ang)**3
  func(9,ie)=dreal(f_response(f_epsilon_eff,ie,ifxc))
  func(10,ie)=dimag(f_response(f_epsilon_eff,ie,ifxc))
  func(11,ie)=dreal(f_response(f_epsilon_eff_scalar,ie,ifxc))
  func(12,ie)=dimag(f_response(f_epsilon_eff_scalar,ie,ifxc))
  func(13,ie)=dreal(f_response(f_sigma,ie,ifxc))*ha2ev
  func(14,ie)=dimag(f_response(f_sigma,ie,ifxc))*ha2ev
  func(15,ie)=dreal(f_response(f_sigma_scalar,ie,ifxc))*ha2ev
  func(16,ie)=dimag(f_response(f_sigma_scalar,ie,ifxc))*ha2ev
  func(17,ie)=-dimag(f_response(f_loss,ie,ifxc))
  func(18,ie)=-dimag(f_response(f_loss_scalar,ie,ifxc))
  func(19,ie)=dreal(f_response(f_epsilon_inv_GqGq,ie,ifxc))
  func(20,ie)=dimag(f_response(f_epsilon_inv_GqGq,ie,ifxc))
  write(160,'(20G14.6)')func(1:20,ie)
enddo
deallocate(func)
close(160)
if (wannier_chi0_chi) then
  if (lrtype.eq.0) then
    fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
      &"__A_"//trim(adjustl(c2))//"__wann__.dat"
  else
    fname=trim(qnm)//"_A"//c2//"_s"//c1//".dat"
  endif
  open(160,file=trim(fname),form='formatted',status='replace')
  call write_response_f_header(160,iq,fxca)
  write(160,'("#")')
  write(160,'("# Definition of columns")')
  write(160,'("#   1: energy            [eV]")')
  write(160,'("#   2: -Re chi0_wann(Gq,Gq)       [1/eV/A^3]   ")')
  write(160,'("#   3: -Im chi0_wann(Gq,Gq)       [1/eV/A^3]   ")')
  write(160,'("#   4: -Re chi_wann(Gq,Gq)        [1/eV/A^3]   ")')
  write(160,'("#   5: -Im chi_wann(Gq,Gq)        [1/eV/A^3]   ")')
  write(160,'("#   6:  Re epsilon_eff                         ")')
  write(160,'("#   7:  Im epsilon_eff                         ")')
  write(160,'("#   8:  Re sigma [eV]                          ")')
  write(160,'("#   9:  Im sigma [eV]                          ")')
  write(160,'("#  10:  loss_function                          ")')

  write(160,'("#")')
  allocate(func(10,lr_nw))
  do ie=1,lr_nw
    func(1,ie)=dreal(lr_w(ie))*ha2ev
    func(2,ie)=-dreal(f_response(f_chi0_wann,ie,ifxc))/ha2ev/(au2ang)**3
    func(3,ie)=-dimag(f_response(f_chi0_wann,ie,ifxc))/ha2ev/(au2ang)**3
    func(4,ie)=-dreal(f_response(f_chi_wann,ie,ifxc))/ha2ev/(au2ang)**3
    func(5,ie)=-dimag(f_response(f_chi_wann,ie,ifxc))/ha2ev/(au2ang)**3
    func(6,ie)=dreal(f_response(f_epsilon_eff_wann,ie,ifxc))
    func(7,ie)=dimag(f_response(f_epsilon_eff_wann,ie,ifxc))
    func(8,ie)=dreal(f_response(f_sigma_wann,ie,ifxc))*ha2ev
    func(9,ie)=dimag(f_response(f_sigma_wann,ie,ifxc))*ha2ev
    func(10,ie)=-dimag(f_response(f_loss_wann,ie,ifxc))
    write(160,'(10G14.6)')func(1:10,ie)
  enddo
  deallocate(func)
  close(160)
endif
return
end
