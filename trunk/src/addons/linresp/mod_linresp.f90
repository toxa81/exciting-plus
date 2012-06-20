module mod_linresp
implicit none

! dimension for q-vectors
integer, parameter :: dim_q=2

! type of linear response calculation
!   0 : charge response
!   1 : magnetic response (not yet implemented)
integer lrtype
data lrtype/0/

! temporary array used in genchi0blh
complex(8), allocatable :: megqblh2(:,:)

! interval of bands to include in chi0
real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/

! interval of bands to exclude from chi0
real(8) chi0_exclude_bands(2)
data chi0_exclude_bands/100.1d0,-100.1d0/

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)

! number of energy-mesh points
integer lr_nw
data lr_nw/201/
! first energy point (Ha)
real(8) lr_w0
data lr_w0/0.d0/
! last energy point (Ha)
real(8) lr_w1
data lr_w1/1.d0/
! energy step
real(8) lr_dw
! energy mesh
complex(8), allocatable :: lr_w(:)
! broadening parameter (Ha)
real(8) lr_eta
data lr_eta/0.01d0/
! inverse temperature for the matsubara frequency in eV^-1
real(8) lr_beta
data lr_beta/30.d0/
! .true. if imaginary frequency mesh is required
logical timgw
data timgw/.false./
! first imaginary frequency
real(8) lr_iw0
data lr_iw0/0.d0/
! last imaginary frequency
real(8) lr_iw1
data lr_iw1/80.d0/

real(8) fxca0
data fxca0/0.d0/
real(8) fxca1
data fxca1/0.d0/
integer nfxca
data nfxca/1/
integer fxctype
data fxctype/0/

! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! high-level switch: .true. if chi0 should be multiplied by 2
logical wannier_chi0_afm
data wannier_chi0_afm/.false./

! indices of response functions in global array f_response(:,:,:)
integer, parameter :: f_chi0                 = 1
integer, parameter :: f_chi                  = 2
integer, parameter :: f_chi_scalar           = 3
integer, parameter :: f_chi_pseudo_scalar    = 4
integer, parameter :: f_epsilon_matrix_GqGq  = 5
integer, parameter :: f_epsilon_scalar_GqGq  = 6
integer, parameter :: f_inv_epsilon_inv_GqGq = 7
integer, parameter :: f_epsilon_eff          = 8
integer, parameter :: f_epsilon_eff_scalar   = 9
integer, parameter :: f_sigma                = 10
integer, parameter :: f_sigma_scalar         = 11
integer, parameter :: f_loss                 = 12
integer, parameter :: f_loss_scalar          = 13
integer, parameter :: f_chi0_wann            = 14
integer, parameter :: f_chi_wann             = 15
integer, parameter :: f_epsilon_eff_wann     = 16
integer, parameter :: f_sigma_wann           = 17
integer, parameter :: f_loss_wann            = 18
integer, parameter :: f_epsilon_inv_GqGq     = 19

integer, parameter :: nf_response            = 19
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: u4(:,:,:,:)
logical screenu4
data screenu4/.true./

complex(8), allocatable :: gw_self_energy(:,:,:)

contains

subroutine genchi0blh(ikloc,ngq,w,chi0w)
use modmain
use mod_nrkp
use mod_expigqr
implicit none
! arguments
integer, intent(in) :: ikloc
integer, intent(in) :: ngq
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngq,ngq)
! local variables
logical l1
integer i,ist1,ist2,ik,jk,ig
real(8) t1,t2
complex(8), allocatable :: wt(:)
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblh(ikloc)))
wt(:)=zzero
do i=1,nmegqblh(ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
    if (abs(t1).gt.1d-6) then
      t2=sign(scissor,t1)
      wt(i)=t1/(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)-t2+w)
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngq
  megqblh2(1:nmegqblh(ikloc),ig)=dconjg(megqblh(1:nmegqblh(ikloc),ig,ikloc))*wt(1:nmegqblh(ikloc))
enddo
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngq,ngq,nmegqblh(ikloc),zone,megqblh(1,1,ikloc),nstsv*nstsv,&
  &megqblh2(1,1),nstsv*nstsv,zone,chi0w(1,1),ngq)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt)
return
end subroutine

end module
