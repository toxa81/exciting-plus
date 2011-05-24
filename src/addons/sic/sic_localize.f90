subroutine sic_localize(fout,wanprop)
use modmain
use mod_sic
implicit none
integer, intent(in) :: fout
real(8), intent(out) :: wanprop(nwanprop,sic_wantran%nwan)
!
integer iter,i,j,n,n1,j1,n2,j2,ispn,ir,ikloc
real(8) tot_diff
complex(8), allocatable :: um(:,:),um1(:,:),um0(:,:)
complex(8), allocatable :: wlm_new(:,:)
!
allocate(um(sic_wantran%nwan,sic_wantran%nwan))
allocate(um0(sic_wantran%nwan,sic_wantran%nwan))
allocate(um1(sic_wantran%nwan,sic_wantran%nwan))
allocate(wlm_new(lmmaxwan,sic_wantran%nwan))
if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("on-site localization")')
  write(fout,'(80("="))')
endif

um0=zzero
do i=1,sic_wantran%nwan
  um0(i,i)=zone
enddo

do iter=1,sic_niter_u0
  if (iter.gt.1) then
    call sic_genvme_dotp(.true.)
    tot_diff=0.d0
    call sic_exp_grad_u((/0.d0,0.d0,0.d0/),sic_u0_eps,tot_diff,um)
    if (wproc) then
      write(fout,'("  total localization error : ",G18.10)')tot_diff
      !write(fout,'("    SIC potential energy : ",G18.10)')&
      !  sum(0.5d0*wanprop(wp_vha,:))+sum(wanprop(wp_exc,:))
    endif
    um1=zzero
    do j1=1,sic_wantran%nwan
      do j2=1,sic_wantran%nwan
        do j=1,sic_wantran%nwan
          um1(j1,j2)=um1(j1,j2)+um0(j1,j)*um(j,j2)
        enddo
      enddo
    enddo
    um0(:,:)=um1(:,:)
    do ispn=1,nspinor
      do ir=1,s_nr
        wlm_new=zzero
        do j=1,sic_wantran%nwan 
          do j1=1,sic_wantran%nwan
            wlm_new(:,j)=wlm_new(:,j)+um(j1,j)*s_wlm(:,ir,ispn,j1)
          enddo
        enddo
        s_wlm(:,ir,ispn,:)=wlm_new(:,:)
      enddo
    enddo    
  endif
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)  
    if (sic_apply(n).eq.2) then
      call timer_start(t_sic_wan_pot)
      call sic_genpot(s_wlm(1,1,1,j),s_wvlm(1,1,1,j),wanprop(1,j))
      call timer_stop(t_sic_wan_pot)
    endif
  enddo
enddo !iter
deallocate(um,um1,wlm_new)
allocate(um(nwantot,nwantot),um1(nwantot,nwantot))
um=zzero
do i=1,nwantot
  um(i,i)=zone
enddo
do j1=1,sic_wantran%nwan
  n1=sic_wantran%iwan(j1)
  do j2=1,sic_wantran%nwan
    n2=sic_wantran%iwan(j2)
    um(n1,n2)=um0(j1,j2)
  enddo
enddo
! update k-independent part of sic_wan_umtrx
do ikloc=1,nkptnrloc
  um1=zzero
  do n1=1,nwantot
    do n2=1,nwantot
      do n=1,nwantot
        um1(n1,n2)=um1(n1,n2)+sic_wan_umtrx(n1,n,ikloc)*um(n,n2)
      enddo
    enddo
  enddo
  sic_wan_umtrx(:,:,ikloc)=um1(:,:)
enddo
deallocate(um,um0,um1)
return
end subroutine

