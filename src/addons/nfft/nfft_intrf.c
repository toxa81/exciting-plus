#ifdef _NFFT_

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "nfft3.h"

nfsft_plan plan;

void nfsft_init_(int* lmax_, int* ntp_, double* tp) 
{
  int lmax=*lmax_;
  int ntp=*ntp_;
  int j;
  double th,ph;
  
  nfsft_precompute(lmax, 1000.0, 0U, 0U);
  //nfsft_init(&plan,lmax,lmmax);
  //nfsft_init_advanced(&plan, lmax, ntp, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
  //  NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED);
  nfsft_init_guru(&plan, lmax, ntp, NFSFT_MALLOC_X | NFSFT_MALLOC_F |
    NFSFT_MALLOC_F_HAT | NFSFT_NORMALIZED, PRE_PHI_HUT | PRE_PSI | FFTW_INIT | FFT_OUT_OF_PLACE, 6);

  printf(" K : %i\n",plan.plan_nfft.K);

  for (j=0; j<ntp; j++) {
    th=tp[2*j];
    ph=tp[2*j+1];
    if (ph>=0.0 && ph<M_PI) {
      plan.x[2*j]=ph/2.0/M_PI;
    } else if (ph>=M_PI && ph<2.0*M_PI) {
      plan.x[2*j]=ph/2.0/M_PI-1.0;
    } else {
      printf("Error(nfsft_init): phi=%f is out of range [0,2*Pi) \n",ph);
      exit(0);
    }
    if (th>=0.0 && th<=M_PI) {
      plan.x[2*j+1]=th/2.0/M_PI;
    } else {
      printf("Error(nfsft_init): theta=%f is out of range [0,Pi] \n",th);
      exit(0);
    }
  }
  nfsft_precompute_x(&plan);
}

void nfsft_ft_(int* lmax, int* ntp, int* nr, double complex* flm, 
  double complex* ftp) 
{
  int ir,ilm,itp,j,l,m;
  double phase;
  
  ilm=0;
  itp=0;
  for (ir=0; ir<*nr; ir++) {
    for (l=0; l<=*lmax; l++) {
      for (m=-l; m<=l; m++) {
        phase=1.0;
        if (m>0 && m%2==1) phase=-1.0;
        plan.f_hat[NFSFT_INDEX(l,m,&plan)]=phase*flm[ilm++];
      }
    }
    nfsft_trafo(&plan);
    //ndsft_trafo(&plan);
    for (j=0; j<*ntp; j++) ftp[itp++]=plan.f[j];
  }
} 

/*void nfsft_test_forward_transform_(int* lmax, int* ntp, int* nr, 
  double complex* flm, double complex* ftp) 
{
  int k,n,j,l,m;
  double phase;
  //int lmmax=pow(*lmax+1,2);
  int ir,ilm,jlm;
  
  ilm=0;
  for (ir=0; ir<*nr; ir++) {
    jlm=ilm;
    for (l=0; l<=*lmax; l++) {
      for (m=-l; m<=l; m++) {
        phase=1.0;
        //if (n>0 && n%2==1) phase=-1.0;
        plan.f_hat[NFSFT_INDEX(l,m,&plan)] = phase*flm[ilm++];
        if (ir == 0) {
          printf("l=%i m=%i, f_hat=%f %f\n",l,m,creal(plan.f_hat[NFSFT_INDEX(l,m,&plan)]),
            cimag(plan.f_hat[NFSFT_INDEX(l,m,&plan)]));
        }
      }
    }
    //ndsft_trafo(&plan);
    //ndsft_adjoint(&plan);
    nfsft_trafo(&plan);
    nfsft_adjoint(&plan);
    for (l=0; l<=*lmax; l++) {
      for (m=-l; m<=l; m++) {
        flm[jlm++] = plan.f_hat[NFSFT_INDEX(l,m,&plan)]*4*M_PI/(*ntp);
      }
    }
  }
} */

/*void nfsft_bt_(int* lmax, int* ntp, int* nr, double complex* ftp, 
  double complex* flm) 
{
  int ir,ilm,itp,j,l,m;
  double phase;
  
  itp=0;
  ilm=0;
  for (ir=0; ir<*nr; ir++) {
    for (j=0; j<*ntp; j++) plan.f[j]=ftp[itp++];
    nfsft_adjoint(&plan);
    //ndsft_adjoint(&plan);
    for (l=0; l<=*lmax; l++) {
      for (m=-l; m<=l; m++) {
        phase=1.0;
        //if (m>0 && m%2==1) phase=-1.0;
        flm[ilm++]=phase*plan.f_hat[NFSFT_INDEX(l,m,&plan)]*4*M_PI/(*ntp);
      }
    }
  }
}*/

void nfsft_bt_(int* lmax, int* ntp, int* nr, double complex* ftp, 
  double complex* flm) 
{
  int ir,ilm,itp,j,l,m,it;
  double phase;
  solver_plan_complex iplan;
  unsigned int method=CGNR;

  printf("M_total : %i\n",plan.M_total);
  printf("N_total : %i\n",plan.N_total);

  solver_init_advanced_complex(&iplan,(mv_plan_complex*)(&plan),method);
  //solver_init_complex(&iplan,(mv_plan_complex*)(&plan));

  itp=0;
  ilm=0;
  for (ir=0; ir<*nr; ir++) {
    for (j=0; j<*ntp; j++) iplan.y[j]=ftp[itp++];
    for (j=0; j<plan.N_total; j++) iplan.f_hat_iter[j]=0;
    solver_before_loop_complex(&iplan);
    //printf("first residual r=%e\n",iplan.dot_r_iter);

    for (it=0; it<26; it++) {
      printf("alpha : %f\n",iplan.alpha_iter);
      printf("beta : %f\n",iplan.beta_iter);
      solver_loop_one_step_complex(&iplan);
      printf("residual r=%e\n",iplan.dot_r_iter);
    }

    //nfsft_adjoint(&plan);
    //ndsft_adjoint(&plan);

    for (l=0; l<=*lmax; l++) {
      for (m=-l; m<=l; m++) {
        phase=1.0;
        //if (m>0 && m%2==1) phase=-1.0;
        flm[ilm++]=phase*iplan.f_hat_iter[NFSFT_INDEX(l,m,&plan)]; //*4*M_PI/(*ntp);
      }
    }
  }
  solver_finalize_complex(&iplan);
}

/*void nfsft_test_backward_transform_(int *lmax, int* ntp, int* nr, double complex* f_hat, double complex* f) {
  int k,n,j;
  double phase;
  int ir,itp,ilm,jtp;
  printf("[nfsft_test_backward_transform]\n");
  printf("  lmax :  %i\n",*lmax);
  printf("  ntp  :  %i\n",*ntp);
  printf("  nr   :  %i\n",*nr);
  printf("  N_total :  %i\n",plan.N_total);
  printf("  M_total :  %i\n",plan.M_total);
  
  itp=0;
  ilm=0;
  for (ir=0; ir<*nr; ir++) {
    jtp=itp;
    for (j=0; j<*ntp; j++) plan.f[j]=f[itp++];
    nfsft_adjoint(&plan);
    for (k=0; k<=*lmax; k++) {
      for (n=-k; n<=k; n++) {
        phase=1.0;
        if (n>0 && n%2==1) phase=-1.0;
        plan.f_hat[NFSFT_INDEX(k,n,&plan)]=phase*plan.f_hat[NFSFT_INDEX(k,n,&plan)]*4*M_PI/(*ntp);
      }
    }
    nfsft_trafo(&plan);
    for (j=0; j<*ntp; j++) f[jtp++]=plan.f[j];
  }
}*/

#endif
