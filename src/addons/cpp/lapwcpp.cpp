#include <iostream>
#include <algorithm>
#include <complex>
#include <stdarg.h>
#include "tensor.h"
#include "lapw.h"

#define FORTFUNC(x) x##_

int natmtot;
int nspecies;
int lmaxvr;
int lmmaxvr;
int lmaxapw;
int lmmaxapw;
int apwordmax;
int ngkmax;
int ngvec;
int ngrtot;
int nlomax;

Geometry geometry;

tensor<int,1> ias2is;
tensor<int,1> ias2ia;
tensor<int,2> intgv;
tensor<int,2> ivg;
tensor<int,3> ivgig;
tensor<int,2> idxlm(t_index(0, 50), t_index(-50, 50));
tensor<std::complex<double>,3> gntyry;

extern "C" void FORTFUNC(lapw_load_global)(
    int *natmtot_,
    int *nspecies_,
    int *lmaxvr_,
    int *lmaxapw_,
    int *apwordmax_,
    int *ngkmax_,
    int *ngvec_,
    int *ngrtot_,
    int *nlomax_,
    int *ias2is_,
    int *ias2ia_,
    int *intgv_,
    int *ivg_,
    int *ivgig_,
    std::complex<double> *gntyry_)
{
    natmtot = *natmtot_;
    nspecies = *nspecies_;
    lmaxvr = *lmaxvr_;
    lmmaxvr = pow(lmaxvr + 1, 2);
    lmaxapw = *lmaxapw_;
    lmmaxapw = pow(lmaxapw + 1, 2);
    apwordmax = *apwordmax_;
    ngkmax = *ngkmax_;
    ngvec = *ngvec_;
    ngrtot = *ngrtot_;
    nlomax = *nlomax_;
    ias2is = tensor<int,1>(ias2is_, natmtot);
    ias2ia = tensor<int,1>(ias2ia_, natmtot);
    intgv = tensor<int,2>(intgv_, 3, 2);
    ivgig = tensor<int,3>(ivgig_, t_index(intgv(0, 0), intgv(0, 1)), 
                                  t_index(intgv(1, 0), intgv(1, 1)),
                                  t_index(intgv(2, 0), intgv(2, 1)));
    ivg = tensor<int,2>(ivg_, 3, ngrtot);
    gntyry = tensor<std::complex<double>,3>(gntyry_, lmmaxvr, lmmaxapw, lmmaxapw);
    
    for (int i = 0; i < natmtot; i++) 
    {
        ias2is(i) -= 1;
        ias2ia(i) -= 1;
    }
    
    for (int i = intgv(0, 0); i <= intgv(0, 1); i++)
        for (int j = intgv(1, 0); j <= intgv(1, 1); j++)
            for (int k = intgv(2, 0); k <= intgv(2, 1); k++)
                ivgig(i, j, k) -= 1;
    
    for (int i = 0; i < nspecies; i++)
        geometry.species.push_back(Species());
    
    for (int i = 0; i < natmtot; i++)
        geometry.atoms.push_back(Atom(&geometry.species[ias2is(i)]));
};

extern "C" void FORTFUNC(lapw_load_species)(
    int *is,
    int *nlorb,
    int *lorbl)
{
    for (int i = 0; i < *nlorb; i++) 
        geometry.species[*is-1].lo_descriptors.push_back(radial_l_channel_descriptor(lorbl[i]));
}

extern "C" void FORTFUNC(lapw_init)()
{
    for (int l = 0, lm = 0; l <= 50; l++) 
        for (int m = -l; m <= l; m++) idxlm(l, m) = lm++;
        
    for (unsigned int ias = 0, offset = 0; ias < geometry.atoms.size(); ias++)
    {
        geometry.atoms[ias].local_orbital_offset = offset;
        for (unsigned int ilo = 0; ilo < geometry.atoms[ias].species->lo_descriptors.size(); ilo++)
        {
            int l = geometry.atoms[ias].species->lo_descriptors[ilo].l;
            for (int m = -l; m <= l; m++, offset++) 
                geometry.atoms[ias].idxlo.push_back(local_orbital_index(l, m, idxlm(l, m), ilo));
        }
    }
}

extern "C" void FORTFUNC(setup_fv_hmlt_v1)(
    int *ngp,
    int *ldh,
    int *ncolh,
    int *igpig_,
    double *vgpc_,
    std::complex<double> *veffig_,
    std::complex<double> *cfunig_,
    double *hlolo_,
    std::complex<double> *h_) 
{
    tensor<double,2> vgpc(vgpc_, 3, ngkmax);
    tensor<int,1> igpig(igpig_, ngkmax);
    tensor<std::complex<double>,1> veffig(veffig_, ngvec);
    tensor<std::complex<double>,1> cfunig(cfunig_, ngrtot);
    tensor<std::complex<double>,2> h(h_, *ldh, *ncolh);
    tensor<double,4> hlolo(hlolo_, lmmaxvr, nlomax, nlomax, natmtot);

    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        for (unsigned int j = 0; j < geometry.atoms[ias].idxlo.size(); j++) 
        {
            int l = geometry.atoms[ias].idxlo[j].l;
            int m = geometry.atoms[ias].idxlo[j].m;
            int lm = geometry.atoms[ias].idxlo[j].lm;
            int ilo = geometry.atoms[ias].idxlo[j].ilo;
            for (unsigned int i = 0; i <= j; i++)
            {
                int l1 = geometry.atoms[ias].idxlo[i].l;
                int m1 = geometry.atoms[ias].idxlo[i].m;
                int lm1 = geometry.atoms[ias].idxlo[i].lm;
                int ilo1 = geometry.atoms[ias].idxlo[i].ilo;
                std::complex<double> zsum(0,0);
                
                int m3 = abs(m1 - m); // only +/- m3 coefficients are non-zero
                int l3 = (l + l1 + m3) % 2 + m3; // starting condition for l3: mod(l+l1+l3,2)==0 and m3<=l3 
                while (l3 <= std::min(lmaxvr, l + l1)) 
                {
                    if (m3 == 0) 
                    {
                        int lm3 = l3 * l3 + l3;
                        zsum += gntyry(lm3, lm1, lm) * hlolo(lm3, ilo1, ilo, ias);
                    } 
                    else 
                    {
                        int lm3 = l3 * l3 + l3 - m3;
                        zsum += gntyry(lm3, lm1, lm) * hlolo(lm3, ilo1, ilo, ias);
                        lm3 += 2*m3;
                        zsum += gntyry(lm3, lm1, lm) * hlolo(lm3, ilo1, ilo, ias);
                    }
                    l3 += 2;
                }
                h(i + *ngp + geometry.atoms[ias].local_orbital_offset,
                  j + *ngp + geometry.atoms[ias].local_orbital_offset) += zsum;


            }
        }
    }

 
            
    
    /*
    int ig1, ig2, iv[3];
    for (int j = 0; j < *ngp; j++) // loop over columns
    {
        ig1 = igpig(j) - 1;
        for (int i = 0; i <= j; i++) // for each column loop over rows
        {
            ig2 = igpig(i) - 1;
            for (int k = 0; k < 3; k++) iv[k] = ivg(k, ig2) - ivg(k, ig1);
            int ig = ivgig(iv[0], iv[1], iv[2]);
            double t1 = 0.5 * (vgpc(0, i) * vgpc(0, j) + 
                               vgpc(1, i) * vgpc(1, j) + 
                               vgpc(2, i) * vgpc(2, j));
            h(i, j) += veffig(ig) + t1 * cfunig(ig);
        }
    }*/
  
}

