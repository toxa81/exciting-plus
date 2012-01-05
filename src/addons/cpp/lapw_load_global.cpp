#include "lapw.h"

extern "C" void FORTFUNC(lapw_load_global)(int *natmtot_,
                                           int *nspecies_,
                                           int *lmaxvr_,
                                           int *lmaxapw_,
                                           int *apwordmax_,
                                           int *nrmtmax_,
                                           int *ngkmax_,
                                           int *ngvec_,
                                           int *ngrtot_,
                                           int *nlomax_,
                                           int *ias2is_,
                                           int *intgv_,
                                           int *ivg_,
                                           int *ivgig_,
                                           std::complex<double> *cfunig_,
                                           std::complex<double> *gntyry_,
                                           int *nstfv_,
                                           int *nstsv_,
                                           int *nmatmax_,
                                           int *nrfmtmax_,
                                           int *ordrfmtmax_,
                                           double *evaltol_)
{
    p.natmtot = *natmtot_;
    p.nspecies = *nspecies_;
    p.lmaxvr = *lmaxvr_;
    p.lmmaxvr = pow(p.lmaxvr + 1, 2);
    p.lmaxapw = *lmaxapw_;
    p.lmmaxapw = pow(p.lmaxapw + 1, 2);
    p.apwordmax = *apwordmax_;
    p.nrmtmax = *nrmtmax_;
    p.ngkmax = *ngkmax_;
    p.ngvec = *ngvec_;
    p.ngrtot = *ngrtot_;
    p.nlomax = *nlomax_;
    p.nstfv = *nstfv_;
    p.nstsv = *nstsv_;
    p.nmatmax = *nmatmax_;
    p.nrfmtmax = *nrfmtmax_;
    p.ordrfmtmax = *ordrfmtmax_;
    p.evaltol = *evaltol_;
    p.intgv = tensor<int,2>(intgv_, 3, 2);
    p.ivgig = tensor<int,3>(ivgig_, t_index(p.intgv(0, 0), p.intgv(0, 1)), 
                                    t_index(p.intgv(1, 0), p.intgv(1, 1)),
                                    t_index(p.intgv(2, 0), p.intgv(2, 1)));
    p.ivg = tensor<int,2>(ivg_, 3, p.ngrtot);
    p.cfunig.resize(p.ngrtot);
    for (int ig = 0; ig < p.ngrtot; ig++)
        p.cfunig[ig] = cfunig_[ig];
        
    p.gntyry = tensor<std::complex<double>,3>(gntyry_, p.lmmaxvr, p.lmmaxapw, p.lmmaxapw);
    p.L3_gntyry = tensor<std::vector<int>,2>(p.lmmaxapw, p.lmmaxapw);
    
    for (int lm1 = 0; lm1 < p.lmmaxapw; lm1++)
        for (int lm2 = 0; lm2 < p.lmmaxapw; lm2++)
            for (int lm3 = 0; lm3 < p.lmmaxvr; lm3++) 
                if (abs(p.gntyry(lm3, lm1, lm2)) > 1e-14)
                    p.L3_gntyry(lm1, lm2).push_back(lm3);
    
    for (int i = p.intgv(0, 0); i <= p.intgv(0, 1); i++)
        for (int j = p.intgv(1, 0); j <= p.intgv(1, 1); j++)
            for (int k = p.intgv(2, 0); k <= p.intgv(2, 1); k++)
                p.ivgig(i, j, k) -= 1;
    
    for (int i = 0; i < p.nspecies; i++)
        geometry.species.push_back(Species());
    
    for (int i = 0; i < p.natmtot; i++)
        geometry.atoms.push_back(Atom(&geometry.species[ias2is_[i] - 1]));
};


