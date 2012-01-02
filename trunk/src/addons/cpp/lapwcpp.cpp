#include <iostream>
#include <algorithm>
#include <complex>
#include <stdarg.h>
#include <stdint.h>
#include "tensor.h"
#include "lapw.h"
#include "linalg.h"
#include "config.h"

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
int nrmtmax;
int nstfv;
int nstsv;
int nmatmax;
int nrfmtmax;
double evaltol;

Geometry geometry;

tensor<int,1> ias2is;
tensor<int,1> ias2ia;
tensor<int,2> intgv;
tensor<int,2> ivg;
tensor<int,3> ivgig;
tensor<std::complex<double>,3> gntyry;
std::vector< std::vector<int> > gntyry_compact;

int wfmt_size_apw;
int wfmt_size_lo;
int wfmt_size;

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
                                           int *ias2ia_,
                                           int *intgv_,
                                           int *ivg_,
                                           int *ivgig_,
                                           std::complex<double> *gntyry_,
                                           int *nstfv_,
                                           int *nstsv_,
                                           int *nmatmax_,
                                           int *nrfmtmax_,
                                           double *evaltol_)
{
    natmtot = *natmtot_;
    nspecies = *nspecies_;
    lmaxvr = *lmaxvr_;
    lmmaxvr = pow(lmaxvr + 1, 2);
    lmaxapw = *lmaxapw_;
    lmmaxapw = pow(lmaxapw + 1, 2);
    apwordmax = *apwordmax_;
    nrmtmax = *nrmtmax_;
    ngkmax = *ngkmax_;
    ngvec = *ngvec_;
    ngrtot = *ngrtot_;
    nlomax = *nlomax_;
    nstfv = *nstfv_;
    nstsv = *nstsv_;
    nmatmax = *nmatmax_;
    nrfmtmax = *nrfmtmax_;
    evaltol = *evaltol_;
    ias2is = tensor<int,1>(ias2is_, natmtot);
    ias2ia = tensor<int,1>(ias2ia_, natmtot);
    intgv = tensor<int,2>(intgv_, 3, 2);
    ivgig = tensor<int,3>(ivgig_, t_index(intgv(0, 0), intgv(0, 1)), 
                                  t_index(intgv(1, 0), intgv(1, 1)),
                                  t_index(intgv(2, 0), intgv(2, 1)));
    ivg = tensor<int,2>(ivg_, 3, ngrtot);
    gntyry = tensor<std::complex<double>,3>(gntyry_, lmmaxvr, lmmaxapw, lmmaxapw);
    gntyry_compact.resize(lmmaxapw*lmmaxapw);
    for (int lm1 = 0; lm1 < lmmaxapw; lm1++)
        for (int lm2 = 0; lm2 < lmmaxapw; lm2++)
            for (int lm3 = 0; lm3 < lmmaxvr; lm3++) 
                if (abs(gntyry(lm3, lm1, lm2)) > 1e-14)
                    gntyry_compact[lm1 + lm2 * lmmaxapw].push_back(lm3);
    
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

extern "C" void FORTFUNC(lapw_load_species)(int *is_,
                                            int *nlorb,
                                            int *lorbl,
                                            int *apword,
                                            double* rmt,
                                            int *nrmt)
{
    int is = *is_ - 1;
    
    geometry.species[is].rmt = *rmt;
    geometry.species[is].nrmt = *nrmt;
    
    for (int i = 0; i < *nlorb; i++) 
        geometry.species[is].lo_descriptors.push_back(radial_l_channel_descriptor(lorbl[i]));
    
    for (int l = 0; l <= lmaxapw; l++) 
    {
        radial_l_channel_descriptor lch(l);
        for (int io = 0; io < apword[l]; io++)
        {
            radial_solution_descriptor rs;
            lch.radial_solution_descriptors.push_back(rs);
        }
        geometry.species[is].apw_descriptors.push_back(lch);
    }
}

extern "C" void FORTFUNC(lapw_init)()
{
    wfmt_size_apw = 0;
    wfmt_size_lo = 0;
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;
        
        atom->offset_apw = wfmt_size_apw;
        atom->offset_lo = wfmt_size_lo;

        std::vector<int> order(lmaxapw + 1, 0);
        int idxrf = 0;
        
        for (int l = 0; l <= lmaxapw; l++)
        {
            for (unsigned int io = 0; io < species->apw_descriptors[l].radial_solution_descriptors.size(); io++)
            {
                for (int m = -l; m <= l; m++)
                    atom->ci_apw.push_back(mtci(l, m, order[l], idxrf));
                
                order[l]++;
                idxrf++;
            }
        }
        wfmt_size_apw += atom->ci_apw.size();
        
        for (unsigned int ilo = 0; ilo < species->lo_descriptors.size(); ilo++)
        {
            int l = species->lo_descriptors[ilo].l;
            for (int m = -l; m <= l; m++) 
                atom->ci_lo.push_back(mtci(l, m, order[l], idxrf, ilo));

            order[l]++;
            idxrf++;
        }
        wfmt_size_lo += atom->ci_lo.size();
        
        //std::cout << "Atom : " << ias << " size_apw=" << atom->size_apw << " size_lo=" << atom->size_lo
        //    << " offset_apw=" << atom->offset_apw << " offset_lo=" << atom->offset_lo << std::endl;    
    }
    
    wfmt_size = wfmt_size_apw + wfmt_size_lo;
}

// \sum_{L3} <L1|L3|L2> * v[L3]
/*inline std::complex<double> L3_sum_gntyry(int l1, int m1, int l2, int m2, std::complex<double> *gnt, double *v)
{
    std::complex<double> zsum(0,0);
    int m3 = abs(m1 - m2); // only +/- m3 coefficients are non-zero
    int l3 = (l1 + l2 + m3) % 2 + m3; // starting condition for l3: mod(l+l1+l3,2)==0 and l3>=m3 
    if (m3)
    {
        while (l3 <= std::min(lmaxvr, l1 + l2)) 
        {
            int lm3 = l3 * l3 + l3;
            zsum += gnt[lm3 - m3] * v[lm3 - m3];
            zsum += gnt[lm3 + m3] * v[lm3 + m3];
            l3 += 2;
        }
    }
    else
    {
        while (l3 <= std::min(lmaxvr, l1 + l2)) 
        {
            int lm3 = l3 * l3 + l3;
            zsum += gnt[lm3] * v[lm3];
            l3 += 2;
        }
    }
    return zsum;
}*/

inline void L3_sum_gntyry(int lm1, int lm2, double *v, std::complex<double>& zsum)
{
    for (unsigned int k = 0; k < gntyry_compact[lm1 + lm2 * lmmaxapw].size(); k++)
    {
        int lm3 = gntyry_compact[lm1 + lm2 * lmmaxapw][k];
        zsum += gntyry(lm3, lm1, lm2) * v[lm3];
    }
}

void lapw_seceqnfv_setup(int ngp,
                         int ldh,
                         tensor<int,1>& igpig,
                         tensor<double,2>& vgpc,
                         tensor<std::complex<double>,1>& veffig,
                         tensor<std::complex<double>,1>& cfunig,
                         tensor<std::complex<double>,4>& apwalm,
                         tensor<double,5>& apwfr,
                         tensor<double,3>& apwdfr,
                         tensor<double,6>& haa,
                         tensor<double,5>& hloa,
                         tensor<double,4>& hlolo,
                         tensor<double,3>& oalo,
                         tensor<double,3>& ololo,
                         tensor<std::complex<double>,2>& h,
                         tensor<std::complex<double>,2>& o)
{
    std::vector< std::complex<double> > zv1(ngp);
    std::vector< std::complex<double> > zm1(ngp * wfmt_size_apw);
    std::vector< std::complex<double> > zm2(ngp * wfmt_size_apw);
    

    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j2 = 0; j2 < atom->ci_apw.size(); j2++)
        {
            int l2 = atom->ci_apw[j2].l;
            int io2 = atom->ci_apw[j2].order;
            int lm2 = atom->ci_apw[j2].lm;

            memset(&zv1[0], 0, ngp * sizeof(std::complex<double>));
            
            for (unsigned int j1 = 0; j1 < atom->ci_apw.size(); j1++)
            {
                int l1 = atom->ci_apw[j1].l;
                int io1 = atom->ci_apw[j1].order;
                int lm1 = atom->ci_apw[j1].lm;

                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm2, &haa(0, io1, l1, io2, l2, ias), zsum);

                if (abs(zsum) > 1e-14) 
                    for (int ig = 0; ig < ngp; ig++) 
                        zv1[ig] += zsum * conj(apwalm(ig, io1, lm1, ias));
            }
            
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
            {
                double t1 = 0.5 * pow(species->rmt, 2) * apwfr(species->nrmt - 1, 0, io1, l2, ias) * apwdfr(io2, l2, ias); 
                for (int ig = 0; ig < ngp; ig++) 
                    zv1[ig] += t1 * conj(apwalm(ig, io1, lm2, ias));
            }
            memcpy(&zm1[(atom->offset_apw + j2) * ngp], &zv1[0], ngp * sizeof(std::complex<double>));
            for (int ig = 0; ig < ngp; ig++)
                zm2[(atom->offset_apw + j2) * ngp + ig] = conj(apwalm(ig, io2, lm2, ias));
        }
    }
    std::complex<double> zone(1,0);
    std::complex<double> zzero(0,0);
    zgemm<gemm_worker>(0, 2, ngp, ngp, wfmt_size_apw, zone, &zm1[0], ngp, &zm2[0], ngp, zzero, &h(0,0), ldh);
    zgemm<gemm_worker>(0 ,2, ngp, ngp, wfmt_size_apw, zone, &zm2[0], ngp, &zm2[0], ngp, zzero, &o(0,0), ldh);

    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j = 0; j < atom->ci_lo.size(); j++) // loop over columns (local-orbital block) 
        {
            int l = atom->ci_lo[j].l;
            int lm = atom->ci_lo[j].lm;
            int ilo = atom->ci_lo[j].idxlo;
            
            // apw-lo block of the Hamiltonian
            for (unsigned int j1 = 0; j1 < atom->ci_apw.size(); j1++)
            {
                int l1 = atom->ci_apw[j1].l;
                int io1 = atom->ci_apw[j1].order;
                int lm1 = atom->ci_apw[j1].lm;
                
                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm, &hloa(0, ilo, io1, l1, ias), zsum);
                        
                if (abs(zsum) > 1e-14)
                    for (int ig = 0; ig < ngp; ig++)
                        h(ig, ngp + atom->offset_lo + j) += conj(apwalm(ig, io1, lm1, ias)) * zsum;
            }
            // apw-lo block of the overlap matrix
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l].radial_solution_descriptors.size(); io1++)
                for (int ig = 0; ig < ngp; ig++)
                    o(ig, ngp + atom->offset_lo + j) += conj(apwalm(ig, io1, lm, ias)) * oalo(io1, ilo, ias);

            // lo-lo block of the Hamiltonian 
            for (unsigned int i = 0; i <= j; i++)
            {
                int lm1 = atom->ci_lo[i].lm;
                int ilo1 = atom->ci_lo[i].idxlo;
                
                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm, &hlolo(0, ilo1, ilo, ias), zsum);
    
                h(ngp + atom->offset_lo + i, ngp + atom->offset_lo + j) += zsum;
                
                // lo-lo block of the overlap matrix
                if (lm1 == lm) 
                    o(ngp + atom->offset_lo + i, ngp + atom->offset_lo + j) += ololo(ilo1, ilo, ias);
            }
        }
    }
    
    int iv[3];
    for (int j = 0; j < ngp; j++) // loop over columns
    {
        int ig2 = igpig(j) - 1;
        double v2[3];
        for (int k = 0; k < 3; k++) v2[k] = vgpc(k, j);
        for (int i = 0; i <= j; i++) // for each column loop over rows
        {
            for (int k = 0; k < 3; k++) iv[k] = ivg(k, igpig(i) - 1) - ivg(k, ig2);
            int ig = ivgig(iv[0], iv[1], iv[2]);
            double t1 = 0.5 * (vgpc(0, i) * v2[0] + vgpc(1, i) * v2[1] + 
                               vgpc(2, i) * v2[2]);
            h(i, j) += veffig(ig) + t1 * cfunig(ig);
            o(i, j) += cfunig(ig);
        }
    }
}

extern "C" void FORTFUNC(lapw_seceqnfv)(
    int32_t *ngp_,
    int32_t *ldh_,
    int32_t *ncolh,
    int32_t *igpig_,
    double *vgpc_,
    std::complex<double> *veffig_,
    std::complex<double> *cfunig_,
    std::complex<double> *apwalm_,
    double *apwfr_,
    double *apwdfr_,
    double *haa_,
    double *hloa_,
    double *hlolo_,
    double *oalo_,
    double *ololo_,
    std::complex<double> *h_,
    std::complex<double> *o_,
    double *evalfv_,
    std::complex<double> *z_) 
{
    int ngp = *ngp_;
    int ldh = *ldh_;
    tensor<double,2> vgpc(vgpc_, 3, ngkmax);
    tensor<int,1> igpig(igpig_, ngkmax);
    tensor<std::complex<double>,1> veffig(veffig_, ngvec);
    tensor<std::complex<double>,1> cfunig(cfunig_, ngrtot);
    tensor<std::complex<double>,2> h(h_, ldh, *ncolh);
    tensor<std::complex<double>,2> o(o_, ldh, *ncolh);
    tensor<std::complex<double>,2> z(z_, nmatmax, nstfv);
    tensor<std::complex<double>,4> apwalm(apwalm_, ngkmax, apwordmax, lmmaxapw, natmtot);
    tensor<double,5> apwfr(apwfr_, nrmtmax, 2, apwordmax, lmaxapw + 1, natmtot);
    tensor<double,3> apwdfr(apwdfr_, apwordmax, lmaxapw + 1, natmtot);
    tensor<double,6> haa(haa_, lmmaxvr, apwordmax, lmaxapw + 1, apwordmax, lmaxapw + 1, natmtot);
    tensor<double,5> hloa(hloa_, lmmaxvr, nlomax, apwordmax, lmaxapw + 1, natmtot);
    tensor<double,4> hlolo(hlolo_, lmmaxvr, nlomax, nlomax, natmtot);
    tensor<double,3> oalo(oalo_, apwordmax, nlomax, natmtot);
    tensor<double,3> ololo(ololo_, nlomax, nlomax, natmtot);

    lapw_seceqnfv_setup(ngp, ldh, igpig, vgpc, veffig, cfunig, apwalm, apwfr,
        apwdfr, haa, hloa, hlolo, oalo, ololo, h, o);

    return;
       
    tensor<std::complex<double>,2> h1;
    tensor<std::complex<double>,2> o1;

    if (check_evecfv) 
    {
        h1 = h;
        o1 = o;
    }

    zhegv<lapack_worker>(ldh, nstfv, evaltol, &h(0, 0), &o(0, 0), evalfv_, 
        &z(0, 0), nmatmax);

    if (check_evecfv) 
    {
        std::complex<double> zone(1, 0);
        std::complex<double> zzero(0, 0);
    
        // use o and h as temporary arrays
        for (int i = 0; i < nstfv; i++)
            for (int j = 0; j < ldh; j++)
                o(j, i) = -evalfv_[i] * z(j, i);
    
        zhemm<blas_worker>(0, 0, ldh, nstfv, zone, &h1(0, 0), ldh, &z(0, 0), 
            nmatmax, zzero, &h(0, 0), ldh);
        zhemm<blas_worker>(0, 0, ldh, nstfv, zone, &o1(0, 0), ldh, &o(0, 0),
            ldh, zone, &h(0, 0), ldh); 
    
        double L2norm = 0.0;
        for (int i = 0; i < nstfv; i++)
            for (int j = 0; j < ldh; j++) L2norm += abs(h(j, i) * conj(h(j, i)));
        
        std::cout << "check evecfv :" << std::endl
                  << "  number of bands = " << nstfv << std::endl
                  << "  total L2 diff = " << sqrt(L2norm) << std::endl;
    }
}









