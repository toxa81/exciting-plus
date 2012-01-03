#include "lapw.h"

Geometry geometry;
Parameters p;

/*inline void L3_sum_gntyry(int lm1, int lm2, double *v, std::complex<double>& zsum)
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
                         tensor<std::complex<double>,2>& apwalmc,
                         tensor<double,5>& apwfr,
                         tensor<double,3>& apwdfr,
                         tensor<double,6>& haa,
                         tensor<double,5>& hloa,
                         tensor<double,4>& hlolo,
                         tensor<double,4>& hmltrad,
                         tensor<double,3>& oalo,
                         tensor<double,3>& ololo,
                         tensor<double,3>& ovlprad,
                         tensor<std::complex<double>,2>& h,
                         tensor<std::complex<double>,2>& o)
{
    std::vector< std::complex<double> > zv(ngp);
    tensor<std::complex<double>,2> zm(ngp, wfmt_size_apw);

    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j2 = 0; j2 < atom->ci_apw.size(); j2++)
        {
            int l2 = atom->ci_apw[j2].l;
            int io2 = atom->ci_apw[j2].order;
            int lm2 = atom->ci_apw[j2].lm;
            int idxrf2 = atom->ci_apw[j2].idxrf;

            memset(&zv[0], 0, ngp * sizeof(std::complex<double>));
            
            for (unsigned int j1 = 0; j1 < atom->ci_apw.size(); j1++)
            {
                int lm1 = atom->ci_apw[j1].lm;
                int idxrf1 = atom->ci_apw[j1].idxrf;
                
                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm2, &hmltrad(0, idxrf1, idxrf2, ias), zsum);

                if (abs(zsum) > 1e-14) 
                    for (int ig = 0; ig < ngp; ig++) 
                        zv[ig] += zsum * apwalmc(ig, atom->offset_apw + j1); 
            }
            
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
            {
                double t1 = 0.5 * pow(species->rmt, 2) * apwfr(species->nrmt - 1, 0, io1, l2, ias) * apwdfr(io2, l2, ias); 
                for (int ig = 0; ig < ngp; ig++) 
                    zv[ig] += t1 * apwalmc(ig, atom->offset_apw + atom->ci_apw_by_lmo(lm2, io1)); //conj(apwalm(ig, io1, lm2, ias));
            }
            memcpy(&zm(0, atom->offset_apw + j2), &zv[0], ngp * sizeof(std::complex<double>));
        }
    }
    std::complex<double> zone(1, 0);
    std::complex<double> zzero(0, 0);
    zgemm<gemm_worker>(0, 2, ngp, ngp, wfmt_size_apw, zone, &zm(0, 0), ngp, &apwalmc(0, 0), ngp, zzero, &h(0, 0), ldh);
    zgemm<gemm_worker>(0 ,2, ngp, ngp, wfmt_size_apw, zone, &apwalmc(0, 0), ngp, &apwalmc(0, 0), ngp, zzero, &o(0, 0), ldh);

    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j = 0; j < atom->ci_lo.size(); j++) // loop over columns (local-orbital block) 
        {
            int l = atom->ci_lo[j].l;
            int lm = atom->ci_lo[j].lm;
            int ilo = atom->ci_lo[j].idxlo;
            int idxrf = atom->ci_lo[j].idxrf;
            
            // apw-lo block of the Hamiltonian
            for (unsigned int j1 = 0; j1 < atom->ci_apw.size(); j1++) // loop over rows
            {
                int lm1 = atom->ci_apw[j1].lm;
                int idxrf1 = atom->ci_apw[j1].idxrf;
                
                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm, &hmltrad(0, idxrf1, idxrf, ias), zsum);
                        
                if (abs(zsum) > 1e-14)
                    for (int ig = 0; ig < ngp; ig++)
                        h(ig, ngp + atom->offset_lo + j) += zsum * apwalmc(ig, atom->offset_apw + j1);
            }
            // apw-lo block of the overlap matrix
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l].radial_solution_descriptors.size(); io1++)
                for (int ig = 0; ig < ngp; ig++)
                    o(ig, ngp + atom->offset_lo + j) += ovlprad(io1, ilo, ias) * apwalmc(ig, atom->offset_apw + atom->ci_apw_by_lmo(lm, io1)); //conj(apwalm(ig, io1, lm, ias)) * oalo(io1, ilo, ias);

            // lo-lo block of the Hamiltonian 
            for (unsigned int j1 = 0; j1 <= j; j1++)
            {
                int lm1 = atom->ci_lo[j1].lm;
                int ilo1 = atom->ci_lo[j1].idxlo;
                int idxrf1 = atom->ci_lo[j1].idxrf;
                
                std::complex<double> zsum(0, 0);
                L3_sum_gntyry(lm1, lm, &hmltrad(0, idxrf1, idxrf, ias), zsum);
    
                h(ngp + atom->offset_lo + j1, ngp + atom->offset_lo + j) += zsum;
                
                // lo-lo block of the overlap matrix
                if (lm1 == lm) 
                    o(ngp + atom->offset_lo + j1, ngp + atom->offset_lo + j) += ovlprad(apwordmax + ilo1, ilo, ias);
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
*/

/*extern "C" void FORTFUNC(lapw_seceqnfv)(int32_t *ngp_,
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
                                        double *hmltrad_,
                                        double *oalo_,
                                        double *ololo_,
                                        double *ovlprad_,
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
    tensor<double,4> hmltrad(hmltrad_, lmmaxvr, nrfmtmax, nrfmtmax, natmtot);
    tensor<double,3> oalo(oalo_, apwordmax, nlomax, natmtot);
    tensor<double,3> ololo(ololo_, nlomax, nlomax, natmtot);
    tensor<double,3> ovlprad(ovlprad_, apwordmax + nlomax, nlomax, natmtot);

    tensor<std::complex<double>,2> apwalm_compact(ngp,wfmt_size_apw);
    reorder_apwalm(ngp, apwalm_, apwalm_compact);

    lapw_seceqnfv_setup(ngp, ldh, igpig, vgpc, veffig, cfunig, apwalm, apwalm_compact, apwfr,
        apwdfr, haa, hloa, hlolo, hmltrad, oalo, ololo, ovlprad, h, o);

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


*/






