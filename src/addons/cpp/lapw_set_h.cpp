#include "lapw.h"

void lapw_set_h(kpoint& kp, tensor<complex16,2>& capwalm, tensor<complex16,2>& h)
{
    timer t("lapw_set_h");
    
    memset(&h(0, 0), 0, h.size() * sizeof(complex16));

    tensor<complex16,2> zm(kp.ngk, p.size_wfmt_apw);

#pragma omp parallel default(shared)
{
    std::vector<complex16> zv(kp.ngk);
#pragma omp for
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        // precompute apw block
        for (unsigned int j2 = 0; j2 < species->size_ci_apw; j2++)
        {
            memset(&zv[0], 0, kp.ngk * sizeof(complex16));
            
            int lm2 = species->ci[j2].lm;
            int idxrf2 = species->ci[j2].idxrf;
            
            for (unsigned int j1 = 0; j1 < species->size_ci_apw; j1++)
            {
                int lm1 = species->ci[j1].lm;
                int idxrf1 = species->ci[j1].idxrf;
                
                complex16 zsum(0, 0);
                p.L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf1, idxrf2, ias), zsum);

                if (abs(zsum) > 1e-14) 
                    for (unsigned int ig = 0; ig < kp.ngk; ig++) 
                        zv[ig] += zsum * capwalm(ig, atom->offset_apw + j1); 
            }
            
            int l2 = species->ci[j2].l;
            int io2 = species->ci[j2].order;
            
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
            {
                double t1 = 0.5 * pow(species->rmt, 2) * p.apwfr(species->nrmt - 1, 0, io1, l2, ias) * p.apwdfr(io2, l2, ias); 
                for (unsigned int ig = 0; ig < kp.ngk; ig++) 
                    zv[ig] += t1 * capwalm(ig, atom->offset_apw + species->ci_by_lmo(lm2, io1));
            }
            memcpy(&zm(0, atom->offset_apw + j2), &zv[0], kp.ngk * sizeof(complex16));
        }

        for (unsigned int j2 = 0; j2 < species->size_ci_lo; j2++) // loop over columns (local-orbital block) 
        {
            int lm2 = species->ci_lo[j2].lm;
            int idxrf2 = species->ci_lo[j2].idxrf;
            
            // apw-lo block
            for (unsigned int j1 = 0; j1 < species->size_ci_apw; j1++) // loop over rows
            {
                int lm1 = species->ci[j1].lm;
                int idxrf1 = species->ci[j1].idxrf;
                
                complex16 zsum(0, 0);
                p.L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf2, idxrf1, ias), zsum);
                        
                if (abs(zsum) > 1e-14)
                    for (unsigned int ig = 0; ig < kp.ngk; ig++)
                        h(ig, kp.ngk + atom->offset_lo + j2) += zsum * capwalm(ig, atom->offset_apw + j1);
            }

            // lo-lo block 
            for (unsigned int j1 = 0; j1 <= j2; j1++)
            {
                int lm1 = species->ci_lo[j1].lm;
                int idxrf1 = species->ci_lo[j1].idxrf;
                
                complex16 zsum(0, 0);
                p.L3_sum_gntyry(lm1, lm2, &p.hmltrad(0, idxrf1, idxrf2, ias), zsum);
    
                h(kp.ngk + atom->offset_lo + j1, kp.ngk + atom->offset_lo + j2) += zsum;
            }
        }
    } // ias
}
    zgemm<gemm_worker>(0, 2, kp.ngk, kp.ngk, p.size_wfmt_apw, zone, &zm(0, 0), zm.size(0), 
        &capwalm(0, 0), capwalm.size(0), zzero, &h(0, 0), h.size(0));

    //int iv[3];
    for (unsigned int j2 = 0; j2 < kp.ngk; j2++) // loop over columns
    {
        //int ig2 = kp.idxg[j2];
        double v2[3];
        for (int k = 0; k < 3; k++) v2[k] = kp.vgkc(k, j2);
        for (unsigned int j1 = 0; j1 <= j2; j1++) // for each column loop over rows
        {
            //for (int k = 0; k < 3; k++) iv[k] = p.ivg(k, kp.idxg[j1]) - p.ivg(k, ig2);
            //int ig = p.ivgig(iv[0], iv[1], iv[2]);
            int ig = idxG12(kp, j1, j2);
            double t1 = 0.5 * (kp.vgkc(0, j1) * v2[0] + kp.vgkc(1, j1) * v2[1] + kp.vgkc(2, j1) * v2[2]);
            h(j1, j2) += p.veffig(ig) + t1 * p.cfunig[ig];
        }
    }
}


