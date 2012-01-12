#include "lapw.h"

void lapw_set_o(kpoint& kp, tensor<complex16,2>& capwalm, tensor<complex16,2>& o)
{
    timer t("lapw_set_o");
    
    memset(&o(0, 0), 0, o.size() * sizeof(complex16));

    zgemm<gemm_worker>(0 ,2, kp.ngk, kp.ngk, p.size_wfmt_apw, zone, &capwalm(0, 0), capwalm.size(0), 
        &capwalm(0, 0), capwalm.size(0), zzero, &o(0, 0), o.size(0));

//#pragma omp parallel for default(shared)
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j2 = 0; j2 < species->size_ci_lo; j2++) // loop over columns (local-orbital block) 
        {
            int l2 = species->ci_lo[j2].l;
            int lm2 = species->ci_lo[j2].lm;
            int order2 = species->ci_lo[j2].order;
            
            // apw-lo block 
            for (unsigned int io1 = 0; io1 < species->apw_descriptors[l2].radial_solution_descriptors.size(); io1++)
                for (unsigned int ig = 0; ig < kp.ngk; ig++)
                    o(ig, kp.ngk + atom->offset_lo + j2) += p.ovlprad(l2, io1, order2, ias) * capwalm(ig, atom->offset_apw + species->ci_by_lmo(lm2, io1)); 

            // lo-lo block
            for (unsigned int j1 = 0; j1 <= j2; j1++)
            {
                int lm1 = species->ci_lo[j1].lm;
                int order1 = species->ci_lo[j1].order;
                if (lm1 == lm2) 
                    o(kp.ngk + atom->offset_lo + j1, kp.ngk + atom->offset_lo + j2) += p.ovlprad(l2, order1, order2, ias);
            }
        }
    }
    
    //int iv[3];
//#pragma omp parallel for default(shared)
    for (unsigned int j2 = 0; j2 < kp.ngk; j2++) // loop over columns
    {
        //int ig2 = kp.idxg[j2];
        for (unsigned int j1 = 0; j1 <= j2; j1++) // for each column loop over rows
        {
            //int ig = idxG12(&kp, int j1, int j2);
            //for (int k = 0; k < 3; k++) iv[k] = p.ivg(k, kp.idxg[j1]) - p.ivg(k, ig2);
            //int ig = p.ivgig(iv[0], iv[1], iv[2]);
            //o(j1, j2) += p.cfunig[ig];
            o(j1, j2) += p.cfunig[idxG12(kp, j1, j2)];
        }
    }
}

