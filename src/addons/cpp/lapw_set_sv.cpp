#include "lapw.h"

void lapw_set_sv(lapw_wave_functions& wf, double *evalfv_, tensor<complex16,2>& h)
{
    unsigned int ngk = wf.kp->ngk;

    memset(&h(0, 0), 0, h.size() * sizeof(complex16));

    if (p.ndmag > 0)
    {
        tensor<complex16,2> wfb(p.size_wfmt + ngk, p.nstfv);
        memset(&wfb(0, 0), 0, wfb.size() * sizeof(complex16));
        
        for (unsigned int ias = 0; ias < p.natmtot; ias++)
        {
            int sz = geometry.atoms[ias].species->ci.size();
            tensor<complex16,2> zm(sz, sz);
            memset(&zm(0, 0), 0, zm.size() * sizeof(complex16));
            for (int j2 = 0; j2 < sz; j2++)
            {
                int lm2 = geometry.atoms[ias].species->ci[j2].lm;
                int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
                for (int j1 = 0; j1 <= j2; j1++)
                {
                    int lm1 = geometry.atoms[ias].species->ci[j1].lm;
                    int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;
                    p.L3_sum_gntyry(lm1, lm2, &p.beffrad(0, idxrf1, idxrf2, ias, 0), zm(j1, j2));
                    // conjugate to call zhemm with upper triangular part
                    zm(j1, j2) = conj(zm(j1, j2));
                }
            }
            /* Warning!!! 
               This is the canonical code (zm should not be conjugated and computed for both upper an lower parts).
               
               zhemm call should be checked very carefully

            for (int i = 0; i < p.nstfv; i++)
                for (int j2 = 0; j2 < sz; j2++)
                    for (int j1 = 0; j1 < sz; j1++)
                        wfb(geometry.atoms[ias].offset_wfmt + j2, i) += zm(j1,j2) * wf.scalar_wf(geometry.atoms[ias].offset_wfmt + j1, i);
            */
 
            zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0), sz, &wf.scalar_wf(geometry.atoms[ias].offset_wfmt, 0), 
                wf.scalar_wf.size(0), zzero, &wfb(geometry.atoms[ias].offset_wfmt, 0), wfb.size(0));
            
        }
        
        std::vector<complex16> zfft(p.ngrtot);
        for (unsigned int i = 0; i < p.nstfv; i++)
        {
            memset(&zfft[0], 0, p.ngrtot * sizeof(complex16));
            for (unsigned int ig = 0; ig < ngk; ig++) 
                zfft[wf.kp->idxgfft[ig]] = wf.scalar_wf(p.size_wfmt + ig, i);
                                        
            lapw_fft(1, &zfft[0]);
                       
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
              zfft[ir] *= (p.beffir(ir) * p.cfunir[ir]);
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ngk; ig++) 
                wfb(p.size_wfmt + ig, i) = zfft[wf.kp->idxgfft[ig]];
        }

        // compute <wf_i | (b * wf_j)>
        zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
            wf.scalar_wf.size(0), &wfb(0, 0), wfb.size(0), zzero, &h(0, 0), h.size(0));

        // copy to dn-dn block and change sign
        for (unsigned int i = 0; i < p.nstfv; i++)
            for (unsigned int j = 0; j < p.nstfv; j++) 
                h(j + p.nstfv, i + p.nstfv) = -h(j, i);

        unsigned int i = 0;
        for (unsigned int ispn = 0; ispn < 2; ispn++)
        {
            for (unsigned int ist = 0; ist < p.nstfv; ist++)
            {
                h(i, i) += evalfv_[ist];
                i++;
            }
        }
    }
}


