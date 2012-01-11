#include "lapw.h"

void lapw_set_sv(lapw_wave_functions& wf, double *evalfv_, tensor<complex16,2>& h)
{
    unsigned int ngk = wf.kp->ngk;

    memset(&h(0, 0), 0, h.size() * sizeof(complex16));

    int nhwf;
    if (p.ndmag == 0) nhwf = 1;
    if (p.ndmag == 1) nhwf = 2;
    if (p.ndmag == 3) nhwf = 3;

    // product of the second-variational hamiltonian and a wave-function
    tensor<complex16,3> hwf(p.size_wfmt + ngk, p.nstfv, nhwf);
    memset(&hwf(0, 0, 0), 0, hwf.size() * sizeof(complex16));

    // compute product of magnetic field and wave-function 
    if (p.ndmag > 0)
    {
        for (unsigned int ias = 0; ias < p.natmtot; ias++)
        {
            int sz = geometry.atoms[ias].species->ci.size();
            
            tensor<complex16,3> zm(sz, sz, p.ndmag);
            memset(&zm(0, 0, 0), 0, zm.size() * sizeof(complex16));

            for (int j2 = 0; j2 < sz; j2++)
            {
                int lm2 = geometry.atoms[ias].species->ci[j2].lm;
                int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
                
                for (unsigned int i = 0; i < p.ndmag; i++)
                {
                    for (int j1 = 0; j1 <= j2; j1++)
                    {
                        int lm1 = geometry.atoms[ias].species->ci[j1].lm;
                        int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;
                    
                        p.L3_sum_gntyry(lm1, lm2, &p.beffrad(0, idxrf1, idxrf2, ias, i), zm(j1, j2, i));
                    }
                }
            }
            
            /* Warning!!! 
               This is the canonical code (zm should be computed for both upper an lower parts).

            for (int i = 0; i < p.nstfv; i++)
                for (int j2 = 0; j2 < sz; j2++)
                    for (int j1 = 0; j1 < sz; j1++)
                        wfb(geometry.atoms[ias].offset_wfmt + j2, i) += zm(j1,j2) * wf.scalar_wf(geometry.atoms[ias].offset_wfmt + j1, i);
            */
            // compute B_z*|wf_i>
            zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 0), sz, &wf.scalar_wf(geometry.atoms[ias].offset_wfmt, 0), 
                wf.scalar_wf.size(0), zzero, &hwf(geometry.atoms[ias].offset_wfmt, 0, 0), hwf.size(0));
            if (p.ndmag == 3)
            {
                zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 1), sz, &wf.scalar_wf(geometry.atoms[ias].offset_wfmt, 0), 
                    wf.scalar_wf.size(0), zzero, &hwf(geometry.atoms[ias].offset_wfmt, 0, 2), hwf.size(0));
                zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 2), sz, &wf.scalar_wf(geometry.atoms[ias].offset_wfmt, 0), 
                    wf.scalar_wf.size(0), -zi, &hwf(geometry.atoms[ias].offset_wfmt, 0, 2), hwf.size(0));
            }
        }
        
        std::vector<complex16> wfr(p.ngrtot);
        std::vector<complex16> zfft(p.ngrtot);
        for (unsigned int i = 0; i < p.nstfv; i++)
        {
            memset(&wfr[0], 0, p.ngrtot * sizeof(complex16));
            for (unsigned int ig = 0; ig < ngk; ig++) 
                wfr[wf.kp->idxgfft[ig]] = wf.scalar_wf(p.size_wfmt + ig, i);
                                        
            lapw_fft(1, &zfft[0]);
                       
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                zfft[ir] = wfr[ir] * p.beffir(ir, 0) * p.cfunir[ir];
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ngk; ig++) 
                hwf(p.size_wfmt + ig, i, 0) = zfft[wf.kp->idxgfft[ig]];

            if (p.ndmag == 3)
            {
                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    zfft[ir] = wfr[ir] * p.beffir(ir, 1) * p.cfunir[ir];
                                                                   
                lapw_fft(-1, &zfft[0]);
                
                for (unsigned int ig = 0; ig < ngk; ig++) 
                    hwf(p.size_wfmt + ig, i, 2) = zfft[wf.kp->idxgfft[ig]];

                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    zfft[ir] = wfr[ir] * p.beffir(ir, 2) * p.cfunir[ir];
                                                                   
                lapw_fft(-1, &zfft[0]);
                
                for (unsigned int ig = 0; ig < ngk; ig++) 
                    hwf(p.size_wfmt + ig, i, 2) -= zi * zfft[wf.kp->idxgfft[ig]];
            }
        }

        // copy -B_z|wf>
        for (unsigned int i = 0; i < p.nstfv; i++)
            for (unsigned int j = 0; j < p.size_wfmt + ngk; j++)
                hwf(j, i, 1) = -hwf(j, i, 0);

        // compute <wf_i | (h * wf_j)> for up-up block
        zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
            wf.scalar_wf.size(0), &hwf(0, 0, 0), hwf.size(0), zzero, &h(0, 0), h.size(0));
        
        // compute <wf_i | (h * wf_j)> for dn-dn block
        zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
            wf.scalar_wf.size(0), &hwf(0, 0, 1), hwf.size(0), zzero, &h(p.nstfv, p.nstfv), h.size(0));

        if (p.ndmag == 3)
        {
            // compute <wf_i | (h * wf_j)> for up-dn block
            zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
                wf.scalar_wf.size(0), &hwf(0, 0, 2), hwf.size(0), zzero, &h(0, p.nstfv), h.size(0));
        }

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


