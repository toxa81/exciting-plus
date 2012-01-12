#include "lapw.h"

void b_dot_wf(lapw_wave_functions& wf, tensor<complex16,3>& hwf)
{
    timer t("b_dot_wf");
    
    unsigned int ngk = wf.kp->ngk;
    
    size_t szmax = 0;
    for (unsigned int is = 0; is < geometry.species.size(); is++)
      szmax = std::max(geometry.species[is].ci.size(), szmax);
    
    tensor<complex16,3> zm(szmax, szmax, p.ndmag);
            
    for (unsigned int ias = 0; ias < p.natmtot; ias++)
    {
        int offset = geometry.atoms[ias].offset_wfmt;
        int sz = geometry.atoms[ias].species->ci.size();
        
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
        // compute hwf = hwf + B_z*|wf_j>
        zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 0), zm.size(0), 
            &wf.scalar_wf(offset, 0), wf.scalar_wf.size(0), zone, &hwf(offset, 0, 0), hwf.size(0));
        
        // compute hwf = hwf + (B_x - iB_y)|wf_j>
        if (p.ndmag == 3)
        {
            for (int j2 = 0; j2 < sz; j2++)
            {
                for (int j1 = 0; j1 <= j2; j1++)
                    zm(j1, j2, 0) = zm(j1, j2, 1) - zi * zm(j1, j2, 2);
                
                for (int j1 = j2 + 1; j1 < sz; j1++)
                    zm(j1, j2, 0) = conj(zm(j2, j1, 1)) - zi * conj(zm(j2, j1, 2));
            }
              
            zgemm<blas_worker>(0, 0, sz, p.nstfv, sz, zone, &zm(0, 0, 0), zm.size(0), 
                &wf.scalar_wf(offset, 0), wf.scalar_wf.size(0), zone, &hwf(offset, 0, 2), hwf.size(0));
            
            //zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 1), zm.size(0), 
            //    &wf.scalar_wf(offset, 0), wf.scalar_wf.size(0), zone, &hwf(offset, 0, 2), hwf.size(0));
           
            //zhemm<blas_worker>(0, 0, sz, p.nstfv, zone, &zm(0, 0, 2), zm.size, 
            //    &wf.scalar_wf(offset, 0), wf.scalar_wf.size(0), -zi, &hwf(offset, 0, 2), hwf.size(0));
        }
    }
        
    std::vector<complex16> wfr(p.ngrtot);
    std::vector<complex16> zfft(p.ngrtot);
    for (unsigned int i = 0; i < p.nstfv; i++)
    {
        memset(&wfr[0], 0, p.ngrtot * sizeof(complex16));
        for (unsigned int ig = 0; ig < ngk; ig++) 
            wfr[wf.kp->idxgfft[ig]] = wf.scalar_wf(p.size_wfmt + ig, i);
                                    
        lapw_fft(1, &wfr[0]);
                   
        for (unsigned int ir = 0; ir < p.ngrtot; ir++)
            zfft[ir] = wfr[ir] * p.beffir(ir, 0) * p.cfunir[ir];
                                                           
        lapw_fft(-1, &zfft[0]);
        
        for (unsigned int ig = 0; ig < ngk; ig++) 
            hwf(p.size_wfmt + ig, i, 0) += zfft[wf.kp->idxgfft[ig]];

        if (p.ndmag == 3)
        {
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                zfft[ir] = wfr[ir] * (p.beffir(ir, 1) - zi * p.beffir(ir, 2)) * p.cfunir[ir];
                                                               
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ngk; ig++) 
                hwf(p.size_wfmt + ig, i, 2) += zfft[wf.kp->idxgfft[ig]];
        }
    }

    // copy -B_z|wf>
    for (unsigned int i = 0; i < p.nstfv; i++)
        for (unsigned int j = 0; j < p.size_wfmt + ngk; j++)
            hwf(j, i, 1) = -hwf(j, i, 0);
}

void lapw_set_sv(lapw_wave_functions& wf, double *evalfv_, tensor<complex16,2>& h)
{
    timer t("lapw_set_sv");
    
    unsigned int ngk = wf.kp->ngk;

    memset(&h(0, 0), 0, h.size() * sizeof(complex16));

    int nhwf;
    if (p.ndmag == 0) nhwf = 1; // have only one block, nonmagnetic
    if (p.ndmag == 1) nhwf = 2; // have up-up and dn-dn blocks, collinear
    if (p.ndmag == 3) nhwf = 3; // have up-up, dn-dn and up-dn blocks, general case

    // product of the second-variational hamiltonian and a wave-function
    tensor<complex16,3> hwf(p.size_wfmt + ngk, p.nstfv, nhwf);
    memset(&hwf(0, 0, 0), 0, hwf.size() * sizeof(complex16));

    // compute product of magnetic field and wave-function 
    if (p.ndmag > 0)
        b_dot_wf(wf, hwf);
        
    // compute <wf_i | (h * wf_j)> for up-up block
    zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
        wf.scalar_wf.size(0), &hwf(0, 0, 0), hwf.size(0), zzero, &h(0, 0), h.size(0));
        
    // compute <wf_i | (h * wf_j)> for dn-dn block
    if (p.ndmag != 0)
        zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
            wf.scalar_wf.size(0), &hwf(0, 0, 1), hwf.size(0), zzero, &h(p.nstfv, p.nstfv), h.size(0));

    // compute <wf_i | (h * wf_j)> for up-dn block
    if (p.ndmag == 3)
        zgemm<blas_worker>(2, 0, p.nstfv, p.nstfv, wf.scalar_wf.size(0), zone, &wf.scalar_wf(0, 0),
          wf.scalar_wf.size(0), &hwf(0, 0, 2), hwf.size(0), zzero, &h(0, p.nstfv), h.size(0));

    unsigned int nspn = (p.ndmag == 0) ? 1 : 2;
    for (unsigned int ispn = 0, i = 0; ispn < nspn; ispn++)
        for (unsigned int ist = 0; ist < p.nstfv; ist++, i++)
                h(i, i) += evalfv_[ist];
}


