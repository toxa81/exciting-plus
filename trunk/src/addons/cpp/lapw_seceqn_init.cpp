#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_init)(double *hmltrad_, double *ovlprad_, double *beffrad_,
                                          double *apwfr_, double *apwdfr_, double *beffir_,
                                          complex16 *veffig_)
                                     
{
    p.hmltrad = tensor<double,4>(hmltrad_, p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot);
    p.ovlprad = tensor<double,4>(ovlprad_, p.lmaxapw + 1, p.ordrfmtmax, p.ordrfmtmax, p.natmtot);
    p.apwfr = tensor<double,5>(apwfr_, p.nrmtmax, 2, p.apwordmax, p.lmaxapw + 1, p.natmtot);
    p.apwdfr = tensor<double,3>(apwdfr_, p.apwordmax, p.lmaxapw + 1, p.natmtot);
    p.veffig = tensor<complex16,1>(veffig_, p.ngvec);
    if (p.ndmag > 0) 
    {
        p.beffrad = tensor<double,5>(beffrad_, p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot, p.ndmag);
        p.beffir = tensor<double,2>(beffir_, p.ngrtot, p.ndmag);
        
        if (lapw_diag == full)
        {
            p.beffig.set_dimensions(p.ngvec, p.ndmag);
            p.beffig.allocate();
            std::vector<complex16> zfft(p.ngrtot);
            for (unsigned int i = 0; i < p.ndmag; i++)
            {
                for (unsigned int ir = 0; ir < p.ngrtot; ir++) zfft[ir] = zone * p.beffir(ir, i) * p.cfunir[ir];
                lapw_fft(-1, &zfft[0]);
                for (unsigned int ig = 0; ig < p.ngvec; ig++) p.beffig(ig, i) = zfft[p.igfft[ig]];
            }
        }
    }
}
