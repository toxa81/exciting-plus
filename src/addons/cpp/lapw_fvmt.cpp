#include "lapw.h"

void lapw_fvmt(tensor<std::complex<double>,2>& capwalm,
               tensor<std::complex<double>,2>& zfv,
               tensor<std::complex<double>,2>& fvmt)
{
    int ngp = capwalm.size(0);
    zgemm<gemm_worker>(2, 0, p.wfmt_size_apw, p.nstfv, ngp, zone, &capwalm(0, 0), 
       capwalm.size(0), &zfv(0, 0), zfv.size(0), zzero, &fvmt(0, 0), fvmt.size(0));
    
    for (int j = 0; j < p.nstfv; j++)
    {
        for (int ias = p.natmtot - 1; ias > 0; ias--)
            memmove(&fvmt(geometry.atoms[ias].offset_wfmt, j), 
                    &fvmt(geometry.atoms[ias].offset_apw, j), 
                    geometry.atoms[ias].species->size_ci_apw * sizeof(std::complex<double>));
        
        for (int ias = 0; ias < p.natmtot; ias++)
            if (geometry.atoms[ias].species->size_ci_lo > 0)
                memcpy(&fvmt(geometry.atoms[ias].offset_wfmt + geometry.atoms[ias].species->size_ci_apw, j), 
                       &zfv(ngp + geometry.atoms[ias].offset_lo, j), 
                       geometry.atoms[ias].species->size_ci_lo * sizeof(std::complex<double>));

        memcpy(&fvmt(p.wfmt_size, j), &zfv(0, j), ngp * sizeof(std::complex<double>));
    }
}

