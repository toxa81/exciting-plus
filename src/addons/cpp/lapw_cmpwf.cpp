#include "lapw.h"

extern "C" void FORTRAN(lapw_cmpwf)(int32_t *ngp_, 
                                     std::complex<double> *zvfv_, 
                                     std::complex<double> *apwalm_,
                                     std::complex<double> *wfmt_)
{
    int ngp = *ngp_;

    tensor<std::complex<double>,2> z(zvfv_, p.nmatmax, p.nstfv);
    
    tensor<std::complex<double>,2> capwalm(ngp, p.wfmt_size_apw);
    compact_apwalm(ngp, apwalm_, capwalm);
    
    tensor<std::complex<double>,4> wfmt_elk(wfmt_, p.lmmaxapw, p.ordrfmtmax, p.natmtot, p.nstfv);

    tensor<std::complex<double>,2> wfmt(p.wfmt_size + ngp, p.nstfv);
    lapw_fvmt(capwalm, z, wfmt);
   
    double d = 0;
    for (int ias = 0; ias < p.natmtot; ias++)
    {
        for (unsigned int j = 0; j < geometry.atoms[ias].species->ci.size(); j++)
        {
            int lm = geometry.atoms[ias].species->ci[j].lm;
            int order = geometry.atoms[ias].species->ci[j].order;
            for (int ist = 0; ist < p.nstfv; ist++)
                d += abs(wfmt(geometry.atoms[ias].offset_wfmt + j, ist) - wfmt_elk(lm, order, ias, ist));
        }
    
    }
    std::cout << "Wave-function difference : " << d << std::endl;
}  
