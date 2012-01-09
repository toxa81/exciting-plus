#include "lapw.h"

extern "C" void FORTRAN(lapw_cmpbf)(int *ias_,
                                     double *beffrad_,
                                     std::complex<double> *zb_)
{
    int ias = *ias_ - 1;
    tensor<std::complex<double>,4> zb(zb_, p.lmmaxapw, p.lmmaxapw, p.ordrfmtmax, p.ordrfmtmax);
    tensor<double,5> beffrad(beffrad_, p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot, p.ndmag);

    int sz = geometry.atoms[ias].species->ci.size();
    tensor<std::complex<double>,2> zm(sz, sz);
    memset(&zm(0, 0), 0, sz * sz * sizeof(std::complex<double>));
    double d = 0;
    for (int j2 = 0; j2 < sz; j2++)
    {
        int lm2 = geometry.atoms[ias].species->ci[j2].lm;
        int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
        int order2 = geometry.atoms[ias].species->ci[j2].order;
        for (int j1 = 0; j1 < sz; j1++)
        {
            int lm1 = geometry.atoms[ias].species->ci[j1].lm;
            int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;
            int order1 = geometry.atoms[ias].species->ci[j1].order;
            p.L3_sum_gntyry(lm1, lm2, &beffrad(0, idxrf1, idxrf2, ias, 0), zm(j1, j2));
            d += abs(zm(j1, j2) - zb(lm1, lm2, order1, order2));
        }
    }
    std::cout << "atom : " << ias << " b diff = " << d << std::endl;
}
