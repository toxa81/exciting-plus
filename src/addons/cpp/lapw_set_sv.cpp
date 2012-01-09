#include "lapw.h"

void lapw_set_sv(int ngp, 
                 tensor<int,1>& igpig, 
                 double *beffrad_, 
                 double *beffir_, 
                 tensor<std::complex<double>,2>& wfmt, 
                 double *evalfv_,
                 tensor<std::complex<double>,2>& evecsv)
{
    tensor<double,5> beffrad(beffrad_, p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot, p.ndmag);
    tensor<double,1> beffir(beffir_, p.ngrtot);

    //tensor<std::complex<double>,2> evecsv(p.nstsv, p.nstsv);
    memset(&evecsv(0, 0), 0, p.nstsv * p.nstsv * sizeof(std::complex<double>));

    if (p.spinpol)
    {
        tensor<std::complex<double>,2> wfb(p.wfmt_size + ngp, p.nstfv);
        memset(&wfb(0, 0), 0, (p.wfmt_size + ngp) * p.nstfv * sizeof(std::complex<double>));
        
        for (int ias = 0; ias < p.natmtot; ias++)
        {
            int sz = geometry.atoms[ias].species->ci.size();
            tensor<std::complex<double>,2> zm(sz, sz);
            memset(&zm(0, 0), 0, sz * sz * sizeof(std::complex<double>));
            for (int j2 = 0; j2 < sz; j2++)
            {
                int lm2 = geometry.atoms[ias].species->ci[j2].lm;
                int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
                for (int j1 = 0; j1 < sz; j1++)
                {
                    int lm1 = geometry.atoms[ias].species->ci[j1].lm;
                    int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;
                    p.L3_sum_gntyry(lm1, lm2, &beffrad(0, idxrf1, idxrf2, ias, 0), zm(j1, j2));
                }
            }

            for (int i = 0; i < p.nstfv; i++)
                for (int j2 = 0; j2 < sz; j2++)
                    for (int j1 = 0; j1 < sz; j1++)
                        wfb(geometry.atoms[ias].offset_wfmt + j2, i) += zm(j1,j2) * wfmt(geometry.atoms[ias].offset_wfmt + j1, i);
        }
        
        std::vector< std::complex<double> > zfft(p.ngrtot);
        for (int i = 0; i < p.nstfv; i++)
        {
            memset(&zfft[0], 0, p.ngrtot * sizeof(std::complex<double>));
            for (int ig = 0; ig < ngp; ig++) 
                zfft[p.igfft[igpig(ig) - 1]] = wfmt(p.wfmt_size + ig, i);
                                        
            int dim = 3;
            int dir = 1;
            FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &zfft[0]);
                       
            for (int ir = 0; ir < p.ngrtot; ir++)
              zfft[ir] *= (beffir(ir) * p.cfunir[ir]);
                                                               
            dir = -1;
            FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &zfft[0]);
            
            for (int ig = 0; ig < ngp; ig++) 
                wfb(p.wfmt_size + ig, i) = zfft[p.igfft[igpig(ig) - 1]];

        }

        for (int j1 = 0; j1 < p.nstfv; j1++)
        {
            for (int j2 = 0; j2 < p.nstfv; j2++)
            {
                int i1 = j1;
                int i2 = j2;
                
                if (i1 <= i2) 
                {
                    for (unsigned int k = 0; k < p.wfmt_size + ngp; k++)
                        evecsv(i1, i2) += wfb(k, j2) * conj(wfmt(k, j1));
                }
                
                i1 += p.nstfv;
                i2 += p.nstfv;
                if (i1 <= i2) 
                {
                    for (unsigned int k = p.wfmt_size; k < p.wfmt_size + ngp; k++)
                        evecsv(i1, i2) -= conj(wfb(k, j1)) * wfmt(k, j2);
                }
            }
        }

        int i = 0;
        for (int ispn = 0; ispn < 2; ispn++)
        {
            for (int ist = 0; ist < p.nstfv; ist++)
            {
                evecsv(i, i) += evalfv_[ist];
                i++;
            }
        }
        
        /*std::vector<double> evalsv(p.nstsv);
        zheev<lapack_worker>(p.nstsv, &evecsv(0, 0), &evalsv[0]);
        
        for (int i = 0; i < p.nstsv; i++)
        {
            std::cout << "i = " << i <<" " << evalsv[i] << std::endl;
        } */
    }
}


