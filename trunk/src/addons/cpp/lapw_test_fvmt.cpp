#include "lapw.h"

const int use_fft = 1;

void lapw_test_fvmt(tensor<double,4>& ovlprad, 
                    tensor<std::complex<double>,2>& fvmt,
                    int ngp,
                    tensor<int,1>& igpig)
{
    double maxerr = 0;
    std::vector< std::complex<double> > zfft;
    std::vector< std::complex<double> > v1;
    std::vector< std::complex<double> > v2;
    
    if (use_fft == 1) 
    {
        v1.resize(ngp);
        zfft.resize(p.ngrtot);
    }
    
    if (use_fft == 2) 
    {
        v1.resize(p.ngrtot);
        v2.resize(p.ngrtot);
    }

    for (int j1 = 0; j1 < p.nstfv; j1++)
    {
        if (use_fft == 1)
        {
            memset(&zfft[0], 0, p.ngrtot * sizeof(std::complex<double>));
            for (int ig = 0; ig < ngp; ig++) 
                zfft[p.igfft[igpig(ig) - 1]] = fvmt(p.wfmt_size + ig, j1);
            
            int dim = 3;
            int dir = 1;
            FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &zfft[0]);
            
            for (int ir = 0; ir < p.ngrtot; ir++)
                zfft[ir] *= p.cfunir[ir];
            
            dir = -1;
            FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &zfft[0]);
            
            for (int ig = 0; ig < ngp; ig++) 
                v1[ig] = zfft[p.igfft[igpig(ig) - 1]];
        }
        
        if (use_fft == 2)
        {
            memset(&v1[0], 0, p.ngrtot * sizeof(std::complex<double>));
            for (int ig = 0; ig < ngp; ig++) 
                v1[p.igfft[igpig(ig) - 1]] = fvmt(p.wfmt_size + ig, j1);
            
            int dim = 3;
            int dir = 1;
            FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &v1[0]);
        }
       
        for (int j2 = 0; j2 < p.nstfv; j2++)
        {
            std::complex<double> zsum(0,0);
            for (int ias = 0; ias < p.natmtot; ias++)
            {
                int offset_wfmt = geometry.atoms[ias].offset_wfmt;
                tensor<int,2> *ci_by_lmo = &geometry.atoms[ias].species->ci_by_lmo;

                for (int l = 0; l <= p.lmaxapw; l++)
                {
                    int ordmax = geometry.atoms[ias].species->rfmt_order[l];
                    for (int io1 = 0; io1 < ordmax; io1++)
                        for (int io2 = 0; io2 < ordmax; io2++)
                            for (int m = -l; m <= l; m++)
                                zsum += conj(fvmt(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io1), j1)) *
                                             fvmt(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io2), j2) * 
                                             ovlprad(l, io1, io2, ias);
                }
            }
            
            if (use_fft == 0) 
            {
                int iv[3];
                for (int ig1 = 0; ig1 < ngp; ig1++)
                {
                    for (int ig2 = 0; ig2 < ngp; ig2++)
                    {
                        for (int k = 0; k < 3; k++) iv[k] = p.ivg(k, igpig(ig1) - 1) - p.ivg(k, igpig(ig2) - 1);
                        zsum += conj(fvmt(p.wfmt_size + ig1, j1)) * fvmt(p.wfmt_size + ig2, j2) * p.cfunig[p.ivgig(iv[0], iv[1], iv[2])];
                    }
               }
           }
           if (use_fft == 1)
           {
               for (int ig = 0; ig < ngp; ig++)
                   zsum += conj(v1[ig]) * fvmt(p.wfmt_size + ig, j2);
           }
           
           if (use_fft == 2)
           {
               memset(&v2[0], 0, p.ngrtot * sizeof(std::complex<double>));
               for (int ig = 0; ig < ngp; ig++) 
                   v2[p.igfft[igpig(ig) - 1]] = fvmt(p.wfmt_size + ig, j2);
            
               int dim = 3;
               int dir = 1;
               FORTRAN(zfftifc)(&dim, &p.ngrid[0], &dir, &v2[0]);

               for (int ir = 0; ir < p.ngrtot; ir++)
                   zsum += conj(v1[ir]) * v2[ir] * p.cfunir[ir] / double(p.ngrtot);
           }

           zsum = (j1 == j2) ? zsum - zone : zsum;
           maxerr = std::max(maxerr, abs(zsum));
        }
    }
    std :: cout << "maximum error = " << maxerr << std::endl;
}
