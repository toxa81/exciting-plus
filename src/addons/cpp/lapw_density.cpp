#include "lapw.h"

void sum_uu(tensor<complex16,5>& zdens, tensor<complex16,2>& gtmp, int lm3, int ias, tensor<double,6>& densmt)
{
    int sz = geometry.atoms[ias].species->ci.size();

    for (int j2 = 0; j2 < sz; j2++)
    {
        int lm2 = geometry.atoms[ias].species->ci[j2].lm;
        int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
        
        for (int j1 = 0; j1 < sz; j1++)
        {
            int lm1 = geometry.atoms[ias].species->ci[j1].lm;
            int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;

            densmt(idxrf1, idxrf2, lm3, ias, 0, 0) += real(zdens(j1, j2, 0, 0, ias) * gtmp(lm1, lm2));
        }
    }
}

void sum_dd(tensor<complex16,5>& zdens, tensor<complex16,2>& gtmp, int lm3, int ias, tensor<double,6>& densmt)
{
    int sz = geometry.atoms[ias].species->ci.size();
    
    for (int j2 = 0; j2 < sz; j2++)
    {
        int lm2 = geometry.atoms[ias].species->ci[j2].lm;
        int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
        
        for (int j1 = 0; j1 < sz; j1++)
        {
            int lm1 = geometry.atoms[ias].species->ci[j1].lm;
            int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;

            densmt(idxrf1, idxrf2, lm3, ias, 1, 1) += real(zdens(j1, j2, 1, 1, ias) * gtmp(lm1, lm2));
        }
    }
}

void sum_ud_du(tensor<complex16,5>& zdens, tensor<complex16,2>& gtmp, int lm3, int ias, tensor<double,6>& densmt)
{
    int sz = geometry.atoms[ias].species->ci.size();
    
    for (int j2 = 0; j2 < sz; j2++)
    {
        int lm2 = geometry.atoms[ias].species->ci[j2].lm;
        int idxrf2 = geometry.atoms[ias].species->ci[j2].idxrf;
        
        for (int j1 = 0; j1 < sz; j1++)
        {
            int lm1 = geometry.atoms[ias].species->ci[j1].lm;
            int idxrf1 = geometry.atoms[ias].species->ci[j1].idxrf;

            densmt(idxrf1, idxrf2, lm3, ias, 0, 1) += 2.0 * real(zdens(j1, j2, 0, 1, ias) * gtmp(lm1, lm2));
            densmt(idxrf1, idxrf2, lm3, ias, 1, 0) -= 2.0 * imag(zdens(j1, j2, 0, 1, ias) * gtmp(lm1, lm2));
        }
    }
}

void lapw_density(lapw_wave_functions& wf, tensor<double,6>& densmt, tensor<double,3>& densir, double *occsv_)
{
    timer t("lapw_density");
    
    size_t szmax = 0;
    for (unsigned int is = 0; is < geometry.species.size(); is++)
      szmax = std::max(geometry.species[is].ci.size(), szmax);

    std::vector<int> idxocc;
    std::vector<double> occsv;
    for (unsigned int j = 0; j < p.nstsv; j++)
    {
        double d1 = occsv_[j] * wf.kp->weight;

        if (d1 > 1e-14)
        {
            idxocc.push_back(j);
            occsv.push_back(d1);
        }
    }
    
    tensor<complex16,5> zdens(szmax, szmax, p.nspinor, p.nspinor, p.natmtot);
    tensor<complex16,3> wf1(szmax, idxocc.size(), p.nspinor);
    tensor<complex16,3> wf2(szmax, idxocc.size(), p.nspinor);

    timer *t1 = new timer("lapw_density:zdens");    
    for (unsigned int ias = 0; ias < p.natmtot; ias++)
    {
        int offset = geometry.atoms[ias].offset_wfmt;
        int sz = geometry.atoms[ias].species->ci.size();

        for (unsigned i = 0; i < idxocc.size(); i++)
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                memcpy(&wf1(0, i, ispn), &wf.spinor_wf(offset, ispn, idxocc[i]), sz * sizeof(complex16));
                for (int k = 0; k < sz; k++) wf2(k, i, ispn) = wf1(k, i, ispn) * occsv[i];
            }
        
        for (unsigned int ispn2 = 0; ispn2 < p.nspinor; ispn2++)
            for (unsigned int ispn1 = 0; ispn1 < p.nspinor; ispn1++)
                if ((p.ndmag == 1 && ispn1 == ispn2) || (p.ndmag != 1))
                    zgemm<blas_worker>(0, 2, sz, sz, idxocc.size(), zone, &wf1(0, 0, ispn1), wf1.size(0), 
                        &wf2(0, 0, ispn2), wf2.size(0), zzero, &zdens(0, 0, ispn1, ispn2, ias), zdens.size(0));
    }
    delete t1;
    
    t1 = new timer("lapw_density:densmt");
#pragma omp parallel default(shared)
{
    tensor<complex16,2> gtmp(p.lmmaxapw, p.lmmaxapw);
#pragma omp for
    for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++)
    {
        for (unsigned int lm2 = 0; lm2 < p.lmmaxapw; lm2++)
            for (unsigned int lm1 = 0; lm1 < p.lmmaxapw; lm1++)
                gtmp(lm1, lm2) = p.gntyry(lm3, lm1, lm2);

        for (unsigned int ias = 0; ias < p.natmtot; ias++)
        {
            sum_uu(zdens, gtmp, lm3, ias, densmt);
            
            if (p.ndmag > 0)
                sum_dd(zdens, gtmp, lm3, ias, densmt);
            
            if (p.ndmag == 3)
                sum_ud_du(zdens, gtmp, lm3, ias, densmt);
        }
    }
}    
    delete t1;
    
    t1 = new timer("lapw_density:densir");    
    tensor<complex16,2> zfft(p.ngrtot, p.nspinor);
    for (unsigned int i = 0; i < idxocc.size(); i++)
    {
        memset(&zfft(0, 0), 0, zfft.size() * sizeof(complex16));
        for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
        {
            for (unsigned int ig = 0; ig < wf.kp->ngk; ig++)
                zfft(wf.kp->idxgfft[ig], ispn) = wf.spinor_wf(p.size_wfmt + ig, ispn, idxocc[i]);
            lapw_fft(1, &zfft(0, ispn));
        }
        
        for (unsigned int ir = 0; ir < p.ngrtot; ir++)
            densir(ir, 0, 0) += real(zfft(ir, 0) * conj(zfft(ir, 0))) * occsv[i] / geometry.omega;
        
        if (p.ndmag > 0)
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                densir(ir, 1, 1) += real(zfft(ir, 1) * conj(zfft(ir, 1))) * occsv[i] / geometry.omega;
        
        if (p.ndmag == 3)
        {
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                densir(ir, 0, 1) += 2.0 * real(zfft(ir, 0) * conj(zfft(ir, 1))) * occsv[i] / geometry.omega;
            
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                densir(ir, 1, 0) -= 2.0 * imag(zfft(ir, 0) * conj(zfft(ir, 1))) * occsv[i] / geometry.omega;
        }
    }
    delete t1;
}
