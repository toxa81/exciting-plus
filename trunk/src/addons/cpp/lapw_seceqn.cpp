#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn)(int32_t *ikloc_, complex16 *apwalm_, complex16 *evecfv_, 
                                     double *evalfv_, complex16 *evecsv_, double *evalsv_)
{
    unsigned int ikloc = *ikloc_ - 1;
    unsigned int ngk = p.kpoints[ikloc].ngk;
    unsigned int msize = ngk + p.size_wfmt_lo;
    
    tensor<complex16,2> h(msize, msize);
    tensor<complex16,2> o(msize, msize);
    tensor<complex16,2> evecfv(evecfv_, p.nmatmax, p.nstfv);
    tensor<complex16,2> evecsv(evecsv_, p.nstsv, p.nstsv); 
    
    lapw_wave_functions wf(&p.kpoints[ikloc]);
    wf.pack_apwalm(apwalm_);

    lapw_set_h(p.kpoints[ikloc], wf.apwalm, h);
    lapw_set_o(p.kpoints[ikloc], wf.apwalm, o);
       
    tensor<complex16,2> h1;
    tensor<complex16,2> o1;
    
    if (check_evecfv) 
    {
        h1 = h;
        o1 = o;
    }

    zhegv<lapack_worker>(msize, p.nstfv, p.evaltol, &h(0, 0), &o(0, 0), evalfv_, &evecfv(0, 0), 
        evecfv.size(0));

    if (check_evecfv) 
    {
        // use o and h as temporary arrays
        for (unsigned int i = 0; i < p.nstfv; i++)
            for (unsigned int j = 0; j < msize; j++)
                o(j, i) = -evalfv_[i] * evecfv(j, i);
    
        zhemm<blas_worker>(0, 0, msize, p.nstfv, zone, &h1(0, 0), h1.size(0), &evecfv(0, 0), 
            evecfv.size(0), zzero, &h(0, 0), h.size(0));
        zhemm<blas_worker>(0, 0, msize, p.nstfv, zone, &o1(0, 0), o1.size(0), &o(0, 0),
            o.size(0), zone, &h(0, 0), h.size(0)); 
    
        for (unsigned int i = 0; i < p.nstfv; i++)
        {
            double L2norm = 0.0;
            for (unsigned int j = 0; j < msize; j++) L2norm += abs(h(j, i) * conj(h(j, i)));
            L2norm = sqrt(L2norm);
            if (L2norm > 1e-14)
                std::cout << "eigen-vector : " << i << " |H*z_i - E_i*O*z_i| = " << L2norm <<std::endl;
        }
    }

    wf.generate_scalar(evecfv);
    
    if (check_scalar_wf)
    {
        for (int i = 0; i < 3; i++)
            wf.test_scalar(i);
    }

    if (p.ndmag == 0)
    {
        memset(&evecsv(0, 0), 0, evecsv.size() * sizeof(complex16));
        for (unsigned int i = 0; i < p.nstfv; i++)
        {
            evecsv(i, i) = zone;
            evalsv_[i] = evalfv_[i];
        }
        return;
    }

    lapw_set_sv(wf, evalfv_, evecsv);
  
    if (p.ndmag == 1)
    {
        zheev<lapack_worker>(p.nstfv, &evecsv(0, 0), evecsv.size(0), evalsv_);
        zheev<lapack_worker>(p.nstfv, &evecsv(p.nstfv, p.nstfv), evecsv.size(0), &evalsv_[p.nstfv]);
    } 
    if (p.ndmag == 3)
    {
        zheev<lapack_worker>(p.nstsv, &evecsv(0, 0), evecsv.size(0), evalsv_);
    } 
}



