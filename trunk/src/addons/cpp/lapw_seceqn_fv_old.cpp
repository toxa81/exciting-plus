#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_fv_old)(int32_t *ikloc_, complex16 *apwalm_, complex16 *h_,
                                            complex16 *o_, double *evalfv_, complex16 *z_) 
{
    unsigned int ikloc = *ikloc_ - 1;
    unsigned int ngk = p.kpoints[ikloc].ngk;
    unsigned int msize = ngk + p.size_wfmt_lo;
    
    tensor<std::complex<double>,2> h(h_, msize, msize);
    tensor<std::complex<double>,2> o(o_, msize, msize);
    tensor<std::complex<double>,2> z(z_, p.nmatmax, p.nstfv);
    
    lapw_wave_functions wf(&p.kpoints[ikloc]);
    wf.pack_apwalm(apwalm_);

    lapw_set_h(p.kpoints[ikloc], wf.apwalm, h);
    lapw_set_o(p.kpoints[ikloc], wf.apwalm, o);
         
    return;
       
    /*tensor<std::complex<double>,2> h1;
    tensor<std::complex<double>,2> o1;
    
    std::cout << "h hash = " << h.hash() << std::endl;
    std::cout << "o hash = " << o.hash() << std::endl;
 
    if (check_evecfv) 
    {
        h1 = h;
        o1 = o;
    }

    zhegv<lapack_worker>(ldh, p.nstfv, p.evaltol, &h(0, 0), &o(0, 0), evalfv_, 
        &z(0, 0), p.nmatmax);

    if (check_evecfv) 
    {
        // use o and h as temporary arrays
        for (int i = 0; i < p.nstfv; i++)
            for (int j = 0; j < ldh; j++)
                o(j, i) = -evalfv_[i] * z(j, i);
    
        zhemm<blas_worker>(0, 0, ldh, p.nstfv, zone, &h1(0, 0), ldh, &z(0, 0), 
            p.nmatmax, zzero, &h(0, 0), ldh);
        zhemm<blas_worker>(0, 0, ldh, p.nstfv, zone, &o1(0, 0), ldh, &o(0, 0),
            ldh, zone, &h(0, 0), ldh); 
    
        double L2norm = 0.0;
        for (int i = 0; i < p.nstfv; i++)
            for (int j = 0; j < ldh; j++) L2norm += abs(h(j, i) * conj(h(j, i)));
        
        std::cout << "check evecfv :" << std::endl
                  << "  number of bands = " << p.nstfv << std::endl
                  << "  total L2 diff = " << sqrt(L2norm) << std::endl;
    }

    tensor<std::complex<double>,2> fvmt(p.wfmt_size + ngp, p.nstfv);
    lapw_fvmt(capwalm, z, fvmt);
*/
    //lapw_test_fvmt(ovlprad, fvmt, ngp, igpig);
    
    //lapw_set_sv(ngp, igpig, beffrad_, beffir_, fvmt, evalfv_);
    
    //exit(0);
}



