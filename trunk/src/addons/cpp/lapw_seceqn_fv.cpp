#include "lapw.h"

extern "C" void FORTRAN(lapw_seceqn_fv)(int32_t *ngp_,
                                         int32_t *ldh_,
                                         int32_t *ncolh,
                                         int32_t *igpig_,
                                         double *vgpc_,
                                         std::complex<double> *veffig_,
                                         std::complex<double> *apwalm_,
                                         double *apwfr_,
                                         double *apwdfr_,
                                         double *hmltrad_,
                                         double *ovlprad_,
                                         double *evalfv_,
                                         std::complex<double> *z_) 
{
    int ngp = *ngp_;
    int ldh = *ldh_;
    tensor<double,2> vgpc(vgpc_, 3, p.ngkmax);
    tensor<int,1> igpig(igpig_, p.ngkmax);
    tensor<std::complex<double>,1> veffig(veffig_, p.ngvec);
    tensor<std::complex<double>,2> h(ldh, *ncolh);
    tensor<std::complex<double>,2> o(ldh, *ncolh);
    tensor<std::complex<double>,2> z(z_, p.nmatmax, p.nstfv);
    tensor<double,5> apwfr(apwfr_, p.nrmtmax, 2, p.apwordmax, p.lmaxapw + 1, p.natmtot);
    tensor<double,3> apwdfr(apwdfr_, p.apwordmax, p.lmaxapw + 1, p.natmtot);
    tensor<double,4> hmltrad(hmltrad_, p.lmmaxvr, p.nrfmtmax, p.nrfmtmax, p.natmtot);
    tensor<double,4> ovlprad(ovlprad_, p.lmaxapw + 1, p.ordrfmtmax, p.ordrfmtmax, p.natmtot);
    
    tensor<std::complex<double>,2> capwalm(ngp, p.wfmt_size_apw);
    compact_apwalm(ngp, apwalm_, capwalm);

    //lapw_set_h(ngp, ldh, igpig, vgpc, veffig, capwalm, apwfr, apwdfr, hmltrad, h);
    //lapw_set_o(ngp, ldh, igpig, capwalm, ovlprad, o);

    tensor<std::complex<double>,2> h1;
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
    
    /*tensor<std::complex<double>,2> fvmt(p.wfmt_size + ngp, p.nstfv);
    lapw_fvmt(capwalm, z, fvmt);
    
    for (unsigned int j = 0; j < geometry.atoms[0].species->size_ci_apw; j++)
    {
        //std::cout.precision(10);
        //std::cout.fill(' ');
        //std::cout
        std::cout << "j=" << j << std::fixed << std::setprecision(10) << " wf=" << fvmt(j, 1) << std::endl;
    }*/
}



