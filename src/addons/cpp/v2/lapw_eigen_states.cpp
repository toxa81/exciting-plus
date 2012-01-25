#include "lapw.h"

void lapw_eigen_states::pack_apwalm(complex16 *apwalm_)
{   
    unsigned int ngk = kp->ngk;
    apwalm.set_dimensions(ngk, p.size_wfmt_apw);
    apwalm.allocate();
    mdarray<complex16,4> apwalm_tmp(apwalm_, p.ngkmax, p.apwordmax, p.lmmaxapw, p.natmtot);
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j = 0; j < species->size_ci_apw; j++)
        {
            int io = species->ci[j].order;
            int lm = species->ci[j].lm;
            for (unsigned int ig = 0; ig < ngk; ig++)
                apwalm(ig, atom->offset_apw + j) = conj(apwalm_tmp(ig, io, lm, ias));
        }
    }
}

void lapw_eigen_states::generate_scalar_wf()
{
    timer t("lapw_eigen_states::generate_scalar_wf");
    
    int ngk = kp->ngk;
    
    scalar_wf.set_dimensions(p.size_wfmt + ngk, p.nstfv);
    scalar_wf.allocate();
    
    zgemm<cpu>(2, 0, p.size_wfmt_apw, p.nstfv, ngk, zone, &apwalm(0, 0), ngk, &evecfv(0, 0), 
        evecfv.size(0), zzero, &scalar_wf(0, 0), scalar_wf.size(0));
    
    for (unsigned int j = 0; j < p.nstfv; j++)
    {
        // move apw blocks to proper positions
        for (unsigned int ias = p.natmtot - 1; ias > 0; ias--)
            memmove(&scalar_wf(geometry.atoms[ias].offset_wfmt, j), 
                    &scalar_wf(geometry.atoms[ias].offset_apw, j), 
                    geometry.atoms[ias].species->size_ci_apw * sizeof(complex16));
        
        // copy lo bock from the eigen-vector
        for (unsigned int ias = 0; ias < p.natmtot; ias++)
            if (geometry.atoms[ias].species->size_ci_lo > 0)
                memcpy(&scalar_wf(geometry.atoms[ias].offset_wfmt + geometry.atoms[ias].species->size_ci_apw, j), 
                       &evecfv(ngk + geometry.atoms[ias].offset_lo, j), 
                       geometry.atoms[ias].species->size_ci_lo * sizeof(complex16));
        
        // copy plane-wave block
        memcpy(&scalar_wf(p.size_wfmt, j), &evecfv(0, j), ngk * sizeof(complex16));
    }
}

void lapw_eigen_states::generate_spinor_wf(diagonalization mode)
{
    timer t("lapw_eigen_states::generate_spinor");
    
    int ngk = kp->ngk;
    spinor_wf.set_dimensions(p.size_wfmt + ngk, p.nspinor, p.nstsv);
    spinor_wf.allocate();

    for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
        zgemm<cpu>(0, 0, spinor_wf.size(0), p.nstsv, p.nstfv, zone, &scalar_wf(0, 0), scalar_wf.size(0),
            &evecsv(ispn * p.nstfv, 0), evecsv.size(0), zzero, &spinor_wf(0, ispn, 0), spinor_wf.size(0) * spinor_wf.size(1));
}

void lapw_eigen_states::test_scalar_wf(int use_fft)
{
    unsigned int ngk = kp->ngk;
    
    std::vector<complex16> zfft;
    std::vector<complex16> v1;
    std::vector<complex16> v2;
    
    if (use_fft == 1) 
    {
        v1.resize(ngk);
        zfft.resize(p.ngrtot);
    }
    
    if (use_fft == 2) 
    {
        v1.resize(p.ngrtot);
        v2.resize(p.ngrtot);
    }
    
    double maxerr = 0;

    for (unsigned int j1 = 0; j1 < p.nstfv; j1++)
    {
        if (use_fft == 1)
        {
            memset(&zfft[0], 0, p.ngrtot * sizeof(complex16));
            for (unsigned int ig = 0; ig < ngk; ig++) 
                zfft[kp->idxgfft[ig]] = scalar_wf(p.size_wfmt + ig, j1);
                
            lapw_fft(1, &zfft[0]);
            
            for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                zfft[ir] *= p.cfunir[ir];
            
            lapw_fft(-1, &zfft[0]);
            
            for (unsigned int ig = 0; ig < ngk; ig++) 
                v1[ig] = zfft[kp->idxgfft[ig]];
        }
        
        if (use_fft == 2)
        {
            memset(&v1[0], 0, p.ngrtot * sizeof(complex16));
            for (unsigned int ig = 0; ig < ngk; ig++) 
                v1[kp->idxgfft[ig]] = scalar_wf(p.size_wfmt + ig, j1);
            
            lapw_fft(1, &v1[0]);
        }
       
        for (unsigned int j2 = 0; j2 < p.nstfv; j2++)
        {
            complex16 zsum(0,0);
            for (unsigned int ias = 0; ias < p.natmtot; ias++)
            {
                int offset_wfmt = geometry.atoms[ias].offset_wfmt;
                tensor<int,2> *ci_by_lmo = &geometry.atoms[ias].species->ci_by_lmo;

                for (unsigned int l = 0; l <= p.lmaxapw; l++)
                {
                    int ordmax = geometry.atoms[ias].species->rfmt_order[l];
                    for (int io1 = 0; io1 < ordmax; io1++)
                        for (int io2 = 0; io2 < ordmax; io2++)
                            for (int m = -l; m <= (int)l; m++)
                                zsum += conj(scalar_wf(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io1), j1)) *
                                             scalar_wf(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io2), j2) * 
                                             p.ovlprad(l, io1, io2, ias);
                }
            }
            
            if (use_fft == 0) 
            {
                int iv[3];
                for (unsigned int ig1 = 0; ig1 < ngk; ig1++)
                {
                    for (unsigned int ig2 = 0; ig2 < ngk; ig2++)
                    {
                        for (int k = 0; k < 3; k++) iv[k] = p.ivg(k, kp->idxg[ig1]) - p.ivg(k, kp->idxg[ig2]);
                        int ig3 = p.ivgig(iv[0], iv[1], iv[2]);
                        zsum += conj(scalar_wf(p.size_wfmt + ig1, j1)) * scalar_wf(p.size_wfmt + ig2, j2) * p.cfunig[ig3];
                    }
               }
           }
           if (use_fft == 1)
           {
               for (unsigned int ig = 0; ig < ngk; ig++)
                   zsum += conj(v1[ig]) * scalar_wf(p.size_wfmt + ig, j2);
           }
           
           if (use_fft == 2)
           {
               memset(&v2[0], 0, p.ngrtot * sizeof(complex16));
               for (unsigned int ig = 0; ig < ngk; ig++) 
                   v2[kp->idxgfft[ig]] = scalar_wf(p.size_wfmt + ig, j2);
            
               lapw_fft(1, &v2[0]);

               for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                   zsum += conj(v1[ir]) * v2[ir] * p.cfunir[ir] / double(p.ngrtot);
           }

           zsum = (j1 == j2) ? zsum - zone : zsum;
           maxerr = std::max(maxerr, abs(zsum));
        }
    }
    std :: cout << "maximum error = " << maxerr << std::endl;
}

void lapw_eigen_states::test_spinor_wf(int use_fft)
{
    unsigned int ngk = kp->ngk;
    
    std::vector<complex16> zfft(p.ngrtot);
    std::vector<complex16> v1[p.nspinor];
    
    for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
    {
        if (use_fft == 1) v1[ispn].resize(ngk);
        if (use_fft == 2) v1[ispn].resize(p.ngrtot);
    }
    
    double maxerr = 0;

    for (unsigned int j1 = 0; j1 < p.nstsv; j1++)
    {
        if (use_fft == 1)
        {
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                memset(&zfft[0], 0, p.ngrtot * sizeof(complex16));
                for (unsigned int ig = 0; ig < ngk; ig++) 
                    zfft[kp->idxgfft[ig]] = spinor_wf(p.size_wfmt + ig, ispn, j1);
                    
                lapw_fft(1, &zfft[0]);
                
                for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                    zfft[ir] *= p.cfunir[ir];
                
                lapw_fft(-1, &zfft[0]);
                
                for (unsigned int ig = 0; ig < ngk; ig++) 
                    v1[ispn][ig] = zfft[kp->idxgfft[ig]];
            }
        }
        
        if (use_fft == 2)
        {
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                memset(&v1[ispn][0], 0, p.ngrtot * sizeof(complex16));
                for (unsigned int ig = 0; ig < ngk; ig++) 
                    v1[ispn][kp->idxgfft[ig]] = spinor_wf(p.size_wfmt + ig, ispn, j1);
                
                lapw_fft(1, &v1[ispn][0]);
            }
        }
       
        for (unsigned int j2 = 0; j2 < p.nstsv; j2++)
        {
            complex16 zsum(0,0);
            
            for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
            {
                for (unsigned int ias = 0; ias < p.natmtot; ias++)
                {
                    int offset_wfmt = geometry.atoms[ias].offset_wfmt;
                    tensor<int,2> *ci_by_lmo = &geometry.atoms[ias].species->ci_by_lmo;

                    for (unsigned int l = 0; l <= p.lmaxapw; l++)
                    {
                        int ordmax = geometry.atoms[ias].species->rfmt_order[l];
                        for (int io1 = 0; io1 < ordmax; io1++)
                            for (int io2 = 0; io2 < ordmax; io2++)
                                for (int m = -l; m <= (int)l; m++)
                                    zsum += conj(spinor_wf(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io1), ispn, j1)) *
                                                 spinor_wf(offset_wfmt + (*ci_by_lmo)(idxlm(l, m), io2), ispn, j2) * 
                                                 p.ovlprad(l, io1, io2, ias);
                    }
                }
            }

            if (use_fft == 0) 
            {
                int iv[3];
                for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
                {
                    for (unsigned int ig1 = 0; ig1 < ngk; ig1++)
                    {
                        for (unsigned int ig2 = 0; ig2 < ngk; ig2++)
                        {
                            for (int k = 0; k < 3; k++) iv[k] = p.ivg(k, kp->idxg[ig1]) - p.ivg(k, kp->idxg[ig2]);
                            int ig3 = p.ivgig(iv[0], iv[1], iv[2]);
                            zsum += conj(spinor_wf(p.size_wfmt + ig1, ispn, j1)) * spinor_wf(p.size_wfmt + ig2, ispn, j2) * p.cfunig[ig3];
                        }
                    }    
                }   
            }
            if (use_fft == 1)
            {
                for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
                {
                    for (unsigned int ig = 0; ig < ngk; ig++)
                        zsum += conj(v1[ispn][ig]) * spinor_wf(p.size_wfmt + ig, ispn, j2);
                }
            }
           
            if (use_fft == 2)
            {
                for (unsigned int ispn = 0; ispn < p.nspinor; ispn++)
                {
                    memset(&zfft[0], 0, p.ngrtot * sizeof(complex16));
                    for (unsigned int ig = 0; ig < ngk; ig++) 
                        zfft[kp->idxgfft[ig]] = spinor_wf(p.size_wfmt + ig, ispn, j2);
            
                    lapw_fft(1, &zfft[0]);

                    for (unsigned int ir = 0; ir < p.ngrtot; ir++)
                        zsum += conj(v1[ispn][ir]) * zfft[ir] * p.cfunir[ir] / double(p.ngrtot);
                }
            }

            zsum = (j1 == j2) ? zsum - zone : zsum;
            maxerr = std::max(maxerr, abs(zsum));
        }
    }
    std :: cout << "maximum error = " << maxerr << std::endl;
}
