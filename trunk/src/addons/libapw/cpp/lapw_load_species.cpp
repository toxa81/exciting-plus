#include "lapw.h"

extern "C" void FORTRAN(lapw_load_species)(int *nlorb, int *lorbl, int *apword, double* rmt, int *nrmt, int *llu)
{
    Species *sp = new Species();
    
    sp->rmt = *rmt;
    sp->nrmt = *nrmt;
    sp->lu = *llu;
    
    for (int i = 0; i < *nlorb; i++) 
        sp->lo_descriptors.push_back(radial_l_descriptor(lorbl[i]));
    
    for (unsigned int l = 0; l <= lapw_global.lmaxapw; l++) 
    {
        radial_l_descriptor lch(l);
        for (int io = 0; io < apword[l]; io++)
            lch.radial_solution_descriptors.push_back(radial_solution_descriptor());
        
        sp->apw_descriptors.push_back(lch);
    }
    
    lapw_global.species.push_back(sp);
}


