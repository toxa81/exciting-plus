#include "lapw.h"

extern "C" void FORTFUNC(lapw_init)()
{
    for (unsigned int is = 0; is < geometry.species.size(); is++)
    {
        geometry.species[is].rfmt_order.resize(p.lmaxapw + 1, 0);
        
        for (int l = 0; l <= p.lmaxapw; l++)
        {
            for (unsigned int io = 0; io < geometry.species[is].apw_descriptors[l].radial_solution_descriptors.size(); io++)
            {
                for (int m = -l; m <= l; m++)
                    geometry.species[is].ci.push_back(mtci(l, m, geometry.species[is].rfmt_order[l], geometry.species[is].nrfmt));
                
                geometry.species[is].rfmt_order[l]++;
                geometry.species[is].nrfmt++;
            }
        }
        geometry.species[is].size_ci_apw = geometry.species[is].ci.size();
        
        for (unsigned int ilo = 0; ilo < geometry.species[is].lo_descriptors.size(); ilo++)
        {
            int l = geometry.species[is].lo_descriptors[ilo].l;
            for (int m = -l; m <= l; m++)
                geometry.species[is].ci.push_back(mtci(l, m, geometry.species[is].rfmt_order[l], geometry.species[is].nrfmt, ilo));
            
            geometry.species[is].rfmt_order[l]++;
            geometry.species[is].nrfmt++;
        }
        geometry.species[is].size_ci_lo = geometry.species[is].ci.size() - geometry.species[is].size_ci_apw;
        geometry.species[is].ci_lo = &geometry.species[is].ci[geometry.species[is].size_ci_apw];

        int maxorder = 0;
        for (int l = 0; l <= p.lmaxapw; l++)
            maxorder = std::max(maxorder, geometry.species[is].rfmt_order[l]);

        geometry.species[is].ci_by_lmo = tensor<int,2>(p.lmmaxapw, maxorder);

        for (unsigned int i = 0; i < geometry.species[is].ci.size(); i++)
            geometry.species[is].ci_by_lmo(geometry.species[is].ci[i].lm, geometry.species[is].ci[i].order) = i;
    }

    p.wfmt_size_apw = 0;
    p.wfmt_size_lo = 0;
    p.wfmt_size = 0;
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;
        
        atom->offset_apw = p.wfmt_size_apw;
        atom->offset_lo = p.wfmt_size_lo;
        atom->offset_wfmt = p.wfmt_size;

        p.wfmt_size_apw += species->size_ci_apw;
        p.wfmt_size_lo += species->size_ci_lo;
        p.wfmt_size += species->ci.size();
    }

    assert(p.wfmt_size == (p.wfmt_size_apw + p.wfmt_size_lo));
}

