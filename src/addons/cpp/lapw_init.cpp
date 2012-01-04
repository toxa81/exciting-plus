#include "lapw.h"

extern "C" void FORTFUNC(lapw_init)()
{
    p.wfmt_size_apw = 0;
    p.wfmt_size_lo = 0;
    
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;
        atom->ci_apw_by_lmo = tensor<int,2>(p.lmmaxapw, p.apwordmax);
        
        atom->offset_apw = p.wfmt_size_apw;
        atom->offset_lo = p.wfmt_size_lo;

        std::vector<int> order(p.lmaxapw + 1, 0);
        int idxrf = 0;
        
        for (int l = 0; l <= p.lmaxapw; l++)
        {
            for (unsigned int io = 0; io < species->apw_descriptors[l].radial_solution_descriptors.size(); io++)
            {
                for (int m = -l, lm = l * l; m <= l; m++, lm++)
                {
                    atom->ci_apw.push_back(mtci(l, m, order[l], idxrf));
                    atom->ci.push_back(mtci(l, m, order[l], idxrf));
                    atom->ci_apw_by_lmo(lm, io) = atom->ci_apw.size() - 1; 
                }
                order[l]++;
                idxrf++;
            }
        }
        p.wfmt_size_apw += atom->ci_apw.size();
        
        for (unsigned int ilo = 0; ilo < species->lo_descriptors.size(); ilo++)
        {
            int l = species->lo_descriptors[ilo].l;
            for (int m = -l; m <= l; m++)
            { 
                atom->ci_lo.push_back(mtci(l, m, order[l], idxrf, ilo));
                atom->ci.push_back(mtci(l, m, order[l], idxrf, ilo));
            }
            order[l]++;
            idxrf++;
        }
        p.wfmt_size_lo += atom->ci_lo.size();
    }
    
    p.wfmt_size = p.wfmt_size_apw + p.wfmt_size_lo;
}

