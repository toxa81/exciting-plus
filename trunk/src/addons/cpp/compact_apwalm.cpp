#include "lapw.h"

void compact_apwalm(int ngp, 
                    std::complex<double> *apwalm_, 
                    tensor<std::complex<double>,2>& capwalm)
{
    tensor<std::complex<double>,4> apwalm(apwalm_, p.ngkmax, p.apwordmax, p.lmmaxapw, p.natmtot);
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];
        Species *species = atom->species;

        for (unsigned int j = 0; j < species->size_ci_apw; j++)
        {
            int io = species->ci[j].order;
            int lm = species->ci[j].lm;
            for (int ig = 0; ig < ngp; ig++)
                capwalm(ig, atom->offset_apw + j) = conj(apwalm(ig, io, lm, ias));
        }
    }
}


