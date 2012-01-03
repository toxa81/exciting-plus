#include "lapw.h"

void compact_apwalm(int ngp, 
                    std::complex<double> *apwalm_, 
                    tensor<std::complex<double>,2>& capwalm)
{
    tensor<std::complex<double>,4> apwalm(apwalm_, p.ngkmax, p.apwordmax, p.lmmaxapw, p.natmtot);
    for (unsigned int ias = 0; ias < geometry.atoms.size(); ias++)
    {
        Atom *atom = &geometry.atoms[ias];

        for (unsigned int j = 0; j < atom->ci_apw.size(); j++)
        {
            int io = atom->ci_apw[j].order;
            int lm = atom->ci_apw[j].lm;
            for (int ig = 0; ig < ngp; ig++)
                capwalm(ig, atom->offset_apw + j) = conj(apwalm(ig, io, lm, ias));
        }
    }
}


