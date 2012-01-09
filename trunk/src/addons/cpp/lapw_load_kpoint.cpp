#include "lapw.h"

extern "C" void FORTRAN(lapw_load_kpoint)(int32_t *ngp_, int32_t *igpig_, double *vgpc_)
{
    int ngp = *ngp_;
    tensor<double,2> vgpc(vgpc_, 3, p.ngkmax);
    tensor<int,1> igpig(igpig_, p.ngkmax);

     p.kpoints.push_back(kpoint(ngp));
     for (int ig = 0; ig < ngp; ig++)
     {
         for (int k = 0; k < 3; k++)
             p.kpoints.back().vgkc(k, ig) = vgpc(k, ig);

         p.kpoints.back().idxg[ig] = igpig(ig) - 1;
         p.kpoints.back().idxgfft[ig] = p.igfft[igpig(ig) - 1];
     }
}
