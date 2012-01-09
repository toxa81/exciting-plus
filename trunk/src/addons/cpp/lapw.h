#ifndef __LAPW_H__
#define __LAPW_H__

#include <vector>
#include <string>
#include <map>
#include <complex>
#include <algorithm>
#include <iomanip>
#include "tensor.h"
#include "config.h"
#include "linalg.h"

/*
    inline functions
*/

inline int idxlm(int l, int m)
{
    return l * l + l + m;
}

struct atomic_level 
{
    int n;
    int l;
    int k;
    int occupancy;
};

struct radial_solution_descriptor
{
    int n;
    int l;
    int dme;
    double enu;
    bool auto_enu;
};

class radial_l_channel_descriptor
{
    public:
        radial_l_channel_descriptor()
        {
        }
        
        radial_l_channel_descriptor(int l) : l(l)
        {
        }
    
        int l;
        std::vector<radial_solution_descriptor> radial_solution_descriptors;
};

// muffin-tin combined indices
class mtci
{
    public:
        mtci(int l, int m, int order, int idxrf) : l(l), m(m), order(order), idxrf(idxrf), idxlo(-1)
        {
            lm = idxlm(l, m);
        }
        mtci(int l, int m, int order, int idxrf, int idxlo) : l(l), m(m), order(order), idxrf(idxrf), idxlo(idxlo)
        {
            lm = idxlm(l, m);
        }
        
        int l;
        int m;
        int lm;
        int order;
        int idxrf;
        int idxlo;
};

class Species 
{
    public:

        Species() : size_ci_lo(0), size_ci_apw(0), nrfmt(0)
        {
        }
    
        Species(const std::string& symbol) : symbol(symbol) 
        {
        };
        
        std::string name;
        std::string symbol;
        int number;
        double mass;
        double rmin;
        double rmax;
        double rmt;
        int nrmt;
  
        std::vector<atomic_level> core;
        std::vector<radial_l_channel_descriptor> lo_descriptors;
        std::vector<radial_l_channel_descriptor> apw_descriptors;
        std::vector<mtci> ci;
        mtci *ci_lo;
        unsigned int size_ci_lo;
        unsigned int size_ci_apw;
        tensor<int,2> ci_by_lmo;
        std::vector<int> rfmt_order;
        unsigned int nrfmt;
};

class Atom 
{
    public:
        
        Atom()
        {
        }

        Atom(Species *species) : species(species)
        {
            
        }
        
        double posl[3];
        double posc[3];
        double bfcmt[3];
        int symclass;
        Species *species;
        unsigned int offset_apw;
        unsigned int offset_lo;
        unsigned int offset_wfmt;
};

class kpoint
{
    public:

        kpoint(unsigned int ngk) : ngk(ngk)
        {
            vgkc = tensor<double,2>(3, ngk);
            idxg.resize(ngk);
            idxgfft.resize(ngk);
        }
        
        unsigned int ngk;
    
        tensor<double,2> vgkc;
        std::vector<int> idxg;
        std::vector<int> idxgfft;
};

class Geometry 
{
    public:
    
        std::vector< std::vector<double> > avec;
        double avec_m[9];
        double ainv_m[9];
        std::vector< std::vector<double> > bvec;
        double bvec_m[9];
        double binv_m[9];
        std::vector<Species> species;
        std::map<std::string,Species*> species_by_symbol;
        std::vector<Atom> atoms;
};

class Parameters
{
    public:
        unsigned int ngkmax;
        unsigned int apwordmax;
        unsigned int lmmaxapw;
        unsigned int natmtot;
        unsigned int nspecies;
        unsigned int lmaxvr;
        unsigned int lmmaxvr;
        unsigned int lmaxapw;
        unsigned int ngvec;
        unsigned int ngrtot;
        unsigned int nlomax;
        unsigned int nrmtmax;
        unsigned int nstfv;
        unsigned int nstsv;
        unsigned int nmatmax;
        unsigned int nrfmtmax;
        unsigned int ordrfmtmax;
        double evaltol;
        int ngrid[3];
        bool spinpol;
        unsigned int ndmag;
        std::vector<int> igfft;
        std::vector< std::complex<double> > cfunig;
        std::vector<double> cfunir;
        tensor<int,2> intgv;
        tensor<int,2> ivg;
        tensor<int,3> ivgig;
        tensor<std::complex<double>,3> gntyry;
        tensor<std::vector<int>,2> L3_gntyry;
        std::vector<kpoint> kpoints;
        unsigned int wfmt_size_apw;
        unsigned int wfmt_size_lo;
        unsigned int wfmt_size;
        tensor<double,4> hmltrad;
        tensor<double,4> ovlprad;
        tensor<double,5> beffrad;
        tensor<double,5> apwfr;
        tensor<double,3> apwdfr;
        tensor<double,1> beffir;
        tensor<complex16,1> veffig;

        inline void L3_sum_gntyry(int lm1, int lm2, double *v, std::complex<double>& zsum)
        {
            for (unsigned int k = 0; k < L3_gntyry(lm1, lm2).size(); k++)
            {
                int lm3 = L3_gntyry(lm1, lm2)[k];
                zsum += gntyry(lm3, lm1, lm2) * v[lm3];
            }
        }
};

extern Geometry geometry;

extern Parameters p;

extern std::complex<double> zone;

extern std::complex<double> zzero;

/*
    forward declarations
*/

void compact_apwalm(int ngp_, complex16 *apwalm_, tensor<complex16,2>& capwalm);

void lapw_set_h(kpoint& kp, tensor<complex16,2>& capwalm, tensor<complex16,2>& h);

void lapw_set_o(kpoint& kp, tensor<complex16,2>& capwalm, tensor<complex16,2>& o);

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
                                         std::complex<double> *z_); 

void lapw_fvmt(tensor<std::complex<double>,2>& capwalm,
               tensor<std::complex<double>,2>& zfv,
               tensor<std::complex<double>,2>& fvmt);

void lapw_test_fvmt(tensor<double,4>& ovlprad,
                    tensor<std::complex<double>,2>& fvmt,
                    int ngp,
                    tensor<int,1>& igpig);

void lapw_set_sv(int ngp, 
                 tensor<int,1>& igpig, 
                 double *beffrad_, 
                 double *beffir_, 
                 tensor<std::complex<double>,2>& wfmt, 
                 double *evalfv_,
                 tensor<std::complex<double>,2>& hsv);

extern "C" void FORTRAN(zfftifc)(int *dim, 
                                  int *ngrid, 
                                  int *dir, 
                                  std::complex<double> *data);

#endif // __LAPW_H__
