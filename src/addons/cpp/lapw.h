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
#include "timer.h"

/*
    forward class declarations
*/
struct atomic_level;
struct radial_solution_descriptor;
class radial_l_channel_descriptor;
class mtci;
class Species;
class Atom;
class kpoint;
class Geometry;
class Parameters;
class lapw_wave_functions;

/*
    function declarations
*/
void lapw_set_h(kpoint& kp, tensor<complex16,2>& apwalm, tensor<complex16,2>& h);

void lapw_set_o(kpoint& kp, tensor<complex16,2>& apwalm, tensor<complex16,2>& o);

void lapw_set_sv(lapw_wave_functions& wf, double *evalfv_, tensor<complex16,2>& h);

void lapw_fft(int32_t direction, complex16 *data);

inline int idxlm(int l, int m);

/*
    actual class definitions
*/
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
        tensor<std::vector<complex16>,2> L3_gntyry_data;
        std::vector<kpoint> kpoints;
        unsigned int size_wfmt_apw;
        unsigned int size_wfmt_lo;
        unsigned int size_wfmt;
        tensor<double,4> hmltrad;
        tensor<double,4> ovlprad;
        tensor<double,5> beffrad;
        tensor<double,5> apwfr;
        tensor<double,3> apwdfr;
        tensor<double,2> beffir;
        tensor<complex16,1> veffig;

        inline void L3_sum_gntyry(int lm1, int lm2, double *v, std::complex<double>& zsum)
        {
            for (unsigned int k = 0; k < L3_gntyry(lm1, lm2).size(); k++)
            {
                int lm3 = L3_gntyry(lm1, lm2)[k];
                zsum += gntyry(lm3, lm1, lm2) * v[lm3];
                //zsum += L3_gntyry_data(lm1, lm2)[k] * v[lm3];
            }
            //for (unsigned int lm3 = 0; lm3 < p.lmmaxvr; lm3++)
            //    zsum += gntyry(lm3, lm1, lm2) * v[lm3];
        }
};

class lapw_wave_functions
{
    public:

        lapw_wave_functions(kpoint *kp) : kp(kp)
        {
        }
                
        void pack_apwalm(complex16 *apwalm_);
        
        void generate_scalar(tensor<complex16,2>& evecfv); 
        
        void test_scalar(int use_fft);

        kpoint *kp;
        tensor<complex16,2> apwalm;
        tensor<complex16,2> scalar_wf;

};

/*
    global variables
*/
extern Geometry geometry;
extern Parameters p;
extern complex16 zone;
extern complex16 zzero;
extern complex16 zi;

/*
    inline functions
*/
inline int idxlm(int l, int m)
{
    return l * l + l + m;
}

inline int idxG12(const kpoint& kp, int ig1, int ig2)
{
    //int iv[3];
    //for (int i = 0; i < 3; i++) iv[i] = p.ivg(i, kp->idxg[ig1]) - p.ivg(i, kp->idxg[ig2]);
    return p.ivgig(p.ivg(0, kp.idxg[ig1]) - p.ivg(0, kp.idxg[ig2]), 
                   p.ivg(1, kp.idxg[ig1]) - p.ivg(1, kp.idxg[ig2]),
                   p.ivg(2, kp.idxg[ig1]) - p.ivg(2, kp.idxg[ig2]));
}

#endif // __LAPW_H__
