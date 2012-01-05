#ifndef __LAPW_H__
#define __LAPW_H__

#include <vector>
#include <string>
#include <map>
#include <complex>
#include <algorithm>
#include "tensor.h"
#include "config.h"
#include "linalg.h"

/*
    forward declarations
*/
void compact_apwalm(int ngp, 
                    std::complex<double> *apwalm_, 
                    tensor<std::complex<double>,2>& capwalm);

void lapw_set_h(int ngp,
                int ldh,
                tensor<int,1>& igpig,
                tensor<double,2>& vgpc,
                tensor<std::complex<double>,1>& veffig,
                tensor<std::complex<double>,2>& capwalm,
                tensor<double,5>& apwfr,
                tensor<double,3>& apwdfr,
                tensor<double,4>& hmltrad,
                tensor<std::complex<double>,2>& h);

void lapw_set_o(int ngp,
                int ldo,
                tensor<int,1>& igpig,
                tensor<std::complex<double>,2>& capwalm,
                tensor<double,3>& ovlprad,
                tensor<std::complex<double>,2>& o);

void lapw_fvmt(tensor<std::complex<double>,2>& capwalm,
               tensor<std::complex<double>,2>& zfv,
               tensor<std::complex<double>,2>& fvmt);

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
        int ngkmax;
        int apwordmax;
        int lmmaxapw;
        int natmtot;
        int nspecies;
        int lmaxvr;
        int lmmaxvr;
        int lmaxapw;
        int ngvec;
        int ngrtot;
        int nlomax;
        int nrmtmax;
        int nstfv;
        int nstsv;
        int nmatmax;
        int nrfmtmax;
        int ordrfmtmax;
        double evaltol;
        
        std::vector< std::complex<double> > cfunig;
        tensor<int,2> intgv;
        tensor<int,2> ivg;
        tensor<int,3> ivgig;
        tensor<std::complex<double>,3> gntyry;
        tensor<std::vector<int>,2> L3_gntyry;
        
        unsigned int wfmt_size_apw;
        unsigned int wfmt_size_lo;
        unsigned int wfmt_size;

};

extern Geometry geometry;
extern Parameters p;
extern std::complex<double> zone;
extern std::complex<double> zzero;

#endif // __LAPW_H__
