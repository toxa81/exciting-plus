#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <mra/operator.h>
#include <world/world.h>

extern "C" void elk_init_(void);
extern "C" void elk_wan_val_(int*, int*, double*, double*, double*);
extern "C" void elk_xc_(double*, double*, double*);
extern "C" void elk_load_wann_unk_(int* ik);
extern "C" void elk_wann_unk_val_(int* ik, int* n, int* ispn, double* r_cutoff,
  double* vrc,double* val);

using namespace madness;

//***************************************************************************
template <typename Q, int NDIM>
struct abs_square_op
{
  typedef typename TensorTypeData<Q>::scalar_type resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, resultT d = abs(*_p0); *_p1 = d*d);
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};
//***************************************************************************

//***************************************************************************
template<typename Q, int NDIM>
Function<typename TensorTypeData<Q>::scalar_type,NDIM> abs_square(const Function<Q,NDIM>& func)
{
  return unary_op(func, abs_square_op<Q,NDIM>());
}


World* world;
std::vector<real_function_3d> wan;
std::vector<complex_function_3d> wann_unk;
std::vector<real_function_3d> rho_wan;

template <typename T>                                                                                   
static void xc_e_op(const Key<3>& key, Tensor<T>& t) {
  UNARY_OPTIMIZED_ITERATOR(T, t, double dens=(*_p0)*(*_p0); double exc; double vxc; elk_xc_(&dens,&vxc,&exc); *_p0=exc);                                                    
}

real_function_3d xc_energy(const real_function_3d& rho) {
  real_function_3d result=copy(rho);
  result.reconstruct();
  result.unaryop(&xc_e_op<double>);
  return result;
}  

template <typename T>                                                                                   
static void xc_v_op(const Key<3>& key, Tensor<T>& t) {
  UNARY_OPTIMIZED_ITERATOR(T, t, double dens=(*_p0)*(*_p0); double exc; double vxc; elk_xc_(&dens,&vxc,&exc); *_p0=vxc);                                                    
}

real_function_3d xc_potential(const real_function_3d& rho) {
  real_function_3d result=copy(rho);
  result.reconstruct();
  result.unaryop(&xc_v_op<double>);
  return result;
}  

class Wannier: public FunctionFunctorInterface<double,3> {
public:
  int _n;
  Vector<int,3> _vtl;

  Wannier(int n, Vector<int,3>& vtl) : _n(n), _vtl(vtl) {}

  double operator()(const coord_3d& x) const
  {
    double val[2];
    val[0]=val[1]=0.0;
    double vrc[3];
    vrc[0]=x[0]; vrc[1]=x[1]; vrc[2]=x[2];
    double r_cutoff=10.0;
    int ispn=1;
    int n=_n;
    
    elk_wan_val_(&n,&ispn,&r_cutoff,&vrc[0],&val[0]);
    return val[0];
  }
};

class WannierExc: public FunctionFunctorInterface<double,3> {
public:
  int _n;
  real_function_3d* _fptr;

  WannierExc(real_function_3d* fptr) : _fptr(fptr) {}

  double operator()(const coord_3d& x) const {
    double dens[2];
    double vxc[2];
    double exc;
    
    coord_3d p=x;
    
    dens[0]=0.0;
    double w=((real_function_3d)(*_fptr))(p);
    //dens[0]=w*w;
    //elk_xc_(&dens[0], &vxc[0], &exc);
    exc=dens[0];
    return exc;
  }
};

class Wannier_unk: public FunctionFunctorInterface<double_complex,3> {
public:
  int _n;
  int _ik;
  int _ispn;

  Wannier_unk(int n, int ik, int ispn) : _n(n), _ik(ik), _ispn(ispn) {}

  double_complex operator()(const coord_3d& x) const {
    double val[2];
    val[0]=val[1]=0.0;
    double vrc[3];
    vrc[0]=x[0]; vrc[1]=x[1]; vrc[2]=x[2];
    double r_cutoff=10.0;
    int ispn=1; //_ispn;
    int n=_n;
    int ik=_ik;
    
    elk_wann_unk_val_(&ik,&n,&ispn,&r_cutoff,&vrc[0],&val[0]);
    return double_complex(val[0], val[1]);
  }
};








extern "C" void madness_genpot_(void) {
  // defaults
  int k = 6;
  double thresh = 0.01;
  double L = 30.0;
  double ha2ev = 27.21138386;
    
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(5);
//  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
  
  real_convolution_3d op = CoulombOperator(*world, 0.001, 1e-6);
  
  Vector<int,3> vtl(0);
  
  // init Elk 
  //elk_init_();
  
  for (int n=1;n<=4;n++) {
    real_function_3d w_ = real_factory_3d(*world).functor(real_functor_3d(new Wannier(n, vtl)));
    w_.truncate();
    wan.push_back(w_);
    if (world->rank() == 0) printf("n= %i resolved\n",n);
  }

  for (int n=1;n<=4;n++) {
    real_function_3d w2n = square(wan[n-1]).truncate();
    real_function_3d hpot = op(w2n).truncate();
    for (int m=1;m<=4;m++) {
      double prod = inner(wan[n-1],wan[m-1]);
      real_function_3d w2m = square(wan[m-1]).truncate();
      double U = inner(hpot,w2m);
      if (world->rank() == 0) printf(" n= %i  m= %i  prod= %15.8f  U= %15.8f (Ha)  %15.8f (eV)\n", 
        n, m, prod, U, U*ha2ev);
    }
  }
  print("");
  for (int n=1;n<=4;n++) {
    real_function_3d excn = xc_energy(wan[n-1]);
    excn.truncate();
    real_function_3d wexcn = wan[n-1]*excn;
    for (int m=1;m<=4;m++) {
      double Exc = inner(wexcn,wan[m-1]);
      if (world->rank() == 0) printf(" n= %i  m= %i Exc= %15.8f (Ha)  %15.8f (eV)\n", 
        n, m, Exc, Exc*ha2ev);
    }
  }
  print("");
  for (int n=1;n<=4;n++) {
    real_function_3d vxcn = xc_potential(wan[n-1]);
    vxcn.truncate();
    real_function_3d wvxcn = wan[n-1]*vxcn;
    for (int m=1;m<=4;m++) {
      double Vxc = inner(wvxcn,wan[m-1]);
      if (world->rank() == 0) printf(" n= %i  m= %i Vxc= %15.8f (Ha)  %15.8f (eV)\n", 
        n, m, Vxc, Vxc*ha2ev);
    }
  }

}

extern "C" void madness_init_(void) {
  initialize(0, 0);
  world = new World(MPI::COMM_WORLD);
  startup(*world, 0, 0);
}

extern "C" void madness_init_box_(void) {
  int k = 6;
  double thresh = 0.1;
  double L = 30.0;

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(5);
//  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
}

extern "C" void madness_resolve_wannier_(int* nwantot, int* nkptnr, int* nspinor) {
  rho_wan.clear(); //clear list of WF densities
  for (int n=1;n<=*nwantot;n++) {
    real_function_3d rho_=real_factory_3d(*world);
    rho_wan.push_back(rho_);
  }
  printf(" nwantot : %i\n",*nwantot);
  printf(" nkptnr  : %i\n",*nkptnr);
  printf(" nspinor : %i\n",*nspinor);
  
  for (int ik=1;ik<=*nkptnr;ik++) {
    elk_load_wann_unk_(&ik);
    for (int n=1;n<=*nwantot;n++) {
      for (int ispn=1;ispn<=*nspinor;ispn++) {
        complex_function_3d unk_ = complex_factory_3d(*world).functor(
          complex_functor_3d(new Wannier_unk(n,ik,ispn))
        );
	real_function_3d rho_=abs_square(unk_);
	rho_.truncate();
	rho_wan[n-1]+=rho_;
      }
    }
    world->gop.fence();
    if (world->rank() == 0) printf("ik= %i resolved\n",ik);
  }
  double norm=rho_wan[0].trace();
  print("norm : %15.8f\n",norm);
  
}

extern "C" void madness_getpot_(double* pot) {
  *pot=wan[0](0.0);
}



// int aaamain(int argn,char** argv) {
//   initialize(argn, argv);
//   World world(MPI::COMM_WORLD);
//   if (world.rank() == 0)
//   {
//     print("");
//     print("--------------------------------------------");
//     print("   MADNESS", " multiresolution testsuite");
//     print("--------------------------------------------");
//     print("");
//     print("   number of processors ...", world.size());
//     print("    processor frequency ...", cpu_frequency());
//     print("          configured by ...", MADNESS_CONFIGURATION_USER);
//     print("          configured on ...", MADNESS_CONFIGURATION_HOST);
//     print("          configured at ...", MADNESS_CONFIGURATION_DATE);
//     print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
//     print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
// #ifdef WORLD_WATCHDOG
//     print("               watchdog ...", WATCHDOG_BARK_INTERVAL,
//         WATCHDOG_TIMEOUT);
// #endif
// #ifdef OPTERON_TUNE
//     print("             tuning for ...", "opteron");
// #elif defined(CORE_DUO_TUNE)
//     print("             tuning for ...", "core duo");
// #else
//     print("             tuning for ...", "core2");
// #endif
// #ifdef BOUNDS_CHECKING
//     print(" tensor bounds checking ...", "enabled");
// #endif
// #ifdef TENSOR_INSTANCE_COUNT
//     print("  tensor instance count ...", "enabled");
// #endif
//     print(" ");
//   }
// 
//   try
//   {
//     printf("Starting up the world ... \n");
// 
//     startup(world,argn,argv);
//     if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
//     hartree_potential(world);
//   }
//   catch (const MPI::Exception& e)
//   {
//     //        print(e);
//     error("caught an MPI exception");
//   }
//   catch (const madness::MadnessException& e)
//   {
//     print(e);
//     error("caught a MADNESS exception");
//   }
//   catch (const madness::TensorException& e)
//   {
//     print(e);
//     error("caught a Tensor exception");
//   }
//   catch (const char* s)
//   {
//     print(s);
//     error("caught a string exception");
//   }
//   catch (const std::string& s)
//   {
//     print(s);
//     error("caught a string (class) exception");
//   }
//   catch (const std::exception& e)
//   {
//     print(e.what());
//     error("caught an STL exception");
//   }
//   catch (...)
//   {
//     error("caught unhandled exception");
//   }
// 
//   if (world.rank() == 0)
//     print("entering final fence");
//   world.gop.fence();
//   if (world.rank() == 0)
//     print("done with final fence");
//   if (world.rank() == 0)
//     print("Final tensor instance count", BaseTensor::get_instance_count());
// 
//   world.gop.fence();
//   finalize();
//   return 0;
// }


 
