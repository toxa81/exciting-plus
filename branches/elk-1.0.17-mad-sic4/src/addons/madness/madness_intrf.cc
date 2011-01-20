#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <mra/operator.h>
#include <world/world.h>

extern "C" void elk_load_wann_unk_(int* n);
extern "C" void elk_wan_rho_(double* r_cutoff,double* vrc,double* wrho);

using namespace madness;

World* world;
real_function_3d wan_rho;
real_function_3d wan_hpot;
// good accuracy
const double thresh = 0.0005;
const int ord_k = 8;
// fast debug accuracy
//const double thresh = 0.01;
//const int ord_k = 5;

class WannierRho: public FunctionFunctorInterface<double,3> {
public:

  double operator()(const coord_3d& x) const {
    double val;
    double vrc[3];
    vrc[0]=x[0]; vrc[1]=x[1]; vrc[2]=x[2];
    double r_cutoff=10.0;
    
    elk_wan_rho_(&r_cutoff,&vrc[0],&val);
    return val;
  }
};

extern "C" void madness_init_(void) {
  initialize(0, 0);
  world = new World(MPI::COMM_WORLD);
  printf("Hello from process %i\n", world->rank());
  startup(*world, 0, 0);
}

extern "C" void madness_init_box_(void) {
  double L = 30.0;

  FunctionDefaults<3>::set_k(ord_k);
  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(5);
//  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
}

extern "C" void madness_gen_hpot_(int* n) {
  real_convolution_3d op = CoulombOperator(*world, thresh*0.001, thresh*0.001);
  wan_rho.clear();
  wan_hpot.clear();
  //elk_load_wann_unk_(&n);
  wan_rho = real_factory_3d(*world).functor(real_functor_3d(new WannierRho()));
  wan_hpot = op(wan_rho).truncate();
  double norm = wan_rho.trace();
  if (world->rank() == 0) printf(" n : %i , norm : %15.8f\n",*n,norm);
  world->gop.fence(); 
}

extern "C" void madness_get_hpot_(double *vrc,double* pot) {
  coord_3d x;
  x[0]=vrc[0]; x[1]=vrc[1]; x[2]=vrc[2];
  if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > 10.0) {
    *pot=0.0;
  } else { 
    *pot=wan_hpot(x);
  }
}

extern "C" void madness_get_rho_(double *vrc,double* rho) {
  coord_3d x;
  x[0]=vrc[0]; x[1]=vrc[1]; x[2]=vrc[2];
  if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) > 10.0) {
    *rho=0.0;
  } else { 
    *rho=wan_rho(x);
  }
}
