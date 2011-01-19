#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <mra/operator.h>
#include <world/world.h>

extern "C" void elk_load_wann_unk_(int* n);
extern "C" void elk_wan_rho_(const int* n,double* r_cutoff,double* vrc,double* wrho);

using namespace madness;

World* world;
real_function_3d wan_hpot;
const double thresh = 0.02;
const int ord_k = 5;

class WannierRho: public FunctionFunctorInterface<double,3> {
public:
  int n_;

  WannierRho(int n) : n_(n) {}

  double operator()(const coord_3d& x) const {
    double val;
    double vrc[3];
    vrc[0]=x[0]; vrc[1]=x[1]; vrc[2]=x[2];
    double r_cutoff=10.0;
    
    elk_wan_rho_(&n_,&r_cutoff,&vrc[0],&val);
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

extern "C" void madness_gen_hpot_(int* nwantot) {
  real_convolution_3d op = CoulombOperator(*world, thresh*0.1, thresh*0.1);
  wan_hpot.~real_function_3d();
//   for (int n=1;n<=*nwantot;n++) {  
//     elk_load_wann_unk_(&n);
//     real_function_3d rho = real_factory_3d(*world).functor(real_functor_3d(new WannierRho(n)));
//     real_function_3d hpot = op(rho).truncate();
//     wan_hpot.push_back(hpot);
//     double norm = rho.trace();
//     printf(" n : %i , norm : %15.8f\n",n,norm);
//     world->gop.fence(); 
//   }
}

extern "C" void madness_get_hpot_(double *x,double* pot) {
  //*pot=wan[0](0.0);
}
