#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
//#include <mra/operator.h>
//#include <world/world.h>

extern "C" void elk_load_wann_unk_(int *n);
extern "C" void elk_xc_(double *dens, double *vx, double *vc, double *ex, double *ec);
extern "C" void elk_wan_rho_(double *x, double *rcutoff, double *wrho);
extern "C" void sic_wan_rho_(int *j, double *x, double *rcutoff, double *wrho);

using namespace madness;

World* world;
real_function_3d wan_rho;
real_function_3d wan_hpot;
real_function_3d wan_exc;
// good accuracy
const double thresh = 0.00001;
const int ord_k = 8;
// fast debug accuracy
//const double thresh = 0.01;
//const int ord_k = 5;

class WannierRho: public FunctionFunctorInterface<double,3> 
{
    public:
    
        double operator()(const coord_3d& x) const 
        {
            double val;
            double vrc[3];
            vrc[0] = x[0]; 
            vrc[1] = x[1]; 
            vrc[2] = x[2];
            double rcutoff = 10.0;
            int j = 1;
            elk_wan_rho_(vrc, &rcutoff, &val);
            //sic_wan_rho_(&j, vrc, &rcutoff, &val);
            return val;
        }
};

class WannierExc: public FunctionFunctorInterface<double,3> 
{
    public:
    
        double operator()(const coord_3d& x) const 
        {
            //double d = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
            //if (d > 10.0) return 0.0;
            
            double rho;// = wan_rho(x);
            double ex, ec, vx[2], vc[2];
            double rcutoff = 10.0;
            double vrc[3];
            vrc[0] = x[0]; 
            vrc[1] = x[1]; 
            vrc[2] = x[2];
            elk_wan_rho_(vrc, &rcutoff, &rho);

            elk_xc_(&rho, vx, vc, &ex, &ec);
            return (ex + ec)*rho;
        }
};


extern "C" void madness_init_(void) 
{
    initialize(0, 0);
    world = new World(MPI::COMM_WORLD);
    printf("Hello from process %i\n", world->rank());
    startup(*world, 0, 0);
}

extern "C" void madness_init_box_(void) 
{
    double L = 12.0;

    FunctionDefaults<3>::set_k(ord_k);
    FunctionDefaults<3>::set_thresh(thresh);
//    FunctionDefaults<3>::set_refine(true);
//    FunctionDefaults<3>::set_initial_level(5);
//    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
}

extern "C" void madness_gen_hpot_(int* n) 
{
    if (world->rank() == 0) printf("in madness_gen_hpot()\n");
    real_convolution_3d op = CoulombOperator(*world, thresh * 0.001, thresh * 0.001);
    wan_rho.clear();
    wan_hpot.clear();
    elk_load_wann_unk_(n);
    wan_rho = real_factory_3d(*world).functor(real_functor_3d(new WannierRho()));
    if (world->rank() == 0) 
    {
        double norm = wan_rho.trace();
        printf(" n : %i , norm : %15.8f\n", *n, norm);
    }
    wan_hpot = op(wan_rho).truncate();
    double eh = inner(wan_rho, wan_hpot);
    if (world->rank() == 0) printf("  Hartree energy : %15.8f\n", eh);
    
    wan_exc.clear();
    wan_exc = real_factory_3d(*world).functor(real_functor_3d(new WannierExc()));
    //double exc = inner(wan_rho, wan_exc);
    double exc = wan_exc.trace();
    if (world->rank() == 0) printf("       XC energy : %15.8f\n", exc);
    world->gop.fence(); 
}
extern "C" void madness_gen_xc_(int* n)
{
    if (world->rank() == 0) printf("in madness_gen_xc()\n");
    wan_rho.clear();
    elk_load_wann_unk_(n);
    wan_rho = real_factory_3d(*world).functor(real_functor_3d(new WannierRho()));
    double norm = wan_rho.trace();
    if (world->rank() == 0) printf(" n : %i , norm : %15.8f\n", *n, norm);

    wan_exc.clear();
    wan_exc = real_factory_3d(*world).functor(real_functor_3d(new WannierExc()));
    double exc = inner(wan_rho, wan_exc);
    if (world->rank() == 0) printf(" n : %i , norm : %15.8f , energy : %15.8f\n", *n, norm, exc);
    world->gop.fence(); 
}

/*extern "C" void madness_get_hpot_(double *vrc,double* pot) {
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

*/
