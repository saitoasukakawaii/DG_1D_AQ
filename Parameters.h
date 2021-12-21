//
// Created by chen on 2021/11/15.
//

#ifndef DG_1D_PARAMETERS_H
#define DG_1D_PARAMETERS_H
#include <string>
#include <vector>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define MAX_ITER 2000
#define TOL 1e-20
using std::vector;
namespace constants{
const std::string CO_filename = "input.dat";    // Input flow file at the heart.
const std::string FileName1 = "Geometry/topology55.txt";
const std::string FileName2 = "Geometry/geometry55.txt";

const int N = 3;                                    // N order, N+1 point in a element
const Jacobi1D jac(N);                              // N order element
                                            // 65536 32768
const int tmstps = 32768;                  // The number of timesteps per period.
const int plts   = 2048;                   // Number of plots per period.
const int tmp    = 1024;
const int Per_step = tmstps/tmp;

const double point  = 4.;                   // 1cm divide to 8 point
const double conv   = 1332.20;                  // Conversion from mmHg to SI-units.
const double rho    = 1.055;                    // Density of blood [g/cm^3].
const double mu     = 0.049;                    // Viscosity of blood [g/cm/s].
const double mu_pl  = mu;                       // Viscosity of blood [g/cm/s].
const double nu     = mu/rho;

//         Tper   = 1.00,                       // The period of one heart beat [s].
                                                // not set as const for changing sake.


const double ff1    = 1.99925e+07;              // artery
const double ff2    = -22.5267;
const double ff3    =  865251;
const double fa1    = 1.89e+07;                 // small artery
const double fa2    = -18.34;
const double fa3    = 3.53e+06;

const double alpha  = 1.1;                      // velocity profile coefficients.

const double Lr     = 1.0;                      // Characteristic radius of the
                                                // vessels in the tree [cm].
const double Lr2    = Lr*Lr;                    // The squared radius [cm2].
const double Lr3    = Lr*Lr2;                   // The radius to the third power [cm^3].
const double g      = 981.0;                    // The gravitational force [cm/s^2].
const double u      = 10.0;                     // The characteristic velocity [cm/s].
const double q      = 10.0*Lr2;                 // The characteristic flow [cm^3/s].
const double Fr2    = u*u/g/Lr;                 // The squared Froudes number.
const double Re     = q*rho/mu/Lr;          // Reynolds number.
// Period = Tper*q/Lr3,           // The dimension-less period.
const double Fcst   = -M_PI*2*alpha/(alpha-1)/Re;
                                                // coeff of friction. correspond to alpha
// k      = Period/tmstps,        // Length of a timestep.
// Deltat = Period/plts,          // Interval between each point plottet.
// p0     = 55.0/rho/g/Lr*conv,   // Ensures a certain diastolic pressure.
const double p0     = 55.*conv;   // Ensures a certain diastolic pressure.
}


namespace RK4{
    const unsigned int N_t = 5;
    double const rk4a[N_t] = { 0.0,
                               -567301805773.0/1357537059087.0,
                               -2404267990393.0/2016746695238.0,
                               -3550918686646.0/2091501179385.0,
                               -1275806237668.0/842570457699.0};

    double const rk4b[N_t] = {1432997174477.0/9575080441755.0,
                              5161836677717.0/13612068292357.0,
                              1720146321549.0/2090206949498.0,
                              3134564353537.0/4481467310338.0,
                              2277821191437.0/14882151754819.0};

    double const rk4c[N_t] = {0.0,
                              1432997174477.0/9575080441755.0,
                              2526269341429.0/6820363962896.0,
                              2006345519317.0/3224310063776.0,
                              2802321613138.0/2924317926251.0};
}

#endif //DG_1D_PARAMETERS_H
