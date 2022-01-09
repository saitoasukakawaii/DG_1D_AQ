//
// Created by chen on 2021/11/9.
//

#ifndef DG_1D_ELEMENT_H
#define DG_1D_ELEMENT_H
#include "Jacobi1D.h"
#include "Parameters.h"
#include <memory>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
using constants::Lr;
using constants::rho;
using constants::g;
using constants::Fr2;
using constants::mu;
using constants::alpha;
using constants::F_c;

class Element {
public:
    int      N;                     // number of polynomial
    int      Np;                    // number of Gauss point
    double   x_l, x_r;              // left right vortex
    double   L;                     // 长度
    double   Min_h;
    double   Courant;               // local CFL
    double   rx;
    double   a;                     // coef for dr0dx, log(rb/rt)

    vector<double>  VX,                    // x-cord of quadrature point r
            A, A0,                // area
            r0, dr0dx,                   // 0-pressure radius
            EH_r0, dfdr0,
            Q,
                                    // _A for continuous, _U for momentum
            F_A, F_U,             // Flux term
    // S_A == 0, dont need
            S_U,                   // source term
            RHS_A, RHS_U,         // all value add to right
            k_A, k_U,             // used for Runge-Kutta(LSERK), as tmp
                                    // dont need reset for the rk4a[0]=0
            P, phi,                     // pressure
            c,c0;                 // wave speed

//    double   FA1, FA2,              // upwind flux 1 for left edge, 2 for right edge
//             FU1, FU2;              // A for continuous U for momentum
    double   dFA1, dFA2,                // use for numerical flux
             dFU1, dFU2;                //
    double   W1L, W1R,              // riemann L for left backward, R for right forward.
             W2L, W2R;              // 1 for left edge, 2 for right edge

    Element() = default;
    Element(const double &_x_l, const double &_x_r, const int &_N);
//    Element(Element &&el);
    ~Element() = default;
    void MeshGen1D();        // meshing, get all physic x
    void Set_A0();           // set A0 and A
    void Set_EH_r0(const double &ff1, const double &ff2, const double &ff3);
    void Set_dfdr0(const double &ff1, const double &ff2);
    void Set_c0();           // set c0 and c
    void Set_c();            // set new c, need new A
    void Set_P();
    void Set_phi();
//    inline void Set_Q() {
//        for (int i = 0; i < Np; i++) {
//            Q[i] = A[i] * U[i];
//        }
//    }
    void local_CFL();
    double Get_c(const int &i, const double &A_);
    double Get_P(const int &i, const double &A_);
    double Get_phi(const int &i, const double &A_);
    void Set_F();
    void Set_S();
    double Get_FU1(const double &Q, const double &A) {return Q*Q/A +Get_phi(0, A);};
    double Get_FU2(const double &Q, const double &A) {return Q*Q/A +Get_phi(Np-1, A);};
    void Set_RHS();
    double Get_dpdx(const int &i);
//    double Hp (const int &i);
//    void poschar (const double &theta,
//                             double &qR, double &aR,
//                             double &cR, double &HpR);
//    double Hn (const int &i);
//    void negchar (const double &theta,
//                  double &qS, double &aS,
//                  double &cS, double &HnS);
    inline void Set_W1() {
        W1L = Q[0]/A[0]-4*c[0];
        W1R = Q[0]/A[0]+4*c[0];
    }
    inline void Set_W2() {
        W2L = Q[Np-1]/A[Np-1]-4*c[Np-1];
        W2R = Q[Np-1]/A[Np-1]+4*c[Np-1];
    }
    void Update(const double &a, const double &b, const double &dt);
    inline void ReSetK(){
        for (int i = 0; i < Np; ++i) {
            k_A[i] = 0.;
            k_U[i] = 0.;
        }
    }
};


#endif //DG_1D_ELEMENT_H
