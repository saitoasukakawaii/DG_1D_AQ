//
// Created by chen on 2021/11/11.
//

#ifndef DG_1D_ARTERY_H
#define DG_1D_ARTERY_H

#include "Element.h"
#include "Parameters.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <cassert>
#include <filesystem>
#include <Eigen/Dense>

// define some non-dimension
#ifndef MAX_ITER
#define MAX_ITER 1000
#endif
#ifndef TOL
#define TOL 1e-20
#endif
#ifndef SMALL
#define SMALL 1e-200
#endif
//extern double   conv, rho, mu, mu_pl, nu, Lr, Lr2, Lr3, g, q, Fr2,
//        Re, p0, pmean, tmst, Period, Fcst, alpha, CO, COm,
//        Deltat;
using constants::tmstps;
using constants::Lr;
using constants::rho;
using constants::g;
using constants::N;
using constants::q;
using constants::CO_filename;
using constants::mu_pl;
using constants::Fr2;
using constants::u;

class Artery {

public:
    const int N_e;              // Element number
    vector<Element> el;               // Element

    int      ID;
    Artery   *LD, *RD;          // left and right daughter arteries.
    Artery   *LP, *RP;          // left and right parent arteries.

    // Geometry
    double   L;                 // length of the artery
    double   rt, rb;            // top and bottom radius
    double   ff1,ff2,ff3;       // use for EH/r0
    double   termresist;        // for terminal arteries
    double   pts;               // number of point per 1cm
    double   Courant;           // CFL condition.
    double   Period;
    double   Q0[tmstps+1];                // use for inlet condition
    double   pL[tmstps];                // use for terminal arteries, and conve integral
    double   y[tmstps];

    Artery(const double &length,const double &rtop,const double &rbottom, const int &id,
           Artery *LeftD, Artery *RightD, Artery *LeftP, Artery *RightP,
           const double &rmin, const double &points, const int &init,
           const double &f1,const double &f2,const double &f3,
           const double &fa1,const double &fa2,const double &fa3,
           const double &trmrst, const double &Period);
    ~Artery() = default;

    inline void Set_F() {
        for(int i=0;i<N_e;++i)
        {
            el[i].Set_F();
        }
    }

    inline void Set_S() {
        for(int i=0;i<N_e;++i)
        {
            el[i].Set_S();
        }
    }
    inline void Set_RHS() {
        for(int i=0;i<N_e;++i)
        {
            el[i].Set_RHS();
        }
    }
    inline void Update(const double &a, const double &b, const double &c, const double &dt) {
        for(int i=0;i<N_e;++i)
        {
            el[i].Update(a,b,c,dt);
        }
    }


    double get_radius(const double &x);
    void local_CFL();

    void Inlet_Flux(const double &T);
    void Inter_Flux();
    void Bifur_Flux();
    void Terminal_Flux(const int &n_step, const int &qLnb, const double &dt);
//    void ReSetK(){for (int i=0;i<N_e;++i){el[i].ReSetK();}} // dont need for the rk4a[0]=0
    double Get_Q(const double &T); // inlet flow rate only for ascending aorta

    void printQxt(std::ofstream &fd, const double &t) const;
    void printUxt(std::ofstream &fd, const double &t) const;
    void printAxt(std::ofstream &fd, const double &t) const;
    void printPxt(std::ofstream &fd, const double &t) const;
    void printUx(std::ofstream &fd) const;
    void printX(std::ofstream &fd) const;
    void printR0(std::ofstream &fd) const;
    void printPx(std::ofstream &fd) const;
    void Save_inverse_Z(std::string &filename);
    void Read_inverse_Z(std::string &filename);
    inline void Update_pL(const int &qLnb){
        pL[(qLnb % tmstps)] = el[N_e-1].P[el[N_e-1].Np-1];
    }
    inline void clean_k(){
        for(int i=0;i<N_e;++i){
            el[i].ReSetK();
        }
    }
//    void teminal_boundary();
//    void bif_boundary();
};
void BioFlux(const int &nbrves, Artery *Arteries[],
             const std::set<int>& ID_Bif, const std::set<int>& ID_Out,
             const int &n_step, const int &qLnb, const double &dt);
// solving Riemann problem for interface of two elements
void Riemann(Element &el1, Element &el2, const int &j, const int &N_e);
// solving Riemann problem for the element at the start of the artery
void RiemannStart(Element &start_el, const double &Q_star, const int &ID);
// solving Riemann problem for the element at the end of the artery
void RiemannEnd(Element &end_el, const double &Q_star, const int &ID);
void solver(const int &nbrves, Artery *Arteries[], int &n_step,
            const double &tstart, const double &tend, const double &dt,
            const std::set<int>& ID_Out, const std::set<int>& ID_Bif);
void solverRHS(const int &nbrves, Artery *Arteries[],
               const std::set<int>& ID_Bif, const std::set<int>& ID_Out,
               const int &n_step, const int &qLnb, const double &dt);
double Get_CFL(const int &nbrves, Artery *Arteries[]);

#endif //DG_1D_ARTERY_H
