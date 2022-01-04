//
// Created by chen on 2021/11/9.
//

#include "Element.h"

Element::Element(const double &_x_l, const double &_x_r, const int &_N) :
// const init here
        x_l(_x_l), x_r(_x_r), L(_x_r - _x_l), N(_N), Np(_N + 1), rx(2 / (L)) {
    // init
    if (_x_r < _x_l) throw std::runtime_error("Wrong mesh for negtive mesh, xl should be less than lr!!!\n");
    VX = vector<double>(Np,0.);
    A0 = vector<double>(Np,0.);
    r0 = vector<double>(Np,0.);
    c0 = vector<double>(Np,0.);
    A = vector<double>(Np,0.);
//    Q = vector<double>(Np,10.);
    Q = vector<double>(Np,1.);
    c = vector<double>(Np,0.);
    P = vector<double>(Np,0.);
    phi = vector<double>(Np,0.);
    F_A = vector<double>(Np,0.);
    F_U = vector<double>(Np,0.);
    S_U = vector<double>(Np,0.);
    k_A = vector<double>(Np,0.);
    k_U = vector<double>(Np,0.);
    RHS_A = vector<double>(Np,0.);
    RHS_U = vector<double>(Np,0.);
    RHS_A_old = vector<double>(Np,0.);
    RHS_U_old = vector<double>(Np,0.);
    RHS_A_oldold = vector<double>(Np,0.);
    RHS_U_oldold = vector<double>(Np,0.);
    EH_r0 = vector<double>(Np,0.);
    dr0dx = vector<double>(Np,0.);
    dfdr0 = vector<double>(Np,0.);
    MeshGen1D();
}

//Element::Element(Element &&el)
//{
//
//}
//Element::~Element() {
//    delete[] VX;
//    delete[] A;
//    delete[] A0;
//    delete[] r0;
//    delete[] EH_r0;
//    delete[] U;
////    delete[] Q;
//    delete[] F_A;
//    delete[] F_U;
//    delete[] S_U;
//    delete[] RHS_A;
//    delete[] RHS_U;
//    delete[] k_A;
//    delete[] k_U;
//    delete[] P;
//    delete[] c;
//    delete[] c0;
//}

// Gauss point
void Element::MeshGen1D() {
    using constants::jac;
    for (int i = 0; i < Np; i++) { VX[i] = 0.5 * x_r * (1 + jac.r(i)) + 0.5 * x_l * (1 - jac.r(i)); }

    // min point distance for CFL
    Min_h = VX[1] - VX[0];
    assert(Min_h > 0);
    for (int i = 2; i < Np; i++) {
        double temp = VX[i] - VX[i - 1];
        assert(temp > 0);
        Min_h = std::min(Min_h, temp);
    }
    return;
}

void Element::Set_A0() {
    for (int i = 0; i < Np; ++i) {
        A0[i] = M_PI * r0[i] * r0[i];
        A[i] = A0[i];
    }
}

// olufsen
void Element::Set_EH_r0(const double &ff1, const double &ff2, const double &ff3) {
    for (int i = 0; i < Np; ++i) {
        EH_r0[i] = (ff1 * exp(ff2 * r0[i]) + ff3) * 4 / 3;
    }
}

void Element::Set_dfdr0(const double &ff1, const double &ff2) {
    for (int i = 0; i < Np; ++i) {
        dfdr0[i] = dr0dx[i]*(ff1 * ff2 * exp(ff2 * r0[i])) * 4 / 3;
    }
}

void Element::Set_c0() {
    for (int i = 0; i < Np; ++i) {
        c0[i] = sqrt(0.5 * EH_r0[i] / rho);
        c[i] = c0[i];
    }
}

void Element::Set_c() {
    for (int i = 0; i < Np; ++i) {
        c[i] = sqrt( 0.5*EH_r0[i]*sqrt(A0[i]/A[i]) / rho );         // olufsen
        //
//        c[i] = sqrt(0.5 * EH_r0[i] * sqrt(A[i] / A0[i]) / rho);
    }
}



void Element::Set_P() {
    for (int i = 0; i < Np; i++) {
        P[i] = EH_r0[i]*(1-sqrt(A0[i]/A[i]));             // olufsen
//        P[i] = EH_r0[i] * (sqrt(A[i] / A0[i]) - 1);
    }
}

void Element::Set_phi() {
    for (int i = 0; i < Np; i++) {
        phi[i] = EH_r0[i]*( sqrt(A0[i]*A[i])-A[0] ) / rho;             // olufsen
//       phi[i] = EH_r0[i] * ( pow(A[i],1.5)-pow(A0[i],1.5) ) /rho/3/sqrt(A0[i]);
    }
}

double Element::Get_c(const int &i, const double &A_) {
    assert(i >= 0 && i < Np);
    return sqrt( 0.5*EH_r0[i]*sqrt(A0[i]/A_)/rho );       // olufsen
//    return sqrt(0.5 * EH_r0[i] * sqrt(A_ / A0[i]) / rho);
}


double Element::Get_P(const int &i, const double &A_) {
    assert(i >= 0 && i < Np);
    return EH_r0[i]*(1-sqrt(A0[i]/A_));               // olufsen
//    return EH_r0[i] * (sqrt(A_ / A0[i]) - 1);
}

double Element::Get_phi(const int &i, const double &A_) {
    assert(i >= 0 && i < Np);
    return EH_r0[i]*( sqrt(A0[i]*A_)-A[0] ) / rho;             // olufsen
//    return EH_r0[i] * ( pow(A[i],1.5)-pow(A0[i],1.5) ) /rho/3/sqrt(A0[i]);
}

// need change, length is wrong
void Element::local_CFL() {
    Courant = 64000000.0;
    for (int i = 0; i < Np; i++) {
        double temp = std::min(Min_h / fabs(Q[i]/A[i] - c[i]),
                               Min_h / fabs(Q[i]/A[i] + c[i]));
        Courant = std::min(temp, Courant);
    }
}

void Element::Set_F() {

    Set_c();
    Set_P();
    Set_phi();
    for (int i = 0; i < Np; i++) {
        F_A[i] = Q[i];
    }

    for (int i = 0; i < Np; i++) {
        F_U[i] = Q[i]*Q[i]/A[i] +phi[i];
    }
}

void Element::Set_S() {
    // S_A don't need
    // attention to the minus sign
    // double F_c = -2*M_PI*mu*alpha/rho/(alpha-1);
    double F_c = -22 * M_PI * mu;

    for (int i = 0; i < Np; i++) {
//        S_U[i] = ( F_c * Q[i] / A[i] + dfdr0[i]*dr0dx[i]*( A[i]-2/3*pow(A[i],1.5)/sqrt(M_PI)/r0[i] )
//                + dr0dx[i]*EH_r0[i]*2/3*pow(A[i],1.5)/sqrt(M_PI)/pow(r0[i],2.) ) / rho;
        // olufsen
        S_U[i] = ( F_c * Q[i] / A[i]
                + dfdr0[i]*dr0dx[i]*( sqrt(A[i]*A[0])*2-A[i]-A0[i] )
                + EH_r0[i]*dr0dx[i]*2*( sqrt(A[i]*M_PI)-M_PI*r0[i] ) ) / rho;
    }
}

void Element::Set_F1() {
    FA1 = Q1;
    FU1 = Q1*Q1/A1 + Get_phi(0, A1);
}

void Element::Set_F2() {
    FA2 = Q2;
    FU2 = Q2*Q2/A2 + Get_phi(Np-1, A2);
}

void Element::Set_RHS() {
    // need carefull cope
    auto Lp1 = constants::jac.invM(all, 0);
    auto Lp2 = constants::jac.invM(all, Np - 1);
    double dFA1 = F_A[0] - FA1;
    double dFA2 = F_A[Np - 1] - FA2;
    double dFU1 = F_U[0] - FU1;
    double dFU2 = F_U[Np - 1] - FU2;
    double dFAdx[Np];
    double dFUdx[Np];
    for (int i = 0; i < Np; ++i) {
        dFAdx[i] = 0.;
        dFUdx[i] = 0.;
        for (int j = 0; j < Np; ++j) {
            dFAdx[i] += constants::jac.Dr(i,j) * F_A[j];
            dFUdx[i] += constants::jac.Dr(i,j) * F_U[j];
        }
    }
    for (int i = 0; i < Np; i++) {
        RHS_A_oldold[i] = RHS_A_old[i];
        RHS_U_oldold[i] = RHS_U_old[i];
        RHS_A_old[i] = RHS_A[i];
        RHS_U_old[i] = RHS_U[i];
    }
    for (int i = 0; i < Np; i++) {
        RHS_A[i] = -rx * dFAdx[i] + rx * dFA2 * Lp2(i) - rx * dFA1 * Lp1(i);
        RHS_U[i] = -rx * dFUdx[i] + rx * dFU2 * Lp2(i) - rx * dFU1 * Lp1(i) + S_U[i];
    }
    

}

void Element::Update(const double &a, const double &b, const double &c, const double &dt) {
    for (int i = 0; i < Np; i++) {
        A[i] = A[i] + dt*(a*RHS_A[i]+b*RHS_A_old[i]+c*RHS_A_oldold[i]);
        Q[i] = Q[i] + dt*(a*RHS_U[i]+b*RHS_U_old[i]+c*RHS_U_oldold[i]);
    }
}


// double Tube :: Hp (int i, double Q, double A)
// {
//   return (F(Q,A) - A*dPdx1(i,A)/Fr2)/(-Q/A + c(i,A));
// }

// double Tube :: Hn (int i, double Q, double A)
// {
//     return (F(Q,A) - A*dPdx1(i,A)/Fr2)/(-Q/A - c(i,A));
// }

// void Element::poschar (const double &dt, double &qR, double &aR, double &cR, double &HpR)
// {
//   double ctm1  = c[Np-1];
//   double Hptm1 = Hp (N, Qold[N], Aold[N]);
//   double uR    = Qold[N] / Aold[N];
//   double ch    = (uR + ctm1) * theta;
//
//   if (uR + ctm1 < 0)
//   {
//     throw("uR + ctm1 < 0, flow is not subsonic!\n");
//   }
//
//   qR  = Qold[N] - (Qold[N] - Qold[N-1])*ch;
//   aR  = Aold[N] - (Aold[N] - Aold[N-1])*ch;
//   cR  = ctm1    - (ctm1  - c (N-1,Aold[N-1]))*ch;
//   HpR = Hptm1   - (Hptm1 - Hp(N-1,Qold[N-1],Aold[N-1]))*ch;
// }
//
//void Element::negchar (const double theta, double &qS, double &aS, double &cS, double &HnS)
// {
//     double ctm1  = c[0];
//     double Hntm1 = Hn(0, Qold[0], Aold[0]);
//     double uS    = Qold[0]/Aold[0];
//     double ch    = (uS - ctm1) * theta;
//
//     if ( ctm1 - uS < 0)
//     {
//         throw("ctm1 - uS < 0, CFL condition violated\n");
//     }
//
//     qS  = Qold[0] + (Qold[0] - Qold[1])*ch;
//     aS  = Aold[0] + (Aold[0] - Aold[1])*ch;
//     cS  = ctm1    + (ctm1  - c (1,Aold[1]))*ch;
//     HnS = Hntm1   + (Hntm1 - Hn(1,Qold[1],Aold[1]))*ch;
// }



