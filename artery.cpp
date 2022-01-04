//
// Created by chen on 2021/11/11.
//

#include "artery.h"
// oulet need complete, and check energy balance and supersonic?

extern "C" void impedance_driver_(int *tmstps, double *Period,
                                  double *ff1, double *ff2, double *ff3,
                                  double *rho, double *mu,
                                  double *r_root, double *rmin,
                                  double *y, double *Lr, double *Fr2, double *q, double *g, double *trmrst);


Artery::Artery(const double &length,const double &rtop,const double &rbottom, const int &id,
                     Artery *LeftD, Artery *RightD, Artery *LeftP, Artery *RightP,
                     const double &rmin, const double &points, const int &init,
                     const double &f1,const double &f2,const double &f3,
                     const double &fa1,const double &fa2,const double &fa3,
                     const double &trmrst, const double &Period):
       L(length), rt(rtop), rb(rbottom), ID(id), LD(LeftD), RD(RightD), LP(LeftP), RP(RightP),
       pts(points), ff1(f1), ff2(f2), ff3(f3), termresist(trmrst), Period(Period), N_e(int(points*length))
{
    double h = 1.0/points;
    double a = log(rb/rt)/L;
    // get the element,
//    el = new Element[N_e];
    el.resize(N_e);
    for(int i=0; i<N_e; ++i)
    {
        // left and right x-axis
        double x_l = i*h;
        double x_r = (i+1)*h;
        el[i] = Element(x_l, x_r, N);
        // calculate the Np, the number of quature points.
        for (int j=0;j<el[i].Np;++j)
        {
            el[i].r0[j] = get_radius(el[i].VX[j]);
            el[i].dr0dx[j] = a*el[i].r0[j];
        }
        el[i].Set_A0();
        el[i].Set_EH_r0(ff1,ff2,ff3);
        el[i].Set_dfdr0(ff1,ff2);
        el[i].Set_c0();
        el[i].a = a;
//        el[i].Set_P();
    }
    local_CFL();
//    std::stringstream ss;
//    std::ofstream ff;
//    ss << "Radius_" << std::setfill('0') << std::setw(2) << ID+1 << ".dat";
//    ff.open(ss.str());
//    for(int i=0; i<N_e; ++i) {
//        for (int j = 0; j < el[i].Np; ++j) {
//            ff << el[i].r0[j] << "\n";
//        }
//    }
//    ff.close();
    if (init == 1)
    {
        std::ifstream fi(CO_filename);
        if(!fi.is_open()) throw std::runtime_error(std::string{"Could not open file: "}+CO_filename);
//        Q0 = VectorXd::Constant(tmstps+1, 0.);
        for (int i=0; i<=tmstps; i++)
        {
            fi >> Q0[i];
            // Q0[i] = Q0[i];    // If the in data have dimensions they should be made
                                // non-dimensional.
//            u0[i] = Q0[i]/A0[i];
        }
        fi.close();
    }

    if (LD == nullptr && RD == nullptr)
    {

        int ATS = constants::tmstps;
        double Per = Period*u/Lr;
        double ffa1 = fa1;
        double ffa2 = fa2;
        double ffa3 = fa3;
        double N_rho = constants::rho;
        double N_mu = mu_pl;
        double r_min = rmin;
        double N_Lr = constants::Lr;
        double N_Fr2 = constants::Fr2;
        double N_q = constants::q;
        double N_g = constants::g;
        double N_ter = termresist;

        std::string filename = "Y_";
        filename += std::to_string(int(1000*rb))+"_"+std::to_string(int(r_min*1000))+"_"+std::to_string(tmstps);

        if (!std::filesystem::exists(filename)){
            std::cout << "Calling f90 subroutines." << std::endl;
            impedance_driver_(&ATS,&Per,&ffa1,&ffa2,&ffa3,&N_rho,&N_mu,&rb,&r_min,y,&N_Lr,&N_Fr2,&N_q,&N_g,&N_ter);
            Save_inverse_Z(filename);
        }else{
            std::cout << "exist inverse resistance file: " << filename << std::endl;
            Read_inverse_Z(filename);
        }


        std::cout << "Finished with f90 subroutines.\n" << std::endl;
        // Initialize the array pL used when determining the convolution
        // in the right boundary condition (see the subroutine bound_right).
//        pL = VectorXd::Constant(tmstps, 0.);
        for (int j=0; j<tmstps; j++)
        {
            pL[j] = 0.0;
        };
    }
}


//Artery::~Artery() { delete[] el;}

double Artery::get_radius(const double &x){
    assert(x>=0 && x<=L);
    // throw std::range_error("x is out of the range!\n");
    double r = rt*pow(rb/rt, x/L);
    return r;
}

void Artery::local_CFL() {
    el[0].local_CFL();
    Courant = el[0].Courant;
    for (int i=1; i<N_e; i++)
    {
        el[i].local_CFL();
        Courant = std::min(el[i].Courant, Courant);
    }
}
// set the upwind flux
// energy balanceann: iteration failed to conv
void Artery::Inter_Flux() {
    // dont deal with interface
    for (int j=1;j<N_e;++j) {
        Riemann(el[j - 1], el[j], j, N_e);
        el[j - 1].Set_F2();
        el[j].Set_F1();
        double dFA = (el[j - 1].FA2 - el[j].FA1);
        double dFU = (el[j - 1].FU2 - el[j].FU1);
//        if ( dFA*dFA > TOL*100000 ) {
//            std::stringstream tmp;
//            tmp << "Error in BioFlux: Q does not balance at the interelement boundary.\n"
//                << "The number of artery is: " << (ID + 1) << ".\n"
//                << "The order of element is: " << (j - 1) << " and " << j << ".\n"
//                << "The number of element is: " << N_e << ".\n";
//            tmp << "The Flux_star at left element is: " << el[j - 1].FA2 << " and " << el[j - 1].FU2 << ",\n"
//                << "The Flux_star at right element is: " << el[j].FA1 << " and " << el[j].FU1 << ".\n"
//                << "the differ of Q and energy is: " << dFA << " and " << dFU << ".\n";
//            throw std::runtime_error(tmp.str());
//        }
//        if ( dFU*dFU > TOL*100000 ) {
//            std::stringstream tmp;
//            tmp << "Error in BioFlux: energy does not balance at the interelement boundary.\n"
//                << "The number of artery is: " << (ID + 1) << ".\n"
//                << "The order of element is: " << (j - 1) << " and " << j << ".\n"
//                << "The number of element is: " << N_e << ".\n";
//            tmp << "The Flux_star at left element is: " << el[j - 1].FA2 << " and " << el[j - 1].FU2 << ",\n"
//                << "The Flux_star at right element is: " << el[j].FA1 << " and " << el[j].FU1 << ".\n"
//                << "the differ of Q and energy is: " << dFA << " and " << dFU << ".\n";
//            throw std::runtime_error(tmp.str());
//        }
    }
}
// subsonic check before flux calculate then check energy balance
void Artery::Bifur_Flux()
{
    // f: function g: delta_x,
    // x: x_new = x_old+delta_x
    // inv_j: inv of f'
    MatrixXd J(6,6);
    VectorXd f(6), dx(6), x(6);
//    double f[6],dx[6],x[6];
//    double J[6][6];
//    double inv_J[6][6];

    double cl, cr1, cr2, k, k1, k2, k3;
    int proceed = 1, iter = 0;

    Element &elP = el[N_e-1];
    Element &eld1 = LD->el[0];
    Element &eld2 = RD->el[0];
    // left element, use right, index: 2
    elP.Set_W2();
    // right element, use left, index: 1
    eld1.Set_W1();
    eld2.Set_W1();
    int N_ = elP.Np-1;
    x[0] = elP.Q[N_];                          // U^P_L
    x[1] = eld1.Q[0];                   // U^d1_R
    x[2] = eld2.Q[0];                   // U^d2_R
    x[3] = elP.A[N_];                          // A^P_L
    x[4] = eld1.A[0];                   // A^d1_R
    x[5] = eld2.A[0];                   // A^d2_R


    if(fabs(x[0]/x[3]) > fabs(elP.Get_c(N_, x[3]))){
        std::stringstream tmp;
        tmp << "Error in Bifurcation: flow is not subsonic. the number of artery is " << (ID+1) << ".\n";
        throw std::runtime_error(tmp.str());
    }
    if(fabs(x[1]/x[4]) > fabs(eld1.Get_c(0, x[4]))){
        std::stringstream tmp;
        tmp << "Error in Bifurcation: flow is not subsonic. the number of left daughter artery is " << (LD->ID+1) << ".\n";
        throw std::runtime_error(tmp.str());
    }
    if(fabs(x[2]/x[5]) > fabs(eld1.Get_c(0, x[5]))){
        std::stringstream tmp;
        tmp << "Error in Bifurcation: flow is not subsonic. the number of right daughter artery is " << (RD->ID+1) << ".\n";
        throw std::runtime_error(tmp.str());
    }

    while((proceed)&&(iter++ < MAX_ITER)){
        // Calculate constraint vector and the wave speed at each vessel.

        cl   = elP.Get_c(N_, x[3]);
        cr1  = eld1.Get_c(0, x[4]);
        cr2  = eld2.Get_c(0, x[5]);
        f[0] = x[0]/x[3] + 4*cl - elP.W2R;
        f[1] = x[1]/x[4] - 4*cr1 - eld1.W1L;
        f[2] = x[2]/x[5] - 4*cr2 - eld2.W1L;
        f[3] = x[0]-x[1]-x[2];
        double P_ALL = 0.5*pow(x[0]/x[3],2.0)*rho+elP.Get_P(N_,x[3]);
        f[4] = P_ALL-0.5*pow(x[1]/x[4],2.0)*rho-eld1.Get_P(N_,x[4]);
        f[5] = P_ALL-0.5*pow(x[2]/x[5],2.0)*rho-eld2.Get_P(N_,x[5]);

        // x(0): U^P_L, x(1): U^d1_R, x(2): U^d2_R;
        // x(3): A^P_L, x(4): A^d1_R, x(5): A^d2_R.
        // Jacobian matrix dfd0u
        /* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
        // x(0~2): Q x(3~5): A
        J(0,0) = 1/x[3];
        J(0,1) = 0;
        J(0,2) = 0;
//        J(0,3) = -x[0]/pow(x[3],2.)-cl/x[3];    // olufsen
        J(0,3) = -x[0]/pow(x[3],2.)+cl/x[3];
        J(0,4) = 0;
        J(0,5) = 0;

        J(1,0) = 0;
        J(1,1) = 1/x[4];
        J(1,2) = 0;
        J(1,3) = 0;
//        J(1,4) = -x[1]/pow(x[4],2.)+cr1/x[4];   // olufsen
        J(1,4) = -x[1]/pow(x[4],2.)-cr1/x[4];
        J(1,5) = 0;

        J(2,0) = 0;
        J(2,1) = 0;
        J(2,2) = 1/x[5];
        J(2,3) = 0;
        J(2,4) = 0;
//        J(2,5) = -x[2]/pow(x[5],2.)+cr2/x[5];   // olufsen
        J(2,5) = -x[2]/pow(x[5],2.)-cr2/x[5];

        J(3,0) =  1;
        J(3,1) = -1;
        J(3,2) = -1;
        J(3,3) = 0;
        J(3,4) = 0;
        J(3,5) = 0;

        J(4,0) =  rho*x[0]/pow(x[3],2.);
        J(4,1) = -rho*x[1]/pow(x[4],2.);
        J(4,2) = 0;
        J(4,3) = -rho*(pow(x[0],2.)/pow(x[3],3.)-cl*cl/x[3]);
        J(4,4) =  rho*(pow(x[1],2.)/pow(x[4],3.)-cr1*cr1/x[4]);
        J(4,5) = 0;

        J(5,0) =  rho*x[0]/pow(x[3],2.);
        J(5,1) = 0;
        J(5,2) = -rho*x[2]/pow(x[5],2.);
        J(5,3) = -rho*(pow(x[0],2.)/pow(x[3],3.)-cl*cl/x[3]);
        J(5,4) = 0;
        J(5,5) =  rho*(pow(x[2],2.)/pow(x[5],3.)-cr2*cr2/x[5]);

        dx = J.partialPivLu().solve(f);
        // Update solution: x_new = x_old - dx
        x = x-dx;
//        for (int j = 0; j < 6; ++j) {
//            x[j] -= dx[j];
//        }

        // Check if the error of the solution is smaller than TOL.
        if(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3]+dx[4]*dx[4]+dx[5]*dx[5] < TOL)      proceed = 0;
//        if( (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]+f[4]*f[4]+f[5]*f[5]) < TOL )      proceed = 0;
    }

    if(iter >= MAX_ITER){
        std::stringstream tmp;
        tmp << "Error in Bifurcation Riemann: iteration failed to converge. the number of the artery is " << (ID+1) << ".\n";
        throw std::runtime_error(tmp.str());
    }


//    double energy_in = 0.0;
//    double energy_out = 0.0;
//    if (x[0] > 0.0) energy_in = x[0]*x[3]*(0.5*rho*x[0]*x[0] + elP.Get_P(N, x[3]));
//    else energy_out = fabs(x[0])*x[3]*(0.5*rho*x[0]*x[0] + elP.Get_P(N, x[3]));
//
//    if (x[1] > 0.0) energy_out = energy_out + x[1]*x[4]*(0.5*rho*x[1]*x[1] + eld1.Get_P(0,x[4]));
//    else energy_in = energy_in + fabs(x[1])*x[4]*(0.5*rho*x[1]*x[1] + eld1.Get_P(0,x[4]));
//
//    if (x[2] > 0.0) energy_out = energy_out + x[2]*x[5]*(0.5*rho*x[2]*x[2] + eld2.Get_P(0,x[5]));
//    else energy_in = energy_in + fabs(x[2])*x[5]*(0.5*rho*x[2]*x[2] + eld2.Get_P(0,x[5]));
//    double denergy = energy_out - energy_in;
//    if(pow((denergy),2.0) > 1000000*TOL){
//
//        std::stringstream tmp;
//        tmp << "The L2 norm is:" << denergy << "\n";
//        tmp << "Error in bifurcation condition: energy per unit time OUT > energy per unit time IN \nat bifurcacion with parent vessel " << (ID+1) << ".\n";
//        tmp << "daughter vessel is: " << (this->LD->ID+1) << " and " << (this->RD->ID+1) << "\n";
//        throw std::runtime_error(tmp.str());
//    }


    eld1.W1R = x[1]/x[4]+4*eld1.Get_c(0,x[4]);
    RiemannStart(eld1,x[1],LD->ID);

    eld2.W1R = x[2]/x[5]+4*eld2.Get_c(0,x[5]);
    RiemannStart(eld2,x[2],RD->ID);

    elP.W2L = x[0]/x[3]-4*elP.Get_c(N_,x[3]);
    RiemannEnd(elP,x[0],ID);
}


void Artery::Inlet_Flux(const double &T) {
    double f, df, x, c, dx;
    int proceed = 1, iter = 0;
    double Q = Get_Q(T);
    double tol = 1e-10;
    Element &start_el = el[0];
    // set left edge riemann inviant
    start_el.Set_W1();
    x = start_el.A[0];
    // init U, A use el.A, el.U => A_start, U_start
    if( ( fabs(start_el.Q[0]/x) > start_el.c[0] ) ){
        throw std::runtime_error("Error in BioFlux: Flow is not subsonic at inlet. \n");
    }
    while((proceed)&&(iter++ < MAX_ITER)){
        // Calculate constraint vector and the wave speed at each vessel.
        c  = start_el.Get_c(0, x);

        f  = Q/x-4*c-start_el.W1L;
//        df = c/x-Q/pow(x,2.);   // olufsen
        df = -c/x-Q/pow(x,2.);
        if( fabs(df) < SMALL ){
            std::stringstream tmp;
            tmp << "Error in Qinflow: zero derivative in the Newton's iteration. Iteration step is: "<< iter << ".\n";
            throw std::runtime_error(tmp.str());
        }
        /* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
        dx = f/df;

        // Update solution: x_new = x_old - dx
        x = x-dx;

        // Check if the error of the solution is smaller than TOL.
        if(fabs(dx) < tol/100) proceed = 0;
//        if(fabs(f) < tol) proceed = 0;
    }

    if(iter >= MAX_ITER){
        std::stringstream tmp;
        tmp << "Error in inlet Riemann: iteration failed to converge. \n"<< iter << ".\n";
        throw std::runtime_error(tmp.str());
//        throw std::runtime_error("Error in inlet Riemann: iteration failed to converge. \n");
    }

    double U_star = Q/x;
    // visual edge
    double local_c = start_el.Get_c(0,x);
    start_el.W1R = U_star + 4*local_c;
    RiemannStart(start_el, Q, ID);
}


// y is non-dimension, P is united
// qLnb have some issues
void Artery::Terminal_Flux(const int &n_step, const int &qLnb, const double &dt)
{
    int proceed = 1, iter = 0;
    int qLnb_1 = qLnb + 1;

    int Np = el[N_e-1].Np;
    Element &end_el = el[N_e-1];
    end_el.Set_W2();
    // Make sure that qLnb_1 runs in the interval [0:tmstps-1].
    if (qLnb_1 == (int) tmstps)
    {
        qLnb_1 = 0;
    }

    // In order to make a relation between P(x_L, t+dt) and Q(x_L, t+dt), and
    // we need to extract the term involving
    // y[0] (see mathematical derivation).
    // This term corresponds to the peripheral
    // The remaining terms in the convolution present at the boundary,
    // see mathematical derivation.
    double pterms = 0.0;


    double Unit_Y = q/rho/g/Lr;
    if (n_step > constants::tmstps)
    {
        for (int m=1; m<tmstps; m++)
        {
            int pindex  = (qLnb_1 + tmstps - m) % tmstps;
            pterms  = pterms  + (pL[pindex])*y[m];
        }
        pterms  = dt*pterms*Unit_Y;
//        pterms  = dt*pterms;
    }

    if( fabs(end_el.Q[Np-1]/end_el.A[Np-1]) > end_el.c[Np-1]){
        std::stringstream tmp;
        tmp << "Error in SmallTree outflow: flow is not subsonic.\n"
            << "Current element is: " << ID+1 << ".\n";
        throw std::runtime_error(tmp.str());
    }
    // The value of the function and the derivative is initialized to 0.
    double f,df,x,dx,c;
    x = end_el.A[Np-1];
    // y is non-dimension, need dimension
    while ((proceed)&&(iter++ < MAX_ITER))
    {
        c = end_el.Get_c(Np-1, x);
        double Q = dt*y[0]*Unit_Y*end_el.Get_P(Np-1, x)+pterms;
        f = Q/x+4*c-end_el.W2R;
//        df = dt*y[0]*Unit_Y*c*c*rho/pow(x,2.)-Q/pow(x,2.)-c/x;  // olufsen
        df = dt*y[0]*Unit_Y*c*c*rho/pow(x,2.)-Q/pow(x,2.)+c/x;
//        f = Q+4*c*x-end_el.W2R*x;
//        df = dt*y[0]*Unit_Y*c*c*rho/x+3*c-end_el.W2R*x;  // olufsen
        if( fabs(df) < SMALL ){
            throw std::runtime_error("Error in SmallTree outflow: zero derivative in the Newton's iteration.\n");
        }
        dx = f/df;
        x = x-dx;
        if (x <= 0.0)
        {
            std::cout << "WARNING (arteries.C): Bound_right: A was negative A = "
                      << x << " time = " << n_step*dt <<  " L = " << L << std::endl;
            x = end_el.A[Np-1]; // Bound xr away from zero.
        }
        if(dx*dx< TOL)      proceed = 0;
//        if(f*f < TOL)      proceed = 0;
    }
    // Solutions are applied, and right boundary and the intermediate array QL
    // are updated.
    if(iter >= MAX_ITER){
        std::stringstream tmp;
        tmp << "Error in terminal Riemann: iteration failed to converge. " << std::endl;
        tmp << "A is: " << x
                  << "; f is: " << f
                  << "; df is: " << df
                  << "; iter is: " << iter
                  << "; time is: " << dt*n_step
                  << std::endl;
        throw std::runtime_error(tmp.str());
    }
    double Qstar = (dt*y[0]*Unit_Y*end_el.Get_P(Np-1, x)+pterms);
    end_el.W2L = (Qstar/x-4*end_el.Get_c(Np-1,x));
    RiemannEnd(end_el,Qstar,ID);
}

void Artery::printX(std::ofstream &fd) const {
    double l = 0;
    fd.precision(15);
    fd << std::fixed;
    for (int i=0; i<N_e; i++)
    {
        int NumberOfPoint = el[i].Np;
        for (int j=0; j<NumberOfPoint; j++)
        {
            fd << el[i].VX[j] << "\n";
        }
    }
}

void Artery::printR0(std::ofstream &fd) const {
    double l = 0;
    fd.precision(15);
    fd << std::fixed;
    for (int i=0; i<N_e; i++)
    {
        int NumberOfPoint = el[i].Np;
        for (int j=0; j<NumberOfPoint; j++)
        {
            fd << el[i].r0[j] << "\n";
        }
    }
}

//void Artery::printUx(std::ofstream &fd) const {
//    double l = 0;
//    fd.precision(15);
//    fd << std::fixed;
//    for (int i=0; i<N_e; i++)
//    {
//        int NumberOfPoint = el[i].Np;
//        for (int j=0; j<NumberOfPoint; j++)
//        {
//            fd << el[i].U[j] << "\n";
//        }
//    }
//}
void Artery::printPx(std::ofstream &fd) const {
    double l = 0;
    fd.precision(15);
    fd << std::fixed;
    for (int i=0; i<N_e; i++)
    {
        int NumberOfPoint = el[i].Np;
        for (int j=0; j<NumberOfPoint; j++)
        {
            fd << el[i].P[j] << "\n";
        }
    }
}

void Artery::printQxt(std::ofstream &fd, const double &t) const {
    double l = 0;
    fd.precision(10);
    fd << std::fixed;
    fd << t << " " << l << " " << el[0].Q[0] << std::endl;
    l+=el[0].L;
    for (int i=1; i<N_e; i++)
    {
        fd << t << " " << l << " " << 0.5*(el[i-1].Q[ el[i-1].Np-1 ]+el[i].Q[0]) << std::endl;
        l+=el[i].L;
    }
    fd << t << " " << l << " " << el[N_e-1].Q[el[N_e-1].Np-1] << std::endl;
}

//void Artery::printUxt(std::ofstream &fd, const double &t) const {
//    double l = 0;
//    fd.precision(10);
//    fd << std::fixed;
//    fd << t << " " << l << " " << el[0].U[0] << std::endl;
//    l+=el[0].L;
//    for (int i=1; i<N_e; i++)
//    {
//        fd << t << " " << l << " " << 0.5*(el[i-1].U[ el[i-1].Np-1 ]+el[i].U[0]) << std::endl;
//        l+=el[i].L;
//    }
//    fd << t << " " << l << " " << el[N_e-1].U[el[N_e-1].Np-1] << std::endl;
//}

void Artery::printAxt(std::ofstream &fd, const double &t) const {
    double l = 0;
    fd.precision(10);
    fd << std::fixed;
    fd << t << " " << l << " " << el[0].A[0] << std::endl;
    l+=el[0].L;
    for (int i=1; i<N_e; i++)
    {
        fd << t << " " << l << " " << 0.5*(el[i-1].A[ el[i-1].Np-1 ]+el[i].A[0]) << std::endl;
        l+=el[i].L;
    }
    fd << t << " " << l << " " << el[N_e-1].A[el[N_e-1].Np-1] << std::endl;
}

void Artery::printPxt(std::ofstream &fd, const double &t) const {
    double l = 0;
    fd.precision(10);
    fd << std::fixed;
    fd << t << " " << l << " " << el[0].P[0] << std::endl;
    l+=el[0].L;
    for (int i=1; i<N_e; i++)
    {
        fd << t << " " << l << " " << 0.5*(el[i-1].P[ el[i-1].Np-1 ]+el[i].P[0]) << std::endl;
        l+=el[i].L;
    }
    fd << t << " " << l << " " << el[N_e-1].P[el[N_e-1].Np-1] << std::endl;
}

double Artery::Get_Q(const double &T){
    double q0   = 0.0;
    double qmax = 500;
    double a    = 2.0/3.0;
    double fi, qnew;
    if(T<=a) {
        fi = 3 * M_PI * T - sqrt(2);
        qnew = q0 + qmax * (0.251 + 0.290 * (cos(fi) + 0.97 * cos(2 * fi) + 0.47 * cos(3 * fi) + 0.14 * cos(4 * fi)));
    }else {
        fi = 3 * M_PI * a - sqrt(2);
        qnew = q0 + qmax * (0.251 + 0.290 * (cos(fi) + 0.97 * cos(2 * fi) + 0.47 * cos(3 * fi) + 0.14 * cos(4 * fi)));
    }
    return qnew;
}
//Qin_H = 4.0*atan(1.0)*0.0126*0.0126*(0.20617+0.37759*sin(2*PI*t/T+0.59605)
//        +0.2804*sin(4*PI*t/T-0.35859)+0.15337*sin(6*PI*t/T-1.2509)
//        -0.049889*sin(8*PI*t/T+1.3921)+0.038107*sin(10*PI*t/T-1.1068)
//        -0.041699*sin(12*PI*t/T+1.3985)-0.020754*sin(14*PI*t/T+0.72921)
//        +0.013367*sin(16*PI*t/T-1.5394)-0.021983*sin(18*PI*t/T+0.95617)
//        -0.013072*sin(20*PI*t/T-0.022417)+0.0037028*sin(22*PI*t/T-1.4146)
//        -0.013973*sin(24*PI*t/T+0.77416)-0.012423*sin(26*PI*t/T-0.46511)
//        +0.0040098*sin(28*PI*t/T+0.95145)-0.0059704*sin(30*PI*t/T+0.86369)
//        -0.0073439*sin(32*PI*t/T-0.64769)+0.0037006*sin(34*PI*t/T+0.74663)
//        -0.0032069*sin(36*PI*t/T+0.85926)-0.0048171*sin(38*PI*t/T-1.0306)
//        +0.0040403*sin(40*PI*t/T+0.28009)-0.0032409*sin(42*PI*t/T+1.202)
//        -0.0032517*sin(44*PI*t/T-0.93316)+0.0029112*sin(46*PI*t/T+0.21405)
//        -0.0022708*sin(48*PI*t/T+1.1869)-0.0021566*sin(50*PI*t/T-1.1574)
//        +0.0025511*sin(52*PI*t/T-0.12915)-0.0024448*sin(54*PI*t/T+1.1185)
//        -0.0019032*sin(56*PI*t/T-0.99244)+0.0019476*sin(58*PI*t/T-0.059885)
//        -0.0019477*sin(60*PI*t/T+1.1655)-0.0014545*sin(62*PI*t/T-0.85829)
//        +0.0013979*sin(64*PI*t/T+0.042912)-0.0014305*sin(66*PI*t/T+1.2439)
//        -0.0010775*sin(68*PI*t/T-0.79464)+0.0010368*sin(70*PI*t/T-0.0043058)
//        -0.0012162*sin(72*PI*t/T+1.211)-0.00095707*sin(74*PI*t/T-0.66203)
//        +0.00077733*sin(76*PI*t/T+0.25642)-0.00092407*sin(78*PI*t/T+1.3954)
//        -0.00079585*sin(80*PI*t/T-0.49973));

// save the time field resistance of small structure tree
void Artery::Save_inverse_Z(std::string &filename){
    std::ofstream s(filename);

    if (!s.is_open()) {
        throw std::runtime_error("Couldn't open file... "+filename +" when saving the "+filename+"!\n");
    }else{
        std::cout << "Saving File " << filename << "......" << std::endl;
    }

    s.precision(17);
    s << std::fixed;
    for(int i=0;i<tmstps;++i){
        s << y[i] << "\n";
    }
    s.close();
    if (s.fail()) {
        throw std::runtime_error(filename +" close failed when saving the file!"+"\n");
    }else{
        std::cout << "Saving " << filename << " File is OK!" << std::endl;
    }
}

// read the time field resistance of small structure tree from cache
void Artery::Read_inverse_Z(std::string &filename){
    std::ifstream s(filename);
    if (!s.is_open()) {
        throw std::runtime_error("Couldn't open file... "+filename +" when reading the "+filename+"!\n");
    }else{
        std::cout << "Reading File " << filename << "......" << std::endl;
    }
    for(int i=0;i<tmstps;++i){
        s >> y[i];
    }

    s.close();
    if (s.fail()) {
        throw std::runtime_error(filename +" close failed!\n");
    }else{
        std::cout << "Reading " << filename << " File is OK!" << std::endl;
    }
}



//B ioFlux(nbrves, Arteries, ID_Bif, ID_Out, qLnb,  dt,  rk4c);
void BioFlux(const int &nbrves, Artery *Arteries[],
             const std::set<int>& ID_Bif, const std::set<int>& ID_Out,
             const int &n_step, const int &qLnb, const double &dt)
{
    // Q0[i], i should get from mod
    // Q0 change with time can modified by add function not read from file
    // the first element of inlet, and dont need check Q and P,
    Arteries[0]->Inlet_Flux((qLnb)*dt);
    for (auto i: ID_Bif) {
        Arteries[i]->Bifur_Flux();
    }
    for (auto i: ID_Out) {
        Arteries[i]->Terminal_Flux(n_step, qLnb, dt);
    }
    for (int i=0;i<nbrves;++i)
    {
        Arteries[i]->Inter_Flux();
    }

}
//
//
void RiemannStart(Element &start_el, const double &Q_star, const int &ID){
    // then use 4 value to solve the A1,A2,U1,U2
    double cl,cr,k1,k2,k;
//    double x[4], f[4], dx[4];
//    double inv_J[4][4];
    VectorXd x(4), dx(4), f(4);
    MatrixXd J(4,4);
    int proceed =1,iter=0;
    if( (fabs(start_el.Q[0]/start_el.A[0]) > start_el.c[0]) )
    {
        std::string Error = "Error in IntletFlux: Flow is not subsonic in domain: "+std::to_string(ID+1)+".\n";
        throw std::runtime_error(Error);
    }
    x[0] = start_el.Q[0];
    x[1] = 2*Q_star-start_el.Q[0];
    x[2] = start_el.A[0];
    x[3] = start_el.A[0];
    while((proceed)&&(iter++ < MAX_ITER)){
        // Calculate constraint vector and the wave speed at each vessel.

        cl   = start_el.Get_c(0, x[2]);
        cr   = start_el.Get_c(0, x[3]);
        // x(0): UL, x(1): UR, x(2): AL, x(3): AR
        // Inverse Jacobian matrix dfdv

        f[0] = x[0]/x[2] + 4*cl - start_el.W1R;
        f[1] = x[1]/x[3] - 4*cr - start_el.W1L;
        f[2] = x[0]-x[1];
        f[3] = 0.5*rho*(x[0]/x[2],2.)+start_el.Get_P(0,x[2])
              -0.5*rho*(x[1]/x[3],2.)-start_el.Get_P(0,x[3]);

        // x(0): UL, x(1): UR, x(2): AL, x(3): AR
        // Inverse Jacobian matrix dfdv

        J(0,0) = 1/x[2];
        J(0,1) = 0;
//        J(0,2) = -x[0]/pow(x[2],2.)-cl/x[2];    // olufsen
        J(0,2) = -x[0]/pow(x[2],2.)+cl/x[2];
        J(0,3) = 0;

        J(1,0) = 0;
        J(1,1) = 1/x[3];
        J(1,2) = 0;
//        J(1,3) = -x[1]/pow(x[3],2.)+cr/x[3];    // olufsen
        J(1,3) = -x[1]/pow(x[3],2.)-cr/x[3];

        J(2,0) =  1;
        J(2,1) = -1;
        J(2,2) = 0;
        J(2,3) = 0;

        J(3,0) =  rho*x[0]/pow(x[2],2.);
        J(3,1) = -rho*x[1]/pow(x[3],2.);
        J(3,2) = -rho*( pow(x[0],2.)/pow(x[2],3.)-pow(cl,2.)/x[2] );
        J(3,3) =  rho*( pow(x[1],2.)/pow(x[3],3.)-pow(cr,2.)/x[3] );

        /* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
        dx = J.partialPivLu().solve(f);
        // Update solution: x_new = x_old - dx
        x = x-dx;

        // Check if the error of the solution is smaller than TOL.
//        if( dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3] < TOL )      proceed = 0;
        if( (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]) < TOL )      proceed = 0;
    }
    start_el.A1 = x[3];
    start_el.Q1 = x[1];
    start_el.Set_F1();
}

void RiemannEnd(Element &end_el, const double &Q_star, const int &ID){
    // then use 4 value to solve the A1,A2,U1,U2
    double cl,cr,k1,k2,k;
//    double x[4], f[4], dx[4];
//    double inv_J[4][4];
    VectorXd x(4), dx(4), f(4);
    MatrixXd J(4,4);
    int proceed =1,iter=0;
    int Np = end_el.Np;
    if( (fabs(end_el.Q[Np-1]/end_el.A[Np-1]) > end_el.c[Np-1]) )
    {
        std::string Error = "Error in outletFlux: Flow is not subsonic in domain: "+std::to_string(ID+1)+".\n";
        throw std::runtime_error(Error);
    }
    x[0] = end_el.Q[Np-1];
    x[1] = 2*Q_star-end_el.Q[Np-1];
    x[2] = end_el.A[Np-1];
    x[3] = end_el.A[Np-1];
    while((proceed)&&(iter++ < MAX_ITER)){
        // Calculate constraint vector and the wave speed at each vessel.

        cl   = end_el.Get_c(Np-1, x[2]);
        cr   = end_el.Get_c(Np-1, x[3]);

        // x(0): UL, x(1): UR, x(2): AL, x(3): AR
        // Inverse Jacobian matrix dfdv

        f[0] = x[0]/x[2] + 4*cl - end_el.W2R;
        f[1] = x[1]/x[3] - 4*cr - end_el.W2L;
        f[2] = x[0]-x[1];
        f[3] = 0.5*rho*(x[0]/x[2],2.)+end_el.Get_P(Np-1,x[2])
              -0.5*rho*(x[1]/x[3],2.)-end_el.Get_P(Np-1,x[3]);

        // x(0): UL, x(1): UR, x(2): AL, x(3): AR
        // Inverse Jacobian matrix dfdv

        J(0,0) = 1/x[2];
        J(0,1) = 0;
//        J(0,2) = -x[0]/pow(x[2],2.)-cl/x[2];    // olufsen
        J(0,2) = -x[0]/pow(x[2],2.)+cl/x[2];
        J(0,3) = 0;

        J(1,0) = 0;
        J(1,1) = 1/x[3];
        J(1,2) = 0;
//        J(1,3) = -x[1]/pow(x[3],2.)+cr/x[3];    // olfusen
        J(1,3) = -x[1]/pow(x[3],2.)-cr/x[3];

        J(2,0) =  1;
        J(2,1) = -1;
        J(2,2) = 0;
        J(2,3) = 0;

        J(3,0) =  rho*x[0]/pow(x[2],2.);
        J(3,1) = -rho*x[1]/pow(x[3],2.);
        J(3,2) = -rho*( pow(x[0],2.)/pow(x[2],3.)-pow(cl,2.)/x[2] );
        J(3,3) =  rho*( pow(x[1],2.)/pow(x[3],3.)-pow(cr,2.)/x[3] );

        /* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
        dx = J.partialPivLu().solve(f);
        // Update solution: x_new = x_old - dx
        x = x-dx;

        // Check if the error of the solution is smaller than TOL.
//        if( dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3] < TOL )      proceed = 0;
        if( (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]) < TOL )      proceed = 0;
    }
    end_el.A2 = x[2];
    end_el.Q2 = x[0];
    end_el.Set_F2();
}

void Riemann(Element &el1, Element &el2, const int &j, const int &N_e)
{
    // f: function g: delta_x,
    // x: x_new = x_old+delta_x
    // inv_j: inv of f'
//    long double f[4],dx[4],x[4];
//    long double inv_J[4][4];
    VectorXd f(4), dx(4), x(4);
    MatrixXd J(4,4);
    long double cl, cr, k1, k2, k;
    int proceed = 1, iter = 0;
    // left element, use right edge use W2R
    el1.Set_W2();
    // right element, use left edge use W1L
    el2.Set_W1();
    int N_ = el1.Np-1;
    // Check if flow is subsonic
    if( (fabs(el1.Q[N_]/el1.A[N_]) > el1.c[N_]) || (fabs(el2.Q[0]/el2.A[0]) > el2.c[0]) )
    {
        std::stringstream tmp;
        tmp << "Error in InterFlux: Flow is not subsonic. \n"
            << "Element is: " << j-1 << " and " << j << ".\n"
            << "The number of element is: " << N_e << ".\n";
        throw std::runtime_error(tmp.str());
    }
    x[0] = el1.Q[N_];
    x[1] = el2.Q[0];
    x[2] = el1.A[N_];
    x[3] = el2.A[0];
    while((proceed)&&(iter++ < MAX_ITER)){
        // Calculate constraint vector and the wave speed at each vessel.

        cl   = el1.Get_c(N_, x[2]);
        cr   = el2.Get_c(0, x[3]);

        f[0] = x[0]/x[2] + 4*cl - el1.W2R;
        f[1] = x[1]/x[3] - 4*cr - el2.W1L;
        f[2] = x[0]-x[1];
        f[3] = 0.5*rho*(x[0]/x[2],2.)+el1.Get_P(N_,x[2])
              -0.5*rho*(x[1]/x[3],2.)-el2.Get_P(0,x[3]);

        // x(0): UL, x(1): UR, x(2): AL, x(3): AR
        // Inverse Jacobian matrix dfdv

        J(0,0) = 1/x[2];
        J(0,1) = 0;
//        J(0,2) = -x[0]/pow(x[2],2.)-cl/x[2];  // olufsen
        J(0,2) = -x[0]/pow(x[2],2.)+cl/x[2];
        J(0,3) = 0;

        J(1,0) = 0;
        J(1,1) = 1/x[3];
        J(1,2) = 0;
//        J(1,3) = -x[1]/pow(x[3],2.)+cr/x[3];  // olufsen
        J(1,3) = -x[1]/pow(x[3],2.)-cr/x[3];

        J(2,0) =  1;
        J(2,1) = -1;
        J(2,2) = 0;
        J(2,3) = 0;

        J(3,0) =  rho*x[0]/pow(x[2],2.);
        J(3,1) = -rho*x[1]/pow(x[3],2.);
        J(3,2) = -rho*( pow(x[0],2.)/pow(x[2],3.)-pow(cl,2.)/x[2] );
        J(3,3) =  rho*( pow(x[1],2.)/pow(x[3],3.)-pow(cr,2.)/x[3] );

        /* Solve the linear system by inverting analyticaly the Jacobian: g = (dfdv)^(-1)*f */
        dx = J.partialPivLu().solve(f);
        // Update solution: x_new = x_old - dx
        x = x-dx;

        // Check if the error of the solution is smaller than TOL.
        if( (dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3]) < TOL )      proceed = 0;
//        if( (f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3]) < TOL )      proceed = 0;
    }

    if(iter >= MAX_ITER){
        std::stringstream tmp;
        tmp << "Error in Riemann: iteration failed to converge. \n"
            << "The L2 norm is: " << dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]+dx[3]*dx[3] << ".\n"
//                << "The L2 norm is: " << f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3] << ".\n"
            << "Element is: " << j-1 << " and " << j << ".\n";
        throw std::runtime_error(tmp.str());
    }
    el1.Q2 = double(x[0]);              // left element U*L
    el2.Q1 = double(x[1]);              // right element U*R
    el1.A2 = double(x[2]);              // left element A*L
    el2.A1 = double(x[3]);              // right element A*R
}
//


void solver(const int &nbrves, Artery *Arteries[], int &n_step,
            const double &tstart, const double &tend, const double &dt,
            const std::set<int>& ID_Out, const std::set<int>& ID_Bif)
{

    // Runge-Kutta coefficient.
//    using namespace RK4;
    // start time
    double t = tstart;
    // used for terminal boundary and inlet boundary, the n_step in a Period
    // qLnb to interval [0, tmstps-1]
    int qLnb = n_step % constants::tmstps;

    while (t < tend)
    {
        // Check that the step we take is valid. If this error occurs when forcing
        // a constant step-size the program must be terminated.
        if (t+dt > tend)
        {
            std::cout << "ERROR (arteries.C): Step-size changed, t+k="
                      << t+dt << ", tend=" << tend << ", dt=" << tend-t
                      << ", dt_old=" << dt << std::endl;
            throw std::runtime_error("Step-size changed!");
        }


        // Check that the CFL-condition applies.
        // CFL condition mean small h the k should also be small.
        double current_CFL = Get_CFL(nbrves, Arteries);
        double courant = dt/current_CFL;
        if (courant > 1)
        {
            std::cout << "Step: " << n_step << ", Time: " << t << ", CFL: " << courant << std::endl;
            throw std::runtime_error("arteries.C, Step-size too large CFL-condition violated.");
        }
        if(!(n_step % constants::Per_step)){
            std::cout << "Step: " << n_step << ", Time: " << t << ", CFL: " << courant << std::endl;
//            for (int i = 0; i < nbrves; ++i) {
//                std::ofstream output;
//                std::stringstream tmp;
//                tmp << "Artery_" << std::setfill('0') << std::setw(2) << i+1
//                    << "_P_" << std::setfill('0') << std::setw(8) << n_step+1 << ".dat";
//                std::string Uname = tmp.str();
//                output.open(Uname);
//                Arteries[i]->printPx(output);
//                output.close();
//            }
        }
        for(int j=0;j<RK4::N_t;++j) {
            try {
//                solverRHS(nbrves, Arteries, ID_Bif, ID_Out, n_step, qLnb, dt, 0);
                solverRHS(nbrves, Arteries, ID_Bif, ID_Out, n_step, qLnb, dt);
//                for (int i = 0; i < nbrves; ++i) {
//                    // calculate the right hand side of ODE
//                    // update
//                    Arteries[i]->Update(rk4a[j], rk4b[j], dt);
////                    Arteries[i]->Update(0, 1, dt);
//                }
                for (auto i: ID_Out){
                    Arteries[i]->Update_pL(qLnb+1);
                }

            } catch (std::exception &e) {
                std::cerr << "----------------------------------------------------"
                          << std::endl
                          << std::endl
                          << "Some error happen at Run time step: " << n_step
                          << ", time is: " << t << ".\n CFL: " << courant << ", Runge-Kutta step is: " << j << ".\n";
                std::cerr << e.what() << "\n"
                          << "----------------------------------------------------"
                          << std::endl;
                exit(0);
            } catch (...) {
                std::cerr << std::endl
                          << std::endl
                          << "----------------------------------------------------"
                          << std::endl;
                std::cerr << "Unknown exception at Run time step: " << n_step
                          << ", time is: " << t << ".\n CFL: " << courant << ", Runge-Kutta step is: " << j << ".\n"
                          << "Aborting!" << std::endl
                          << "----------------------------------------------------"
                          << std::endl;

                exit(0);
            }
        }
//        for (int i = 0; i < nbrves; ++i) {
//            //
//            Arteries[i]->clean_k();
//        }
        ++n_step;
        t = t + dt;
        qLnb = n_step % tmstps;
    }

}
////
void solverRHS(const int &nbrves, Artery *Arteries[],
               const std::set<int>& ID_Bif, const std::set<int>& ID_Out,
               const int &n_step, const int &qLnb, const double &dt)
{
    // \frac{\partial F}{\partial x}
    // upwind Flux-local Flux
    // source
    for (int i=0;i<nbrves;++i)
    {
        Arteries[i]->Set_F(); // set c, P and then set F
        Arteries[i]->Set_S(); // set source term
    }
    BioFlux(nbrves, Arteries, ID_Bif, ID_Out, n_step, qLnb, dt);
    for (int i=0;i<nbrves;++i)
    {
        Arteries[i]->Set_RHS();
    }
}
//
double Get_CFL(const int &nbrves, Artery *Arteries[])
{
    Arteries[0]->local_CFL();
    double Courant = Arteries[0]->Courant;
    for (int i=1; i<nbrves; i++)
    {
        Arteries[i]->local_CFL();
        Courant = std::min(Arteries[i]->Courant, Courant);
    }
    return Courant;
}


