#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <chrono>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <filesystem>
#include "Jacobi1D.h"
#include "Element.h"
#include "artery.h"
#include <memory>
#include "tools.h"
#include "Parameters.h"
#include "test.h"
using namespace Eigen;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::set;
using std::string;
using std::ofstream;
using namespace chrono;

using constants::u;
using constants::Lr;
using constants::tmstps;
using constants::plts;
using constants::FileName1;
using constants::FileName2;
using constants::ff1;
using constants::ff2;
using constants::ff3;
using constants::fa1;
using constants::fa2;
using constants::fa3;
using constants::point;

extern "C"  void impedance_init_driver_(int *tmstps);

int main() {

//    double Tper   = 1.0;                    // The period of one heart beat [s].
    double Period = 1.0;    //Tper*u/Lr;              // The dimension-less period.
    double dt     = Period/tmstps;          // Length of a timestep.
    double Deltat = Period/plts;            // Interval between each point
    double tstart, tend, finaltime;
    double rm4    = 0.03;                           // cm

    int    nbrves = 55;             // The number of vessels in the network.
    int    n_step = 0;                      // solver time counter.
    // get the id of the topology
    std::cout << "----------------------------------------------------\n\n";
    std::cout << "The Gauss-Lobatto-Legendre point is: \n";
    std::cout << constants::jac.r << std::endl;
    std::cout << "The space first derivative matrix is: \n";
    std::cout << constants::jac.Dr << std::endl;
    std::cout << "The transfer matrix between physic space and legendre space is: \n";
    std::cout << constants::jac.V << std::endl;
    std::cout << "The inverse mass matrix is: \n";
    std::cout << constants::jac.invM << std::endl;
    std::cout << "----------------------------------------------------\n" << std::endl;
    set<int> ID_Out, ID_Bif;
//    ID_Out.insert(0);
    vector<vector<double> > Geometry(nbrves, vector<double> (3));
    vector<vector<int> > Topology(nbrves, vector<int> (5));
    get_ID(FileName1, ID_Out, ID_Bif, Topology);
    get_Geo(FileName2, Geometry);
    // Creat the result dir.
    string Path = "result";
    std::filesystem::path dir{Path};
    if(!std::filesystem::exists(dir)) std::filesystem::create_directory(dir);
    // result file name
    vector<string> nameP(nbrves, Path+string("/p"));
    vector<string> nameQ(nbrves, Path+string("/q"));
    vector<string> nameA(nbrves, Path+string("/A"));
    for(int i=0;i<nbrves;i++)
    {
        auto j=i+1;
        nameP[i] += to_string(j);
        nameQ[i] += to_string(j);
        nameA[i] += to_string(j);
    }

    ofstream fp[nbrves];
    ofstream fq[nbrves];
    ofstream fA[nbrves];
    for(int i=0;i<nbrves;i++)
    {
        fp[i].open(nameP[i], ofstream::out);
        if (!fp[i].is_open()) {
            cerr << "Couldn't open file... " << nameP[i] << "!" << endl;
            return -1;
        }else{
            cout << "Open File " << nameP[i] << " is OK!" << endl;
        }
        fq[i].open(nameQ[i], ofstream::out);
        if (!fq[i].is_open()) {
            cerr << "Couldn't open file... " << nameQ[i] << "!" << endl;
            return -1;
        }else{
            cout << "Open File " << nameQ[i] << " is OK!" << endl;
        }
        fA[i].open(nameA[i], ofstream::out);
        if (!fA[i].is_open()) {
            cerr << "Couldn't open file... " << nameA[i] << "!" << endl;
            return -1;
        }else{
            cout << "Open File " << nameA[i] << " is OK!" << endl;
        }
    }

    // read geometry and create the artery class

    auto c1 = system_clock::now();        // Only used when timing the program.
    tstart     = 0.0;            // Starting time.
    finaltime  = 10*Period;       // Final end-time during a simulation.
    tend       = 0*Period;       // Timestep before the first plot-point
    // is reached.
    int AST = constants::tmstps;
    impedance_init_driver_(&AST);
    Artery *Arteries[nbrves];
    for(int i=nbrves-1;i>=0;--i)
    {
        auto L  = Geometry[i][0];
        auto rt = Geometry[i][1];
        auto rb = Geometry[i][2];
        double rmin = (Arteries[Topology[i][1]]==0 && Arteries[Topology[i][2]] ==0 ) ? rm4: 0;

        if (i)
            if (Topology[i][1]==0 && Topology[i][2]==0)
                Arteries[i] = new Artery(L,rt,rb,i,nullptr,nullptr,nullptr,nullptr,rmin,point,0,ff1,ff2,ff3,fa1,fa2,fa3,0.0, Period);
            else
                Arteries[i] = new Artery(L,rt,rb,i,Arteries[Topology[i][1]-1],Arteries[Topology[i][2]-1],nullptr,nullptr,rmin,point,0,ff1,ff2,ff3,fa1,fa2,fa3,0.0, Period);
        else
            Arteries[i] = new Artery(L,rt,rb,i,Arteries[Topology[i][1]-1],Arteries[Topology[i][2]-1],nullptr,nullptr,rmin,point,1,ff1,ff2,ff3,fa1,fa2,fa3,0.0, Period);
    }
//    test_arteries(Arteries[0]);
    for (int i = 0; i < nbrves; ++i) {
        std::ofstream output;
        std::stringstream tmp;
        tmp << "Artery_" << std::setfill('0') << std::setw(2) << i+1
            << "_X.dat";
        std::string Xname = tmp.str();
        output.open(Xname);
        Arteries[i]->printX(output);
        output.close();
    }
    for (int i = 0; i < nbrves; ++i) {
        std::ofstream output;
        std::stringstream tmp;
        tmp << "Artery_" << std::setfill('0') << std::setw(2) << i+1
            << "_r0.dat";
        std::string Xname = tmp.str();
        output.open(Xname);
        Arteries[i]->printR0(output);
        output.close();
    }
//    Arteries[0] = new Artery(50,1.5,0.3,0,nullptr,nullptr,nullptr,nullptr,rm4,point,1,ff1,ff2,ff3,fa1,fa2,fa3,0.0, Period);
//    solver(nbrves, Arteries, n_step,
//           tstart, tend, dt, ID_Out, ID_Bif);


    // Runge-Kutta coefficient.
    using namespace RK4;
    // start time
    double t = tstart;
    // used for terminal boundary and inlet boundary, the n_step in a Period
    // qLnb to interval [0, tmstps-1]
    int qLnb = n_step % constants::tmstps;

    while (t < finaltime)
    {

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
            if(t>=tend) {
                for (int i = 0; i < nbrves; ++i) {
                    Arteries[i]->printPxt(fp[i], t);
                    Arteries[i]->printQxt(fq[i], t);
                    Arteries[i]->printAxt(fA[i], t);
                }
            }
        }
        for(int j=0;j<RK4::N_t;++j) {
            try {
//                solverRHS(nbrves, Arteries, ID_Bif, ID_Out, n_step, qLnb, dt, 0);
                solverRHS(nbrves, Arteries, ID_Bif, ID_Out, n_step, qLnb, dt, rk4c[j]);
                for (int i = 0; i < nbrves; ++i) {
                    // calculate the right hand side of ODE
                    // update
                    Arteries[i]->Update(rk4a[j], rk4b[j], dt);
//                    Arteries[i]->Update(0, 1, dt);
                }
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


    auto c2 = system_clock::now();
    auto duration = duration_cast<microseconds>(c2-c1);
    cout << "Number of arteries: " << nbrves <<endl;
    cout << "CPU-time : "
         << double(duration.count()) * microseconds::period::num / microseconds::period::den
         << "s" << endl;
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    // close the result files.
    for(int i=0;i<nbrves;++i)
    {
        fp[i].close();
        if (fp[i].fail()) {
            cerr << nameP[i] << "failed " << "!" << endl;
            return -1;
        }
        else{
            cout << "Saving File " << nameP[i] << " is OK!" << endl;
        }
        fq[i].close();
        if (fq[i].fail()) {
            cerr << nameQ[i] << "failed " << "!" << endl;
            return -1;
        }
        else{
            cout << "Saving File " << nameQ[i] << " is OK!" << endl;
        }
        fA[i].close();
        if (fA[i].fail()) {
            cerr << nameA[i] << "failed " << "!" << endl;
            return -1;
        }
        else{
            cout << "Saving File " << nameA[i] << " is OK!" << endl;
        }
    }

    return 0;
}
