//
// Created by chen on 2021/12/21.
//

#include <iostream>
#include <fstream>
#include <cmath>
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
using namespace std;

int main(int argc, char *argv[]){
    int tmstps = 32768;
    if(argc == 1) cout << "Using default number of time point: " << tmstps << "." << endl;
    if(argc > 2) cerr << "Use: \"" << argv[0] << "\" or \"" << argv[0] << " " << "Number\""
                 << "\n to run the program." << endl;
    if(argc == 2) {
        tmstps = stoi(argv[1]);
    }
    double Period = 1.;
    double dt = Period/tmstps;

    auto f = [](const double &t){
        double q0   = 0.0;
        double qmax = 500;
        double a    = 2.0/3.0;
        double fi, qnew;
        if(t<=a) {
            fi = 3 * M_PI * t - sqrt(2);
            qnew = q0 + qmax * (0.251 + 0.290 * (cos(fi) + 0.97 * cos(2 * fi) + 0.47 * cos(3 * fi) + 0.14 * cos(4 * fi)));
        }else {
            fi = 3 * M_PI * a - sqrt(2);
            qnew = q0 + qmax * (0.251 + 0.290 * (cos(fi) + 0.97 * cos(2 * fi) + 0.47 * cos(3 * fi) + 0.14 * cos(4 * fi)));
        }
        return qnew;
    };
    ofstream ss("input.dat");
    ss.precision(15);
    ss << std::fixed;
    for(int i=0;i<=tmstps;++i){
        double t = i*dt;
        ss << f(t) << "\n";
    }
    ss.close();
    return 0;
}