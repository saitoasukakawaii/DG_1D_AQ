//
// Created by chen on 2021/11/16.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include "tools.h"
#include <filesystem>
#include <Eigen/Dense>
#include "artery.h"
using namespace Eigen;
using std::endl;
using std::cout;
//using namespace std;

//void get_ID(const string &FileName, set<int>& ID_Out, set<int>& ID_Bif, vector<vector<int> >&Topology);
//void get_Geo(const string &FIleName, vector<vector<double> >&Geometry);
void test_arteries(Artery *Arteries)
{

}
//int main(int argc, char *argv[])
//{
//    for (int i = 0; i < argc; ++i)
//        cout << argv[i] << "\n";
//    string filename;
//    if (argc > 2) throw runtime_error("Too many args.");
//    if (argc == 1)  filename = "Geometry/topology55.txt";
//    else filename = argv[1];
//    string FileName2 = "Geometry/geometry55.txt";
//    vector<vector<double> > Geometry(55, vector<double> (3));
//    vector<vector<int> > Topology(55, vector<int> (5));
//    get_Geo(FileName2, Geometry);
//    for (int i=0;i<55;++i)
//    {
//        cout << Geometry[i][0] << ", "
//             << Geometry[i][1] << ", "
//             << Geometry[i][2] << endl;
//    }
//    set<int> ID_Out, ID_Bif;
//    get_ID(filename, ID_Out, ID_Bif, Topology);
//    for (int i=0;i<55;++i)
//    {
//        cout << Topology[i][0] << ", "
//             << Topology[i][1] << ", "
//                << Topology[i][2] << ", "
//                << Topology[i][3] << ", "
//             << Topology[i][4] << endl;
//    }
//    print_ln(ID_Out);
//    print_ln(ID_Bif);
//    string Path = "result";
//    filesystem::path dir{Path};
//    cout << dir << endl;
//    if(!filesystem::exists(dir)) filesystem::create_directory(dir);
//    ifstream fi("input.dat");
//    VectorXd Q0(16385);
//    for(int i=0;i<16385;++i)
//    {
//        fi >> Q0[i];
//    }
//    print_ln(Q0);
//    return 0;
//}

//void get_ID(const string &FileName, set<int>& ID_Out, set<int>& ID_Bif, vector<vector<int> >&Topology)
//{
//    ID_Out.clear();
//    ID_Bif.clear();
//    int numOfArtery, LD, RD, LP, RP;
//    string line;
//    ifstream VascularTopoloy(FileName);
//    if(!VascularTopoloy.is_open()) throw std::runtime_error(string{"Could not open file: "}+FileName);
//    int i=0;
//    while (getline(VascularTopoloy, line))
//    {
//        stringstream ss(line);
//        ss >> numOfArtery;
//        ss >> LD;
//        ss >> RD;
//        ss >> LP;
//        ss >> RP;
//        numOfArtery = numOfArtery-1;
//        if (LD!=RD) ID_Bif.insert(numOfArtery);
//        else if (LD==0 && RD==0) ID_Out.insert(numOfArtery);
//        Topology[i] = vector{numOfArtery, LD,RD,LP,RP};
//        ++i;
//    }
//    VascularTopoloy.close();
//}
//
//void get_Geo(const string &FileName, vector<vector<double> >&Geometry)
//{
//
//    double numOfArtery, rt, rb;
//    string line;
//    ifstream VascularTopoloy(FileName);
//    if(!VascularTopoloy.is_open()) throw std::runtime_error(string{"Could not open file: "}+FileName);
//    int i=0;
//    while (getline(VascularTopoloy, line))
//    {
//        stringstream ss(line);
//        ss >> numOfArtery;
//        ss >> rt;
//        ss >> rb;
//        numOfArtery = numOfArtery-1;
//        Geometry[i] = vector{numOfArtery,rt,rb};
//        ++i;
//    }
//    VascularTopoloy.close();
//}