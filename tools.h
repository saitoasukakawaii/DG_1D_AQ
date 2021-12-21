//
// Created by chen on 2021/9/16.
//

#ifndef DG_1D_BURGERS_TOOLS_H
#define DG_1D_BURGERS_TOOLS_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>
#include <vector>
#include "tools.h"
#include <filesystem>
using namespace std;

template<typename T>
void print_ln(T u)
{
    for(auto x: u)
        cout << x << "\t";
    cout << endl;
    return ;
}
void get_ID(const string &FileName, set<int>& ID_Out, set<int>& ID_Bif, vector<vector<int> >&Topology)
{
    ID_Out.clear();
    ID_Bif.clear();
    int numOfArtery, LD, RD, LP, RP;
    string line;
    ifstream VascularTopoloy(FileName);
    if(!VascularTopoloy.is_open()) throw std::runtime_error(string{"Could not open file: "}+FileName);
    int i=0;
    while (getline(VascularTopoloy, line))
    {
        stringstream ss(line);
        ss >> numOfArtery;
        ss >> LD;
        ss >> RD;
        ss >> LP;
        ss >> RP;
        Topology[i] = vector{numOfArtery, LD,RD,LP,RP};
        numOfArtery = numOfArtery-1;
        if (LD!=RD) ID_Bif.insert(numOfArtery);
        else if (LD==0 && RD==0) ID_Out.insert(numOfArtery);
        ++i;
    }
    VascularTopoloy.close();
}

void get_Geo(const string &FileName, vector<vector<double> >&Geometry)
{

    double numOfArtery, rt, rb;
    string line;
    ifstream VascularTopoloy(FileName);
    if(!VascularTopoloy.is_open()) throw std::runtime_error(string{"Could not open file: "}+FileName);
    int i=0;
    while (getline(VascularTopoloy, line))
    {
        stringstream ss(line);
        ss >> numOfArtery;
        ss >> rt;
        ss >> rb;
        numOfArtery = numOfArtery;
        Geometry[i] = vector{numOfArtery,rt,rb};
        ++i;
    }
    VascularTopoloy.close();
}
#endif //DG_1D_BURGERS_TOOLS_H
