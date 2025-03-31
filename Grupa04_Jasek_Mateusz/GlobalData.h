#ifndef GLOBALDATA_H
#define GLOBALDATA_H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

struct GlobalData {
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int nN;
    int nE;
    int Npc;
    GlobalData(int sTime, int ssTime, int con, int alf, int to, int iTemp, int den, int sHeat, int n, int e, int np)
        : SimulationTime(sTime), SimulationStepTime(ssTime), Conductivity(con), Alfa(alf), Tot(to),
        InitialTemp(iTemp), Density(den), SpecificHeat(sHeat), nN(n), nE(e),Npc(np) {}
};

GlobalData GloadData();
void displayData(GlobalData newData);
#endif // GLOBALDATA_H
