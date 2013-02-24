#include "energytest.h"
using namespace arma;
EnergyTest::EnergyTest()
{
    vec dt_vec(5);
    dt_vec<<0.0<<0.01<<0.001;
    double dt;
    for(int i=0;i<dt_vec.size();i++){
        string s;
        stringstream out;
        out<<i;
        s= out.str();
        dt = dt_vec[i];
        runFor(dt,300,s);
    }

}

void EnergyTest::runFor(double dt, int numOfTimeSteps, string trial) {
    double L0 = 3.405;
    double r_cut = 3*L0;
    int Nx = 8, Ny = 8, Nz = 8;
    int CellNx, CellNy, CellNz;
    double b = 5.26;
    double T = 100;
    CellNx = Nx*b/r_cut;
    CellNy = CellNx;
    CellNz = CellNx;
    r_cut = Nx*b/(CellNx);
    //double maxTime = 1;
    string element = "Ar";
    string filenamebase = "energyTest";
    //string current_trial_string;
    //stringstream out;
    bool writeVMD = false;
    bool writeMeasurements = true;
    double T_bath = 1;
    bool Andersen = false;
    bool Berendsen = false;

    CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    mySolver->solve(0,numOfTimeSteps, dt,filenamebase+trial, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);

}
