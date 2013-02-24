#include "timetest.h"

TimeTest::TimeTest()
{

    double L0 = 3.405;
    double r_cut = 3*L0;
    int Nx = 8, Ny = 8, Nz = 8;
    int CellNx, CellNy, CellNz;
    double b = 5.26;
    double T = 100;
    CellNx = Nx*b/r_cut;
    //CellNx = 1;
    CellNy = CellNx;
    CellNz = CellNx;
    r_cut = Nx*b/(CellNx);
    string element = "Ar";

    vec time_steps_vec(3);
    time_steps_vec<<100<<500<<1000;
    vec N_vec(3);
    N_vec<<8<<12<<16;
    double numOfTimeSteps;
    bool writeVMD = false;
    bool writeMeasurements = false;
    //int N;
    double dt = 0.01;
    string filenamebase = "test";
    double mass =1;
    int max_time_step = numOfTimeSteps;
    for(int i=0;i<3;i++){
        Nx = N_vec(i);
        CellNx = Nx*b/r_cut;
        //CellNx = 1;
        CellNy = CellNx;
        CellNz = CellNx;
        r_cut = Nx*b/(CellNx);
        Ny= Nx;
        Nz = Nx;
        for(int j=0;j<3;j++){

            numOfTimeSteps = time_steps_vec(j);
            int max_time_step = numOfTimeSteps;
            CellSolver* mySolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element);
            cout<<"N: "<<Nx<<" numOfTimesteps: "<<numOfTimeSteps<<endl;
            mySolver->solve(0,numOfTimeSteps, dt,filenamebase, writeVMD, writeMeasurements);



            Lattice mol = Lattice(Nx,Nx,Nx, element,b, T, mass);
            Verlet_solver* vSolver = new Verlet_solver(mol);
            //vSolver->solve(0,1000,0.01,filenamebase);
            vSolver->solve(0,numOfTimeSteps,dt,filenamebase);
        }
    }
}
