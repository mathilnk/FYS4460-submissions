#include "thermostattest.h"

ThermostatTest::ThermostatTest()
{

    //runBandA_Always();
    //runBUpDown(300,2,200,200);
    //runWithoutThermo();
    //varyTmeasureP(200,1000, 20);
    varyTmeasureP(700,700, 1);


}

void ThermostatTest::runBandA_Always(){
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
    double dt = 0.01;

    double numOfTimeSteps = 300;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string B_filename = "B_thermo";
    string A_filename = "A_thermo";
    bool Berendsen = true;
    bool Andersen = false;
    double T_bath = 100;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    bSolver->thermostatON = true;
    CellSolver* aSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    aSolver->thermostatON = true;

    bSolver->solve(0, numOfTimeSteps, dt, B_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    Andersen = true;
    Berendsen = false;
    aSolver->solve(0, numOfTimeSteps, dt, A_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
}

void ThermostatTest::runBUpDown(double first, double second, double first_timestep, double second_timestep){
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
    double dt = 0.01;

    double numOfTimeSteps;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "UpDown";
    //string filename2 = "UpDown2";
    //uble numOfTimeSteps;


    bool Berendsen = true;
    bool Andersen = false;
    double T_bath;
    int max_time_step = first_timestep+second_timestep;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element);
    bSolver->thermostatON = true;
    numOfTimeSteps = first_timestep;
    T_bath = first;
    bSolver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    cout<<"done wiht first"<<endl;
    T_bath = second;
    numOfTimeSteps = second_timestep;
    bSolver->solve(dt*first_timestep, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    cout<<"done with second"<<endl;
}

void ThermostatTest::runWithoutThermo(){
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
    double dt = 0.01;

    double numOfTimeSteps=1;
    bool writeVMD = false;
    bool writeMeasurements = true;
    string filename1 = "testTo";
    bool Berendsen = false;
    bool Andersen = false;
    double T_bath=0;
    CellSolver* solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element);
    solver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
}
void ThermostatTest::varyTmeasureP(double T_start, double T_stop, int numOfPlots){
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
    double dt = 0.01;

    double numOfTimeSteps;
    bool writeVMD = true;
    bool writeMeasurements = false;
    string filename1 = "pressure_temp";
    //string filename2 = "UpDown2";
    //uble numOfTimeSteps;


    bool Berendsen = true;
    bool Andersen = false;
    //double T_bath;
    vec temps = linspace(T_start,T_stop, numOfPlots);

    double T_bath;
    int max_time_step = numOfPlots*400;
    int numOfTimeStepsNOTHERMO = 330;
    int numOfTimeStepsWITHTHERMO = 70;
    int tot_time = numOfTimeStepsNOTHERMO + numOfTimeStepsWITHTHERMO;

    //CellSolver* Solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element);

    for(int i=0;i<temps.size();i++){
        string s;
        stringstream out;
        out << i;
        s = out.str();
        CellSolver* Solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,tot_time, b, T, r_cut, element);
        T_bath = temps[i];
        Solver->thermostatON = true;
        writeMeasurements = false;
        Solver->solve(dt*i*tot_time, numOfTimeStepsWITHTHERMO, dt, filename1+s, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
        Solver->thermostatON = false;
        writeMeasurements = true;
        Solver->solve(dt*i*tot_time + numOfTimeStepsWITHTHERMO, numOfTimeStepsNOTHERMO, dt, filename1+s, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);


    }
    cout<<temps.size()<<endl;
}
