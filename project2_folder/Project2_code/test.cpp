#include "test.h"

Test::Test()
{

    //runBandA_Always();
    //runBUpDown(300,2,200,200);
    //runWithoutThermo();
    //parallel_efficiency_3();
    //varyTmeasureP(200,1000, 20);
    //varyTmeasureP(350,350, 1);
    //testLoad();
    run_parallel(4, 20);


}

void Test::runBandA_Always(){
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
    int num_pro = 1;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element, num_pro);
    bSolver->thermostatON = true;
    CellSolver* aSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element, num_pro);
    aSolver->thermostatON = true;

    bSolver->solve(0, numOfTimeSteps, dt, B_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    Andersen = true;
    Berendsen = false;
    aSolver->solve(0, numOfTimeSteps, dt, A_filename, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
}

void Test::runBUpDown(double first, double second, double first_timestep, double second_timestep){
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
    int num_pro = 1;
    int max_time_step = first_timestep+second_timestep;
    CellSolver* bSolver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element, num_pro);
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

void Test::runWithoutThermo(){
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

    double numOfTimeSteps=200;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "testTo1";
    bool Berendsen = false;
    bool Andersen = false;
    double T_bath=0;
    int num_pro = 1;
    CellSolver* solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element, num_pro);
    solver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
}
void Test::varyTmeasureP(double T_start, double T_stop, int numOfPlots){
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
    int max_time_step = numOfPlots*200;
    int numOfTimeStepsNOTHERMO = 130;
    int numOfTimeStepsWITHTHERMO = 70;
    int tot_time = numOfTimeStepsNOTHERMO + numOfTimeStepsWITHTHERMO;
    int num_pro = 1;

    //CellSolver* Solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,max_time_step, b, T, r_cut, element);

    for(int i=0;i<temps.size();i++){
        string s;
        stringstream out;
        out << i;
        s = out.str();
        CellSolver* Solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,tot_time, b, T, r_cut, element, num_pro);
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

double Test::run_parallel(int num_pro, double N){
    double L0 = 3.405;
    double r_cut = 3*L0;
    int Nx = 8, Ny = 8, Nz = 8;
    Nx = N, Ny = N, Nz = N;
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

    double numOfTimeSteps=200;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "testTo1";
    bool Berendsen = false;
    bool Andersen = false;
    double T_bath=0;
    CellSolver* solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element, num_pro);
    solver->print_to_screen = false;
    solver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    return solver->time_spent;
}

double Test::run_serial_code(){
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

    double numOfTimeSteps=200;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "testTo1";
    bool Berendsen = false;
    bool Andersen = false;
    double T_bath=0;
    int num_pro = 1;
    CellSolver* solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Ny, Nz,numOfTimeSteps, b, T, r_cut, element, num_pro);
    solver->print_to_screen = false;
    solver->solve(0, numOfTimeSteps, dt, filename1, writeVMD, writeMeasurements, T_bath, Berendsen, Andersen);
    return solver->time_spent;
}

void Test::parallel_efficiency_1(){
    /*
     *keep the number of processors constant and increase the number of atoms
     *processors = 4
     *N = [8x,10x,12x,14x,16x,18x,20x, 22, 24]
     */
    Col<int> N_list(9);
    N_list<<8<<10<<12<<14<<16<<18<<20<<22<<24;
    vec time_spent_list = zeros(9);
    for(int i=0;i<9;i++){
        cout<<"runs for N = "<<N_list[i]<<endl;
        time_spent_list[i] = run_parallel(1, N_list[i]);


    }
    cout<<time_spent_list<<endl;

}
void Test::parallel_efficiency_2(){
    /*
     *Keep the number of atoms constant and increase the number of processors
     *N = 10
     *P = 1,2,3,4,5,6,7,8,9,10;
     */
    Col<int>P_list(10);
    vec time_spent_list = zeros(10);
    vec speed=zeros(10);
    P_list<<1<<2<<3<<4<<5<<6<<7<<8<<9<<10;
    for(int i=0;i<10;i++){
        cout<<"runs for P= "<<P_list[i]<<endl;
        time_spent_list[i] = run_parallel(P_list[i], 10);

    }
    cout<<time_spent_list<<endl;
    for(int i=0;i<10;i++){
        speed[i] = time_spent_list[0]/time_spent_list[i];
        cout<<P_list[i]<<" "<<speed[i]<<endl;
    }


}
void Test::parallel_efficiency_3(){
    /*
     *Increase the number of processors and the number of atoms so that the number of atoms per processors,N/P, is constant
     */
    Col<int>P_list(5);
    Col<int>N_list(5);
    vec time_spent_list = zeros(5);
    vec time_spent_list_P1 = zeros(5);
    vec speed = zeros(5);
    double ratio = 128;
    N_list<<8<<10<<12<<14<<16;
    for(int i=0;i<5;i++){
        P_list[i] = N_list[i]*N_list[i]*N_list[i]/ratio;
    }

    for(int i=0;i<5;i++){
        cout<<"N= "<<N_list[i]<<" P= "<<P_list[i]<<endl;
        time_spent_list[i] = run_parallel(P_list[i], N_list[i]);
        time_spent_list_P1[i] = run_parallel(1, N_list[i]);
        speed[i] = time_spent_list_P1[i]/time_spent_list[i];
    }
    cout<<speed<<endl;
    cout<<time_spent_list<<endl;

}

void Test::testLoad(){
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

    double numOfTimeSteps=200;
    bool writeVMD = true;
    bool writeMeasurements = true;
    string filename1 = "testLoad.xyz";
    bool Berendsen = false;
    bool Andersen = false;
    double T_bath=0;
    int num_pro = 1;
    int max_time = 200;
    string results_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/";
    string file = "cut_cyllinder.xyz";
    CellSolver* solver = new CellSolver(results_dir+file, CellNx, CellNy, CellNz, Nx,Ny,Nz, max_time,b, T, r_cut,"ar",4);
    solver->writeToVMDfile(results_dir+filename1, "argon", element);

}

