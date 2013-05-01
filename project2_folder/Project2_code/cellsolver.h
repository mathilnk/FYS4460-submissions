#ifndef CELLSOLVER_H
#define CELLSOLVER_H

#include<math.h>
#include"cellcontainer.h"
#include"normal.hpp"
#include"omp.h"
#include<libconfig.h++>
#include <armadillo>

class CellSolver
{
public:
    CellSolver(string config_file_name, string experiment_folder, string config_file);
    CellSolver(string filename, int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, int max_time_steps,double b, double T, double r_cut, string element, int num_pro);
    //CellSolver(string filename, int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, int max_time_steps,int start_time_step, double r_cut, double b, int num_pro);
    CellSolver(int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz,int max_time_steps, double b, double T, double r_cut, string element, int num_pro);
    //system parameters
    int CellNx;
    int CellNy;
    int CellNz;
    int CellN;
    vec CellPressures;
    int Nx;
    int Ny;
    int Nz;
    int num_walked_through_wall;
    int N;
    double b;
    double T;
    double mass;
    double k;
    double sig, mean,pi, r_cut, porosity;
    double Lx,Ly,Lz,L0;
    int seed,current_time_step, num_time_steps, real_current_time_step;
    double volume, cellVolume;
    //Conversion factors
    double T0;
    double dt;
    int numOfBins;
    double T_bath;
    bool Andersen, Berendsen;
    bool thermostatON;
    int measure_spatial_pressure_at, num_thermostatON;
    int measurement_frequency;
    string current_time_step_string;
    vector<vec> forces;
    string experiment_name, experiment_folder, result_folder, measurement_filename, next_state_file;
    double Fx;
    int num_movable;
    //initialize
    void initializeContainer();
    void makeLattice();
    void findPosAndMakeAtoms(vec posBase);
    //write to file:
    void writeToVMDfile(string filename, string comment, string element);
    void writeEnergyToFile(string filename);
    void writeMeasurementsToFile(string filename);
    void writeRadialToFile(string filename);

    void loadVMDfile(string filename);
    //Solve_methods
    void solve_one_time_step(double t, double dt, string filename, bool writeVMD,bool writeMeasurements);
    void solve(double t_start, int timesteps, double dt, string filename, bool writeVMD=true,bool writeMeasurements=false, double T_bath = 0, bool Berendsen = false, bool Andersen=false);
    void solve();


    //forces
    void findForces();
    void cleanForces();
    vec force_between(vec r_1, vec r_2, double *pot_energy_thread, double *pressure_sum_thread);
    void measureSpatialPressure();


    //update and measurements
    void updateCells();
    void measureMeanSquareDisplacement();
    double measureTemp();
    double measurePressure();
    double BerendsenThermo();
    void AndersenThermo(Atom *atm);
    void findKinetic();
    void updateKinetic(vec velocity);

    //other
    double gauss(double s, double mean);
    vec findRadial();
    int num_pro;
    bool print_to_screen;
    string writeToFilename;


    CellContainer* myContainer;
    string element;
    int counter;

    double d_energy;

    double pressure_sum;
    double d_pressure;
    double current_displacement;

    //results:
    vec kin_energy;
    vec pot_energy;
    vec temperature;
    vec pressure;
    vec displacement;
    vec radial_distribution;
    double time_spent;

    double count;
    bool slow;

    /////////
    //bools//
    /////////
    bool writeVMD;
    bool writeMeasurements;
};

#endif // CELLSOLVER_H
