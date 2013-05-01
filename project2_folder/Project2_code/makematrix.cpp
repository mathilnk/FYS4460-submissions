#include "makematrix.h"

MakeMatrix::MakeMatrix()
{
//    makePorousMatrix(20,30,);
//    mySolver->myContainer->writeOnlyNOTMovAtomsVMD("/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/one_pore1.xyz", "one pore", "ar");
//    mySolver->myContainer->writeOnlyMovAtomsVMD("/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/one_pore2.xyz", "one pore", "ar");
}



CellSolver* MakeMatrix::runCyllinder(){
    double b = 5.72;
    double T = 0.851*119.74;
    vec centrum = zeros(2,1);
    double R = 20;
    //mySolver = makeAndThermalize(20, b, T);
    string result_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/";
    string startfilename = "thermalized_argon_N20_b5_72_T_0_851.xyz";
    mySolver = makeFromFile(result_dir+startfilename, 20, b, 500);
    double Lx = mySolver->Lx/2;
    centrum<<Lx<<Lx;
    cout<<"Lx: "<<Lx<<endl;
    cout<<"mySolver Lx: "<<mySolver->Lx<<endl;
    cutCyllinder(centrum, R);
    string results_dir = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/";
    string filename = "cut_cyllinder2.xyz";
    //mySolver->myContainer->writeVMDfile(results_dir + filename, "cyllinder", "ar");
    mySolver->myContainer->writeOnlyNOTMovAtomsVMD(results_dir + filename, "arg", "ar");
}

CellSolver* MakeMatrix::makeFromFile(string filename, int N, double b, int max_time){
    double L0 = 3.405;
    double r_cut = 3*L0;
    int CellNx, CellNy, CellNz;
    //double b = 5.26;
    double T = 100;
    int Nx = N;
    CellNx = Nx*b/r_cut;


    CellNy = CellNx;
    CellNz = CellNx;
    r_cut = Nx*b/(CellNx);
    cout<<r_cut<<" "<<CellNy<<endl;
    string element = "Ar";
//    double dt = 0.01;
//    int num_pro = 4;
//    int tot_time = 1000;
//    int thermalize_time = 100;


    bool writeVMD = true;
    bool writeMeasurements =true;
    //string filename1 = "noe";

    bool Berendsen = true;
    bool Andersen = false;
    int start_time_step = 0;
    int max_time_steps = max_time;
    CellSolver* Solver = new CellSolver(filename,CellNx, CellNy, CellNz, Nx,Nx,Nx,max_time_steps, b, T,r_cut,"ar",4);
    return Solver;
}

CellSolver* MakeMatrix::makeAndThermalize(int N,double b, double T){
    double L0 = 3.405;
    double r_cut = 3*L0;
    int CellNx, CellNy, CellNz;
    //double b = 5.26;
    //double T = 100;
    int Nx = N;
    CellNx = Nx*b/r_cut;


    CellNy = CellNx;
    CellNz = CellNx;
    r_cut = Nx*b/(CellNx);
    string element = "Ar";
    double dt = 0.01;
    int num_pro = 4;
    int tot_time = 1000;
    int thermalize_time = 100;


    bool writeVMD = true;
    bool writeMeasurements =true;
    string filename1 = "makeCyllinderMatrix";

    bool Berendsen = true;
    bool Andersen = false;
    CellSolver* Solver = new CellSolver(CellNy, CellNx, CellNz, Nx,Nx, Nx,tot_time, b, T, r_cut, element, num_pro);
    Solver->thermostatON = true;
    Solver->solve(0, thermalize_time, dt, filename1, writeVMD, writeMeasurements, T, Berendsen, Andersen);
    Solver->thermostatON = false;
    return Solver;
}
void MakeMatrix::cutCyllinder(vec centrum, double R){
    double L0 = 3.405;
    double MD_R = R/mySolver->L0;
    vec MD_centrum = centrum;//two entries, x and y
    MD_centrum[0] = MD_centrum[0];
    MD_centrum[1] = MD_centrum[1];
    Atom * atm;
    CellContainer *container = mySolver->myContainer;
    double x,y,radius_sq;
    for(int i=0;i<container->numberOfAtoms;i++){
        atm = container->allAtoms[i];
        x = atm->position[0];
        y = atm->position[1];

        radius_sq = (x-MD_centrum[0])*(x-MD_centrum[0]) + (y-MD_centrum[1])*(y-MD_centrum[1]);
        //cout<<"radius_sq "<<radius_sq<<" R: "<<R*R<<endl;
        if(radius_sq<MD_R*MD_R){
            atm->canMove=true;

        }else{
            atm->canMove = false;
            atm->velocity = zeros(3,1);
        }
    }
    container->updateCells();
}

vec random_vec(double min, double max, int N){
    vec r_vec = zeros(N,1);
    //srand(time(NULL));
    double r;
    for(int i=0;i<N;i++){
        r = (double)rand()/RAND_MAX;
        r_vec[i] = min + r*(max-min);
    }
    return r_vec;
}
double random_double(double min, double max){
    //srand(time(NULL));
    double r = (double) rand()/RAND_MAX;;
    return min + r*(max-min);
}

void MakeMatrix::makePorousMatrix(double r_max, double r_min, int max_time){
    string filename = "thermalized_argon_N20_b5_72_T_0_851.xyz";
    string path = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/";
    mySolver = makeFromFile(path+filename,20,5.72, max_time);
    CellContainer * con = mySolver->myContainer;
    double L0 = mySolver->L0;
    double MD_r_max = r_max/L0;
    double MD_r_min = r_min/L0;
    double average_atom_volume = mySolver->volume/con->numberOfAtoms;
    int remove_count = 0;
    srand(time(NULL));
    for(int j=0;j<100;j++){
        vec pos = random_vec(0, mySolver->Lx, 3);

        double random_r = random_double(MD_r_min, MD_r_max);

        double x,y,z, radius_sq;
        Atom * atm;
        for(int i=0;i<con->numberOfAtoms;i++){
            atm = con->allAtoms[i];
            x = atm->position[0];
            y = atm->position[1];
            z = atm->position[2];
            radius_sq = (x-pos[0])*(x-pos[0]) + (y-pos[1])*(y-pos[1])+ (z - pos[2])*(z-pos[2]);
            if(radius_sq<random_r){
                atm->canMove = false;
                remove_count+=1;
            }
        }
    }
    con->updateCells();
    double pore_volume = average_atom_volume*remove_count;
    double porosity = 1-pore_volume/mySolver->volume;
    cout<<"Porosity of the system: "<<porosity<<endl;
}
