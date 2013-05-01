#include "cellsolver.h"



using namespace std;
using namespace libconfig;
CellSolver::CellSolver(string experiment_name, string experiment_folder, string config_file){
    Config reader;
    this->experiment_name = experiment_name;
    this->experiment_folder = experiment_folder;
    this->result_folder = experiment_folder + "results/";
//string config_file_str = experiment_folder+experiment_name + "_config.cfg";
    const char* config_file_char = config_file.c_str();
    try{
        reader.readFile(config_file_char);
    }
    catch(const FileIOException &fioex){
        std::cerr<<"I/O error while reading file.\nTried reading "<<config_file<<std::endl;
    }
    try{
        ///////////////////
        //read attributes//
        ///////////////////
        int CellN_copy = reader.lookup("CellN");
        this->CellNx = CellN_copy;
        this->CellNy = CellN_copy;
        this->CellNz = CellN_copy;
        this->CellN = CellNx*CellNy*CellNz;
        this->CellPressures = zeros(CellN, 1);
        this->measure_spatial_pressure_at = reader.lookup("measure_spatial_pressure_at");
        this->measurement_frequency = reader.lookup("measurement_frequency");
        this->Fx = reader.lookup("Fx");
        int N_copy = reader.lookup("N");
        this->Nx = N_copy;
        this->Ny = N_copy;
        this->Nz = N_copy;
        this->N = Nx*Ny*Nz;

        T0 = 119.74; //Kelvin
        L0 = 3.405; //AAngstrom
        this->b = reader.lookup("b");
        this->b = this->b/L0;
        this->T = reader.lookup("T_start");

        this->T_bath = reader.lookup("T_bath");
        this->T_bath = this->T_bath/T0;
        //this->T_bath = 1.05;
        this->num_thermostatON = reader.lookup("num_thermostatON");
        this->mass = 1;
        k = 8.617e-5;//ev/K
        this->pi = 3.1415;
        this->seed = 123456754;
        this->Lx = (Nx)*this->b;
        this->Ly = (Ny)*this->b;
        this->Lz = (Nz)*this->b;
        this->r_cut = reader.lookup("r_cut");
        this->r_cut = this->r_cut/L0;
        const char* element_char = reader.lookup("element");
        this->element = element_char;
        this->mean = 0;
        this->sig= sqrt(T/T0);
        this->T = T/T0;
        this->thermostatON = reader.lookup("thermostatON");
        this->writeVMD = reader.lookup("writeVMD");
        this->writeMeasurements = reader.lookup("writeMeasurements");
        const char* writeToFilename_char = reader.lookup("writeToFilename");
        const char* measurement_filename_char = reader.lookup("measurement_file_name");
        this->measurement_filename = measurement_filename_char;
        writeToFilename = writeToFilename_char;
        writeToFilename = this->experiment_name;
        bool makeCrystal = reader.lookup("makeCrystal");


        //string firstStateFile = reader.lookup("firstStateFile");
        string state_file = reader.lookup("state_file");
        const char *next_state_file_char = reader.lookup("next_state_file");
        next_state_file = next_state_file_char;
        //////////////
        //initialize//
        //////////////
        num_walked_through_wall = 0;
        int max_time_steps = reader.lookup("num_time_steps");
        this->temperature = zeros((max_time_steps+1)/measurement_frequency,1);
        this->kin_energy = zeros(max_time_steps+1,1);
        if(makeCrystal){
            temperature[0] = this->T;
            initializeContainer();
            porosity =1;
            this->num_movable = this->myContainer->numberOfAtoms;
        }else{
            cout<<state_file<<endl;
            loadVMDfile(state_file);

        }

//        vec z = zeros(3,1);
//        for(int i=0;i<myContainer->numberOfAtoms;i++){
//            forces.push_back(z);
//        }



        num_time_steps = reader.lookup("num_time_steps");
        dt = reader.lookup("dt");
        Berendsen = reader.lookup("Berendsen");
        Andersen = reader.lookup("Andersen");
        num_pro = reader.lookup("num_pro");
        current_time_step = 0;
        real_current_time_step = reader.lookup("real_current_time_step");
        counter = 0;
        volume = Nx*Ny*Nz*this->b*this->b*this->b;
        cellVolume = volume/CellN;
        d_pressure = 0;
        current_displacement=0;
        numOfBins = CellNx*100;
        radial_distribution = zeros(numOfBins, 1);
        //displacement = 0;

        this->pot_energy = zeros(max_time_steps,1);
        this->pressure = zeros(max_time_steps,1);


        this->displacement = zeros(max_time_steps,1);
        current_time_step_string="0";

        this->num_pro = num_pro;

        this->print_to_screen=true;
        this->time_spent = 0;
        //cout<<"num: "<<myContainer->numberOfAtoms<<endl;

    }
    catch(const SettingNotFoundException &nfex){
        cerr<<"One of the variables not found"<<endl;
    }

}



CellSolver::CellSolver(string filename, int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, int max_time_steps,double b, double T, double r_cut, string element, int num_pro){

    this->CellNx = CellNx;
    this->CellNy = CellNy;
    this->CellNz = CellNz;
    this->CellN = CellNx*CellNy*CellNz;
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;

    T0 = 119.74; //Kelvin
    L0 = 3.405; //AAngstrom
    this->b = b/L0;
    this->T = T;
    this->mass = 1;
    k = 8.617e-5;//ev/K
    this->pi = 3.1415;
    this->seed = 123456754;
    this->Lx = (Nx)*this->b;
    this->Ly = (Ny)*this->b;
    this->Lz = (Nz)*this->b;
    this->r_cut = r_cut/L0;
    this->element = element;
    this->mean = 0;
    this->sig= sqrt(T/T0);
    this->T = T/T0;
    this->thermostatON = false;
    loadVMDfile(filename);
    vec z = zeros(3,1);
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        forces.push_back(z);
    }
    current_time_step = 0;
    counter = 0;
    volume = Nx*Ny*Nz*this->b*this->b*this->b;
    d_pressure = 0;
    current_displacement=0;
    numOfBins = CellNx*100;
    radial_distribution = zeros(numOfBins, 1);
    //displacement = 0;
    this->kin_energy = zeros(max_time_steps,1);
    this->pot_energy = zeros(max_time_steps,1);
    this->pressure = zeros(max_time_steps,1);

    this->temperature = zeros(max_time_steps+1,1);
    this->displacement = zeros(max_time_steps,1);
    this->temperature[0] = this->T;
    current_time_step_string="0";

    this->num_pro = num_pro;

    this->print_to_screen=true;
    this->time_spent = 0;
}


CellSolver::CellSolver(int CellNx, int CellNy, int CellNz, int Nx, int Ny, int Nz, int max_time_steps,double b, double T, double r_cut, string element, int num_pro){
    this->CellNx = CellNx;
    this->CellNy = CellNy;
    this->CellNz = CellNz;
    this->CellN = CellNx*CellNy*CellNz;
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;

    T0 = 119.74; //Kelvin
    L0 = 3.405; //AAngstrom
    this->b = b/L0;
    this->T = T;
    this->mass = 1;
    k = 8.617e-5;//ev/K
    this->pi = 3.1415;
    this->seed = 123456754;
    this->Lx = (Nx)*this->b;
    this->Ly = (Ny)*this->b;
    this->Lz = (Nz)*this->b;
    this->r_cut = r_cut/L0;
    this->element = element;
    this->mean = 0;
    this->sig= sqrt(T/T0);
    this->T = T/T0;
    this->thermostatON = false;

    initializeContainer();

    current_time_step = 0;
    counter = 0;
    volume = Nx*Ny*Nz*this->b*this->b*this->b;
    d_pressure = 0;
    current_displacement=0;
    numOfBins = CellNx*100;
    radial_distribution = zeros(numOfBins, 1);
    //displacement = 0;
    this->kin_energy = zeros(max_time_steps,1);
    this->pot_energy = zeros(max_time_steps,1);
    this->pressure = zeros(max_time_steps,1);

    this->temperature = zeros(max_time_steps+1,1);
    this->displacement = zeros(max_time_steps,1);
    this->temperature[0] = this->T;
    current_time_step_string="0";

    this->num_pro = num_pro;

    this->print_to_screen=true;
    this->time_spent = 0;

}



void CellSolver::solve(){
    /*
      Filename should be without ending. If you want the files to be called test*.xyz, just write filename = test.
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */


    count = 0;
    double start = omp_get_wtime();
    double t_start = 0;
    double time;
    //////////////////////
    //run many timesteps//
    //////////////////////
    //findForces();
    int num_timesteps_with_pressure = 0;
    for(int i=0;i<num_time_steps;i++){
        cout<<"tidssteg: "<<current_time_step<<endl;
        time= t_start + dt*i;
        solve_one_time_step(time, dt, result_folder + writeToFilename, writeVMD, writeMeasurements);
        if(current_time_step<num_thermostatON){
            thermostatON=true;
            cout<<"thermostat is on"<<endl;
        }else{
            cout<<"thermostat is off"<<endl;
            thermostatON=false;
        }
        if(real_current_time_step>measure_spatial_pressure_at){
            measureSpatialPressure();
            num_timesteps_with_pressure++;
        }
        cout<<"number of particles in system: "<<myContainer->numberOfAtoms<<endl;

    }
    //write the spatial pressure to file
    string spatial_pressure_file =  experiment_folder + "/results/" + experiment_name + "_spatial_pressure.txt";
    ofstream pressure_file;
    const char * filename_char = spatial_pressure_file.c_str();
    pressure_file.open(filename_char);
    int help_counter=0;
    for(int j=0;j<CellN;j++){
        pressure_file<<CellPressures[j]/num_timesteps_with_pressure<<" ";
        if(help_counter++==CellNx){
            pressure_file<<"\n";
            help_counter =0;
        }
    }
    cout<<"wrote spatial pressure to "<<spatial_pressure_file<<endl;





    //clock_t stop = clock();
    double stop = omp_get_wtime();
    if(writeMeasurements){
        writeMeasurementsToFile(result_folder + measurement_filename+"_measurements.txt");
        cout<<"wrote measurements to "<<result_folder<<endl;
        //writeRadialToFile(filename+current_time_step_string+"_radial.txt");
    }
    string lastFile = next_state_file;
    myContainer->writeVMDfile(lastFile, "comment", element);
    cout<<"wrote to "<<lastFile<<endl;
    time_spent = stop-start;
    cout<<time_spent<< " s for Cell solver"<<endl;
    int initial_num = N*4;
    double atom_volume = 2*volume/initial_num;
    cout<<"atom_volum in md units:"<<atom_volume<<endl;
    cout<<"number of atoms walked throu wall: "<<num_walked_through_wall<<endl;
    cout<<"U in md units "<<num_walked_through_wall*atom_volume/(dt*num_time_steps)<<endl;

}




void CellSolver::solve(double t_start, int timesteps, double dt, string filename, bool writeVMD, bool writeMeasurements, double T_bath,bool Berendsen, bool Andersen){
    /*
      Filename should be without ending. If you want the files to be called test*.xyz, just write filename = test.
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */
    num_time_steps = timesteps;
    this->dt = dt;
    writeToFilename = filename;
    this->writeVMD = writeVMD;
    this->writeMeasurements = writeMeasurements;
    this->T_bath = T_bath/this->T0;
    this->Andersen = Andersen;
    this->Berendsen = Berendsen;


    solve();
}


void CellSolver::solve_one_time_step(double t, double dt, string filename, bool writeVMD, bool writeMeasurements){
    /*
      Solves for one time step, and returns nothing.
      dt is the time_step_size,
      t is the time
      filename is the base for the filename the lattice is saved to (VMD style).
      if you want the filename to be test*.xyz, write string filename = test.
      */

    /////////////////////////////////////////////////////////////
    //gives the filename an ending, with the current_time_step//
    ///////////////////////////////////////////////////////////
    this->dt = dt;
    //string current_time_step_string;
    stringstream out;
    out<<real_current_time_step;
    current_time_step_string = out.str();
    string filename_end = filename + current_time_step_string + ".xyz";
    //finds the radial distribution function
    //vec g = findRadial();
    ////////cout<<g<<endl;
    //writes to file
    if(writeVMD){
        myContainer->writeVMDfile(filename_end,"comment", element);
    }
//    if(current_time_step>70){
//        findRadial();
//    }
    //writeRadialToFile(filename + current_time_step_string + "_radial.txt");

    ////////cout<<"0 "<<d_energy<<endl;

    ///////////////////////////////////////////////////////////////////////
    //calculates the new velocity and position with the Verlet algorithm//
    /////////////////////////////////////////////////////////////////////

    int n_atoms = myContainer->numberOfAtoms;
    vec v_new(3);
    vec v(3);
    vec r_new(3);
    vec r(3);
    vec dr(3);
    vec v_half(3);
    double dx,dy,dz;
    double Lx = Nx*b;
    double Ly = Ny*b;
    double Lz = Nz*b;
    Cell currentCell;
    //AtomNode* currentNode;
    Atom* atm;
    pressure_sum = 0;



    //////////////////////
    //loop to find r_new//
    //////////////////////

    double maxv = 0;
    double test;
    for(int j=0;j<CellN;j++){
        //Atom* atm = molecule.allAtoms[i];
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atm = currentCell.myAtoms[i];
            if(atm->canMove){
                v = atm->velocity;
                r = atm->position;
                //v_half = v + force_on(r, i)/(2)*dt;
                v_half = v + atm->force/2*dt;
                dr = v_half*dt;
                r_new = r + dr;
                if(r_new[2]<0){
                    num_walked_through_wall-=1;
                }
                else if(r_new[2]>this->Lz){
                    num_walked_through_wall+=1;
                }

                //}
                //else r_new = r;
                ////////cout<<r_new<<endl;
                v_new = v_half; //saves v_half so I dont have to calculate it again in the loop to find v_new
                //cout << v_new << endl;
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Periodic boundary conditions (assume that a particle wont go further than 1000 systems away in the negative direction//
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                r_new[0] = fmod((r_new[0]+1000*Lx), Lx);
                r_new[1] = fmod((r_new[1]+1000*Ly), Ly);
                r_new[2] = fmod((r_new[2]+1000*Lz), Lz);

                //if(atm->canMove){
                atm->position = r_new;
                atm->real_position = atm->real_position + dr;
                atm->velocity = v_new;
                //cout << atm->velocity << endl;
                test = v_new[0]*v_new[0] + v_new[1]*v_new[1] + v_new[2]*v_new[2];
//                if (maxv <test){
//                    maxv = test;
//                    cout << v_new << " " << atm->myindexInAllAtoms << endl;
//                    cout << myContainer->findCellNumberForAtom(atm) << endl;
//                    cout << r_new << r << "ballemor"<< endl;
//                }
            }
            //measureMeanSquareDisplacement(atm->real_position, atm->initial_position);
        }
    }
    //cout << maxv <<" rollfmao " << endl;



    ///////////////////////////////////////////
    //finds the forces with the new positions//
    ///////////////////////////////////////////
    ////////cout<<"before forces"<<endl;
    myContainer->updateCells();
    //cout<<"3"<<endl;
    findForces();
    ////////cout<<"before 2 loop"<<endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //loop to find v_new, now the atoms all have a new position. v_half is saved in atm->velocity//
    ///////////////////////////////////////////////////////////////////////////////////////////////
    double gamma;
    if(thermostatON && Berendsen){
        gamma = BerendsenThermo();
    }else{
        gamma = 1;
        cout<<"thermostat is off"<<endl;
    }

    //////cout<<"T_bath "<<T_bath*T0<<endl;
    for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atm = currentCell.myAtoms[i];
            if(atm->canMove){
                v_half = atm->velocity;
                v_new = v_half + atm->force/2*dt;
                v_new = gamma*v_new;
                if(thermostatON && Andersen){
                    AndersenThermo(atm);
                    cout<<"nope, its on"<<endl;
                }else{
                    atm->velocity = v_new;
                }
//            }else{
//                atm->velocity = zeros(3,1);
//            }
//            if(atm->canMove){
            //d_energy += 0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
            kin_energy[current_time_step+1]+=0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
            }
        }
    }

    double temp = measureTemp();

    measureMeanSquareDisplacement();

    //findRadial();
    if(print_to_screen){
        d_pressure = myContainer->numberOfAtoms*temp/volume + 1/(3*volume)*pressure_sum;
        ////cout<<"pressure sum: "<< 1/(3*volume)*pressure_sum<<endl;
        ///cout<<"pressure constant: "<<myContainer->numberOfAtoms*temp/volume<<endl;
        cout<<"Temperature: "<<temperature[current_time_step+1]*T0<<endl;
        cout<<"Pressure: "<<measurePressure()<<endl;
        cout<<"Kinetic energy: "<<kin_energy[current_time_step+1]<<endl;
        cout<<"Potential energy: "<<pot_energy[current_time_step]<<endl;

        cout<<"displacement: "<<displacement[current_time_step]<<endl;
        cout<<myContainer->numberOfAtoms<<endl;
        cout<<"Porosity: "<<porosity<<endl;
}


    current_time_step++;//update current_time_step;
    real_current_time_step++;

    return;
}


void CellSolver::updateCells(){
    CellContainer* buffer = new CellContainer(CellNx, CellNy, CellNz, r_cut);
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        buffer->addAtom(myContainer->allAtoms[i]);
    }
    myContainer = buffer;
}


void CellSolver::cleanForces(){


//#pragma omp parallel
//    {
//   Cell currentCell;
//   Atom* atm;
//   vec z = zeros(3,1);
//   #pragma omp for
//   for(int j=0;j<CellN;j++){
//        currentCell = myContainer->myCells[j];
//        for(int i=0;i<currentCell.numberOfAtoms;i++){
//            atm = currentCell.myAtoms[i];
//            atm->force = z;
//        }
//    }
//    }
    Atom* atm;
    vec z = zeros(3,1);
    z[2] = this->Fx;

    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm = myContainer->allAtoms[i];
        atm->force = z;
        forces[i] = z;
    }
}


void CellSolver::findForces(){
    cleanForces();
    vector<vec> forces_thread;

    double pot_energy_thread, pressure_sum_thread;
    vec CellPressures_thread;
#pragma omp parallel private(forces_thread, pot_energy_thread, pressure_sum_thread) num_threads(num_pro)
    {
    Atom * atm;
    pot_energy_thread = 0;
    pressure_sum_thread = 0;

    for(int i=0;i<myContainer->numberOfAtoms;i++){
        forces_thread.push_back(zeros(3,1));
    }
    Cell currentCell;
    //AtomNode* currentNode1, *currentNode2;
    Atom *atom_1, *atom_2;
    int index_1,index_2;
    ////////////////////////////////////////////////////////////////
    //first find the forces between the particles in the same cell//
    ////////////////////////////////////////////////////////////////
    //////cout<<"cellN"<<CellN<<endl;
#pragma omp for
    for(int j=0;j<CellN;j++){
        currentCell = myContainer->myCells[j];
        for(int i=0;i<currentCell.numberOfAtoms;i++){
            atom_1 = currentCell.myAtoms[i];
            index_1 = atom_1->myindexInAllAtoms;
            for(int k=i+1;k<currentCell.numberOfAtoms;k++){
                atom_2 = currentCell.myAtoms[k];
                index_2 = atom_2->myindexInAllAtoms;
                if(atom_1->canMove || atom_2->canMove){
                    vec f = force_between(atom_1->position,atom_2->position, &pot_energy_thread, &pressure_sum_thread);
                    forces_thread[index_1] += f;
                    forces_thread[index_2] -=f;
                }
                //cout<<"ting skjer\n"<<k<<endl;
                //////cout<<atom_1->force<<endl;
                //atom_1->force = atom_1->force + f;
                //atom_2->force = atom_2->force - f;

            }
        }
    }



    ///////////////////////////////////////////////////////////////////////
    //Now I have to calculate the forces between atoms in different cells//
    ///////////////////////////////////////////////////////////////////////
    Cell currentCell1, currentCell2;
    vector<int> myNeighbors;
    //////cout<<CellN<<endl;
#pragma omp for
    for(int j=0;j<CellN;j++){
        myNeighbors = myContainer->findMyNeighbors(j);
        //////cout<<myNeighbors.size()<<endl;
        currentCell1 = myContainer->myCells[j];
        for(int i=0;i<myNeighbors.size();i++){
            //////cout<<"crraaaaaaazy"<<endl;
            if(myNeighbors[i]<CellN){
            currentCell2 = myContainer->myCells[myNeighbors[i]];
            for(int k=0;k<currentCell1.numberOfAtoms;k++){
                atom_1 = currentCell1.myAtoms[k];
                index_1 = atom_1->myindexInAllAtoms;
                for(int p=0;p<currentCell2.numberOfAtoms;p++){
                    atom_2 = currentCell2.myAtoms[p];
                    index_2 = atom_2->myindexInAllAtoms;
                    if(atom_1->canMove || atom_2->canMove){
                        vec f = force_between(atom_1->position,atom_2->position,&pot_energy_thread, &pressure_sum_thread);
//                       atom_1->force = atom_1->force + f;
//                      atom_2->force = atom_2->force - f;
                        forces_thread[index_1] +=f;
                        forces_thread[index_2] -=f;
                    }
                }
            }
        }

    }


}
#pragma omp critical
    for(int j =0;j<myContainer->numberOfAtoms;j++){
        forces[j] +=forces_thread[j];
    }
    pressure_sum += pressure_sum_thread;
    pot_energy[current_time_step] += pot_energy_thread;
    }

    Atom *atm;
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm = myContainer->allAtoms[i];
        atm->force = forces[i];
    }



       // cout<<count<<endl;
    for(int i = 0; i<myContainer->numberOfAtoms;i++){
        //ces[i].print("yoyosagga");
        //myContainer->allAtoms[i]->force.print("#yoloswag4lyf");
    }
    //exit(0);
}


void CellSolver::updateKinetic(vec velocity){
    kin_energy[current_time_step] += 0.5*(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
}
void CellSolver::findKinetic(){
    /*should be used before the first timestep. temperature[0] is the initial temperature
     */
    Atom *atm;
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm = myContainer->allAtoms[i];
        updateKinetic(atm->velocity);
    }
}

vec CellSolver::force_between(vec r_1, vec r_2, double *pot_energy_thread, double *pressure_sum_thread){
    /*
      finds the force on r_1 from r_2 and vica versa. the force on r_2 is the same vector, but with different sign.Dimensionless
      */
    vec r = r_1-r_2;
    double Lx = Nx*b;
    double Ly = Ny*b;
    double Lz = Nz*b;
    double new_r_x_min = r[0] - Lx;
    double new_r_y_min = -(Ly - r[1]);
    double new_r_z_min = -(Lz - r[2]);
    double new_r_x_plus = r[0] + Lx;
    double new_r_y_plus = (Ly + r[1]);
    double new_r_z_plus = (Lz + r[2]);




    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Complicated expressions. They choose the components for r(distance between particles) by choosing the smallest of (x_i - x_j + dL) where dL can be (-L,0,L)//
    //Called minimum image convention. (we have periodic boundary conditions.                                                                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    r[0] = (fabs(new_r_x_min)<fabs(new_r_x_plus))*(fabs(new_r_x_min)<fabs(r[0]))*new_r_x_min  +  (fabs(new_r_x_plus)<fabs(new_r_x_min))*(fabs(new_r_x_plus)<fabs(r[0]))*new_r_x_plus    +   (fabs(r[0])<fabs(new_r_x_plus))*(fabs(r[0])<fabs(new_r_x_min))*r[0];
    r[1] = (fabs(new_r_y_min)<fabs(new_r_y_plus))*(fabs(new_r_y_min)<fabs(r[1]))*new_r_y_min  +  (fabs(new_r_y_plus)<fabs(new_r_y_min))*(fabs(new_r_y_plus)<fabs(r[1]))*new_r_y_plus    +   (fabs(r[1])<fabs(new_r_y_plus))*(fabs(r[1])<fabs(new_r_y_min))*r[1];
    r[2] = (fabs(new_r_z_min)<fabs(new_r_z_plus))*(fabs(new_r_z_min)<fabs(r[2]))*new_r_z_min  +  (fabs(new_r_z_plus)<fabs(new_r_z_min))*(fabs(new_r_z_plus)<fabs(r[2]))*new_r_z_plus    +   (fabs(r[2])<fabs(new_r_z_plus))*(fabs(r[2])<fabs(new_r_z_min))*r[2];


    double l = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if(l<0.09){
        cout<<"smalll length: "<<l<<endl;
        r_1.print("r_1");
        r_2.print("r_2");

//cout<<l<<endl;
    }
    double length = l;//std::max(l,0.9);
    //cout<<length<<endl;

    vec f = -24*(1/pow(length,4)-2/pow(length,7))*r;//*length/l;

    if(true){
        *pot_energy_thread += 4*(1/pow(length,6) - 1/pow(length,3));
        *pressure_sum_thread +=(f[0]*r[0] + f[1]*r[1] + f[2]*r[2]);
    }
    //d_energy-=4*(1/pow(length,6) - 1/pow(length,3));
    ////////cout<<4*(1/pow(length,6) - 1/pow(length,3))<<endl;
//    if(current_time_step == 0){
//        vec z = zeros(3,1);
//        return z;
//    }
    return f;
}



void CellSolver::initializeContainer(){
    myContainer = new CellContainer(CellNx, CellNy, CellNz, r_cut);

    makeLattice();

    vec z = zeros(3,1);
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        forces.push_back(z);
    }
}


void CellSolver::makeLattice(){

    vec posBase(3);
    for(int k=0;k<Nz;k++){
        for(int j=0;j<Ny;j++){
            for(int i=0;i<Nx; i++){
                posBase<<i*b<<j*b<<k*b;

                findPosAndMakeAtoms(posBase);


            }
        }
    }

    myContainer->first = false;
}


void CellSolver::writeToVMDfile(string filename, string comment, string element){
    myContainer->writeVMDfile(filename,comment,element);
}


//http://stackoverflow.com/questions/236129/splitting-a-string-in-c by Evan Teran
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


void CellSolver::loadVMDfile(string filename){
    myContainer = new CellContainer(CellNx, CellNy, CellNz, r_cut);
    string experiment = filename;
    this->num_movable=0;
    ifstream myfile (experiment.c_str());
    string line;
    vector<string> contents;
    string name;
    vec3 pos = zeros(3);
    vec3 vel = zeros(3);
    Atom * atm;
    int canMove;;
    int just_not_move;
    just_not_move = 0;
    double kinetic;//, atom_displacement;
    if (myfile.is_open())
    {
        getline(myfile,line);
        int totalAtoms = atoi(line.c_str());
        getline(myfile,line);


        int counter = 0;
        for(int i=0;i<totalAtoms;i++)
        {

            getline (myfile,line);
            contents = split(line, ' ', contents);
            name = contents[0];
            pos(0) = atof(contents[1].c_str());
            pos(1) = atof(contents[2].c_str());
            pos(2) = atof(contents[3].c_str());
            vel(0) = atof(contents[4].c_str());
            vel(1) = atof(contents[5].c_str());
            vel(2) = atof(contents[6].c_str());
            //atom_displacement = atof(contents[7].c_str());

            canMove = atoi(contents[10].c_str());

            atm = new Atom(pos,vel,name);
            atm->canMove = (bool) canMove;
            //atm->myDisplacement = atomt_displacement;

            if(atm->canMove){
                this->num_movable+=1;
                kinetic += 0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));

            }else{
                atm->velocity = zeros(3,1);
            }
            myContainer->addAtom(atm);
            just_not_move +=1-(int)canMove;
            //kinetic += 0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
            this->forces.push_back(zeros(3,1));
            contents.clear();
            counter++;
        }

            myfile.close();
            this->porosity = 1-double(just_not_move)/totalAtoms;
            this->temperature[0] = 2*kinetic/(3*this->num_movable);
            this->kin_energy[0] = kinetic;
            cout<<"temp "<<this->temperature[0]<<endl;
            cout<<"kinetic: "<<kinetic<<endl;
            myfile.close();

            myContainer->updateCells();

            findForces();

                myContainer->first = false;
            }
            else cout << "Unable to open file";




//    myContainer = new CellContainer(CellNx, CellNy, CellNz, r_cut);

//    ifstream myfile;
//    const char* filename_ch = filename.c_str();
//    cout<<"read from "<<filename<<endl;
//    myfile.open(filename_ch);
//    string line, comment, element;
//    string num_of_particles;
//    getline(myfile,num_of_particles);
//    getline(myfile,comment);
//    double x,y,z,vx,vy,vz, fx,fy,fz;
//    int canMove;
//    Atom * atm;
//    int all_num=0;
//    int just_not_move = 0;
//    vec pos(3),vel(3), force;
//    double kinetic = 0;
//    while(getline(myfile,line)){
//        istringstream iss(line);
//        iss>>element>>x>>y>>z>>vx>>vy>>vz>>fx>>fy>>fz>>canMove;
//        pos<<x<<y<<z;
//        vel<<vx<<vy<<vz;
//        //vel = zeros(3,1);
//        force<<fx<<fy<<fz;
//        atm = new Atom(pos,vel,element);

//        if(kinetic ==0){
//  //          cout<<pos<<endl;
////cout<<atm->position<<endl;
//        }
//        atm->canMove = (bool) canMove;
////        if(!canMove){
////            cout<<atm->canMove<<" "<<canMove<<endl;
////        }
//        atm->force = force;
//        this->forces.push_back(atm->force);
//        all_num++;
//        just_not_move +=1-(int)canMove;
//        kinetic += 0.5*(pow(atm->velocity[0],2) + pow(atm->velocity[1],2) + pow(atm->velocity[2],2));
//        this->kin_energy[0] = kinetic;
//        this->temperature[0] =2*kin_energy[0]/(3*myContainer->numberOfAtoms);
//        //cout<<x<<endl;
//        //cout<<atm->canMove<<endl;

//        myContainer->addAtom(atm);
//    }

//    this->porosity = (1-double(just_not_move))/all_num;
//    this->temperature[0] = 2*kinetic/(3*myContainer->numberOfAtoms);
//    myfile.close();
//    cout<<"counted numbers :"<<all_num<<" num in file: "<<num_of_particles<<endl;
//    findForces();

//    myContainer->first = false;
}

void CellSolver::findPosAndMakeAtoms(vec posBase){
    /*
    given a position base x,y,z, this method creates four atoms as described in the exercise.
    The result is a face-centered cubic lattice
    */
    vec posxy(3);
    vec posyz(3);
    vec poszx(3);
    vec velxy(3);
    vec velyz(3);
    vec velzx(3);
    vec velbase(3);

    posxy<<0.5*b<<0.5*b<<0;
    posyz<<0<<0.5*b<<0.5*b;
    poszx<<0.5*b<<0<<0.5*b;
    posxy+=posBase;
    posyz+=posBase;
    poszx+=posBase;

    double s = sig;
    velxy<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velyz<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velzx<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    velbase<<gauss(s, mean)<<gauss(s, mean)<<gauss(s, mean);
    vec velo(3);
    velo<<0<<0<<0;
    Atom *xy = new Atom(posxy,velxy, element);
    Atom *yz = new Atom(posyz,velyz, element);
    Atom *zx = new Atom(poszx,velzx, element);
    Atom *base = new Atom(posBase,velbase, element);

    myContainer->addAtom(xy);

    myContainer->addAtom(yz);

    myContainer->addAtom(zx);

    myContainer->addAtom(base);




}

double CellSolver::measureTemp(){
    double temp = 2*kin_energy[current_time_step+1]/(3*this->num_movable);
    cout<<kin_energy[current_time_step+1]<<endl;
    cout<<myContainer->numberOfAtoms<<endl;
    temperature[current_time_step+1] = temp;
    return temp;
}

double CellSolver::measurePressure(){
    /*
      should rune measureTemp first
      */
    pressure[current_time_step] =  myContainer->numberOfAtoms*temperature[current_time_step]/volume + 1/(3*volume)*pressure_sum;
    cout<<"second part of pressure "<<1/(3*volume)*pressure_sum<<" first part of pressure "<<myContainer->numberOfAtoms*temperature[current_time_step]/volume<<endl;
    pressure_sum = 0;
    return pressure[current_time_step];
}

void CellSolver::measureSpatialPressure(){
    Atom *atm_1, *atm_2;
    Cell currentCell;
    int N;
    double tulle_pot_thread = 0;
    double pressure_sum_tull;
    double cell_kinetic;
    double cell_temp;
    for(int i=0;i<CellN;i++){
        currentCell = myContainer->myCells[i];
        N = currentCell.numberOfAtoms;
        pressure_sum_tull = 0;
        cell_kinetic = 0;
        for(int j=0;j<N;j++){
            atm_1 = currentCell.myAtoms[j];
            double vx = atm_1->velocity[0];
            double vy = atm_1->velocity[1];
            double vz = atm_1->velocity[2];
            cell_kinetic += 0.5*(vx*vx +vy*vy + vz*vz);
            for(int k=j+1;k<N;k++){
                atm_2 = currentCell.myAtoms[k];
                if(atm_1->canMove || atm_2->canMove){
                    vec f = force_between(atm_1->position, atm_2->position,&tulle_pot_thread, &pressure_sum_tull);
                }
            }
        }

        CellPressures[i] += pressure_sum_tull/(3*cellVolume) + cell_kinetic*2/(3*cellVolume);
    }
}

void CellSolver::measureMeanSquareDisplacement(){
    Atom * atm;
    double dis = 0;
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm = myContainer->allAtoms[i];
        dis += (atm->real_position[0] - atm->initial_position[0])*(atm->real_position[0] - atm->initial_position[0]) + (atm->real_position[1] - atm->initial_position[1])*(atm->real_position[1] - atm->initial_position[1]) + (atm->real_position[2] - atm->initial_position[2])*(atm->real_position[2] - atm->initial_position[2]);
        //dis +=atm->myDisplacement;
    }
    dis/=myContainer->numberOfAtoms;
    displacement[current_time_step] = dis;
    //current_displacement = (r_real[0] - r_initial[0])*(r_real[0] - r_initial[0]) + (r_real[1] - r_initial[1])*(r_real[1] - r_initial[1]) + (r_real[2] - r_initial[2])*(r_real[2] - r_initial[2])/myContainer->numberOfAtoms;
    //displacement[current_time_step] += current_displacement;
}
void CellSolver::writeMeasurementsToFile(string filename){
    //cout<<"I am writing"<<endl;
    ofstream file;
    //string newfilename = "/home/mathilde/Dropbox/V2013/FYS4460/FYS4460project2/results/" + filename;
    const char* filename_ch = filename.c_str();
    //////cout<<filename_ch<<endl
    file.open(filename_ch);
    file<<"dt: "<<dt<<" CellNx: "<<CellNx<<" b: "<<b*L0<<" noParticles: "<<myContainer->numberOfAtoms<<endl;
    file<<"Kinetic energy     Potential energy   Total energy   Temperature         Pressure   Average displacement    Radial distribution"<<endl;
    for(int i=0;i<current_time_step;i++){
        file<<kin_energy[i]<<"            "<<pot_energy[i]<<"               "<<kin_energy[i] + pot_energy[i]<<"         "<<temperature[i]*T0<<"         "<<pressure[i]<<"         "<<displacement[i]<<endl;//"           "<<radial_distribution[i]<<endl;
    }
}

void CellSolver::writeEnergyToFile(string filename){
    ofstream file;
    const char* filename_ch = filename.c_str();
    //////cout<<filename_ch<<endl;
    file.open(filename_ch);
    for(int i=0;i<current_time_step;i++){
        file<<kin_energy[i]<<"      "<<pot_energy[i]<<"     "<<kin_energy[i] + pot_energy[i]<<endl;
    }

}

void CellSolver::writeRadialToFile(string filename){
    ofstream file;
    const char * filename_ch = filename.c_str();
    file.open(filename_ch);
    for(int i=0;i<numOfBins;i++){
        file<<radial_distribution[i]<<endl;
    }
}

vec CellSolver::findRadial(){
    /*
      finds the radial distribution g(r)
      */
    //double numOfBins = CellNx*10;
    //vec bins = zeros(numOfBins, 1);
    Atom * atm1;
    Atom * atm2;
    double distance;
    double maxlength = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    int index;
    double binSize = maxlength/numOfBins;
    int maxIndex=0;
    vec r;
    for(int i=0;i<myContainer->numberOfAtoms;i++){
        atm1 = myContainer->allAtoms[i];
        for(int j=i;j<myContainer->numberOfAtoms;j++){
            atm2 = myContainer->allAtoms[j];
            r = atm1->position - atm2->position;

            double new_r_x_min = r[0] - Lx;
            double new_r_y_min = -(Ly - r[1]);
            double new_r_z_min = -(Lz - r[2]);
            double new_r_x_plus = r[0] + Lx;
            double new_r_y_plus = (Ly + r[1]);
            double new_r_z_plus = (Lz + r[2]);




            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //Complicated expressions. They choose the components for r(distance between particles) by choosing the smallest of (x_i - x_j + dL) where dL can be (-L,0,L)//
            //Called minimum image convention. (we have periodic boundary conditions.                                                                                    //
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            r[0] = (fabs(new_r_x_min)<fabs(new_r_x_plus))*(fabs(new_r_x_min)<fabs(r[0]))*new_r_x_min  +  (fabs(new_r_x_plus)<fabs(new_r_x_min))*(fabs(new_r_x_plus)<fabs(r[0]))*new_r_x_plus    +   (fabs(r[0])<fabs(new_r_x_plus))*(fabs(r[0])<fabs(new_r_x_min))*r[0];
            r[1] = (fabs(new_r_y_min)<fabs(new_r_y_plus))*(fabs(new_r_y_min)<fabs(r[1]))*new_r_y_min  +  (fabs(new_r_y_plus)<fabs(new_r_y_min))*(fabs(new_r_y_plus)<fabs(r[1]))*new_r_y_plus    +   (fabs(r[1])<fabs(new_r_y_plus))*(fabs(r[1])<fabs(new_r_y_min))*r[1];
            r[2] = (fabs(new_r_z_min)<fabs(new_r_z_plus))*(fabs(new_r_z_min)<fabs(r[2]))*new_r_z_min  +  (fabs(new_r_z_plus)<fabs(new_r_z_min))*(fabs(new_r_z_plus)<fabs(r[2]))*new_r_z_plus    +   (fabs(r[2])<fabs(new_r_z_plus))*(fabs(r[2])<fabs(new_r_z_min))*r[2];

            distance = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            index = distance/binSize;
            ////////cout<<index<<" "<<binSize<<" "<<numOfBins<<endl;
            if(index>maxIndex){
                maxIndex = index;
            }
            double v = 4*pi*(index*binSize)*(index*binSize)*binSize;
            ////cout<<binSize*L0<<endl;
            radial_distribution[index]+=1/v;

        }
    }
    ////////cout<<"max: "<<maxIndex<<endl;
    return radial_distribution;
}

double CellSolver::gauss(double s, double mean){
    /*
      finds a random number from a gaussian distribution. the seed changes when the method is called
      */

    //return fmod(rand(),(4*s)) - 2*s;
    return r8_normal(mean,s,seed);
}

double CellSolver::BerendsenThermo(){
    double tau = 10*dt;
    double gamma = sqrt(1 + dt/tau*(this->T_bath/this->temperature[current_time_step] -1));
    //cout << "gamma er " << gamma << endl;
    return gamma;

}

void CellSolver::AndersenThermo(Atom * atm){
    double rand_num = (double) rand()/RAND_MAX;
    ////////cout<<rand_num<<endl;
    double tau = dt*10;
    if(rand_num<dt/tau){
        double s = sqrt(this->T_bath);
        vec new_velo(3);
        new_velo<<gauss(s,0)<<gauss(s,0)<<gauss(s,0);
        atm->velocity = new_velo;

    }
}
