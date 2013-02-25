#include "verlet_solver.h"
using namespace std;

Verlet_solver::Verlet_solver(Lattice &molecule){
    /*
      Constructor
      */

    this->molecule = molecule;
    current_time_step = 0;
}

void Verlet_solver::solve(double t_start,int timesteps, double dt, string filename){
    /*
      Filename should be without ending. If you want the files to be called test*.xyz, just write filename = test.
      This method solves for many timesteps, and for each timestep a xyz file is made.
      */

    double time = t_start;
    clock_t start = clock();
    ///////////////////////////////////////////////////
    //before the first timestep: calculate the forces//
    ///////////////////////////////////////////////////
    findForces();
    for(int i=0;i<timesteps;i++){
        time = t_start + dt*i;
        this->molecule = solve_one_time_step(time, dt, filename);
    }
    clock_t stop = clock();
    cout<<(stop-start)/CLOCKS_PER_SEC<< " s for Verlet Solver"<<endl;

}



Lattice Verlet_solver::solve_one_time_step(double t, double dt, string filename){
    /*
      Solves for one time step, and returns the Lattice at the next timestep.
      dt is the time_step_size,
      t is the time
      filename is the base for the filename the lattice is saved to (VMD style).
      if you want the filename to be test*.xyz, write string filename = test.
      */

    /////////////////////////////////////////////////////////////
    //gives the filename an ending, with the current_time_step//
    ///////////////////////////////////////////////////////////
    //cout<<"tidssteg: "<<current_time_step<<endl;
    string current_time_step_string;
    stringstream out;
    out<<current_time_step;
    current_time_step_string = out.str();
    string filename_end = filename + current_time_step_string + ".xyz";
    //writes to file
    molecule.writeVMDfile(filename_end,"comment");

    ///////////////////////////////////////////////////////////////////////
    //calculates the new velocity and position with the Verlet algorithm//
    /////////////////////////////////////////////////////////////////////
    int n_atoms = molecule.numberOfAtoms;
    vec v_new(3);
    vec v(3);
    vec r_new(3);
    vec r(3);
    vec dr(3);
    vec v_half(3);
    double dx,dy,dz;
    double Lx = molecule.Lx;
    double Ly = molecule.Ly;
    double Lz = molecule.Lz;
    double m = molecule.mass;

    //////////////////////
    //loop to find r_new//
    //////////////////////
    for(int i=0;i<n_atoms;i++){
        Atom* atm = molecule.allAtoms[i];
        v = atm->velocity;
        r = atm->position;
        //v_half = v + force_on(r, i)/(2)*dt;
        v_half = v + atm->force/2*dt;
        r_new = r + v_half*dt;
        v_new = v_half; //saves v_half so I dont have to calculate it again in the loop to find v_new

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //Periodic boundary conditions (assume that a particle wont go further than 1000 systems away in the negative direction//
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        r_new[0] = fmod((r_new[0]+1000*Lx), Lx);
        r_new[1] = fmod((r_new[1]+1000*Ly), Ly);
        r_new[2] = fmod((r_new[2]+1000*Lz), Lz);



        atm->position = r_new;
        atm->velocity = v_new;
    }
    ///////////////////////////////////////////
    //finds the forces with the new positions//
    ///////////////////////////////////////////
    findForces();
    ///////////////////////////////////////////////////////////////////////////////////////////////
    //loop to find v_new, now the atoms all have a new position. v_half is saved in atm->velocity//
    ///////////////////////////////////////////////////////////////////////////////////////////////
    for(int i=0;i<n_atoms;i++){
        Atom * atm = molecule.allAtoms[i];
        v_half = atm->velocity;
        //v_new = v_half + force_on(atm->position, i)/(2)*dt;
        v_new = v_half + atm->force/2*dt;
        atm->velocity = v_new;
    }

    current_time_step++;//update current_time_step;
    return molecule;
}


vec Verlet_solver::force_on(vec r, int index){
    /*
      finds the force on one particle from all the other atoms in the lattice
      r is the particle, and index tells us which place in allAtoms r lies.
      */
    vec f = zeros(3,1);
    int n_atoms = molecule.numberOfAtoms;
    for(int i=0;i<n_atoms;i++){
        if(i!=index){
            Atom * atm_other = molecule.allAtoms[i];
            f = f + force_between(r, atm_other->position);

        }

    }

    return f;
}



vec Verlet_solver::force_between(vec r_1, vec r_2){
    /*
      finds the force on r_1 from r_2 and vica versa. the force on r_2 is the same vector, but with different sign.Dimensionless
      */
    vec r = r_1-r_2;
    double Lx = molecule.Lx;
    double Ly = molecule.Ly;
    double Lz = molecule.Lz;
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


    double length = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    length = std::max(0.9, length);
    vec f = -24*(1/pow(length,8)-2/pow(length,14))*r;
    return f;
}

void Verlet_solver::cleanForces(){
    vec z = zeros(3,1);
    Atom * atm;
    for(int i=0; i<molecule.numberOfAtoms;i++){
        atm = molecule.allAtoms[i];
        atm->force = z;
    }
}

 void Verlet_solver::findForces(){
    cleanForces();
    int n_atoms = molecule.numberOfAtoms;
    vec f;
    Atom * atom_1;
    Atom * atom_2;
    for(int i=0;i<n_atoms;i++){
        atom_1 = molecule.allAtoms[i];
        for(int j=i+1;j<n_atoms;j++){
            atom_2 = molecule.allAtoms[j];
            f = force_between(atom_1->position, atom_2->position);
            atom_1->force = atom_1->force + f;
            atom_2->force = atom_2->force - f;

        }
    }

}
