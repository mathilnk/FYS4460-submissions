#include "lattice.h"




Lattice::Lattice(){}
Lattice::Lattice(int Nx, int Ny, int Nz, string element, double b, double T, double mass)
{
    /*
      Constructor. Makes the lattice
      */
    double T0 = 119.74;//K
    double L0 = 3.405;//AAngstrom
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;
    this->numberOfAtoms = N*4;
    this->element = element;
    this->b = b/L0;
    this->T = T;

    //this->k =1.38e-23; //J/K
    k = 8.617e-5;//ev/K
    //this->mass = mass*931e6;//ev/c^2
    //this->s = sqrt(k*T/mass)*1e-5;//aangstrom per femtosekund
    //this->s = 144.0;
    this->mass = 1;
    this->s = sqrt(T/T0);
    cout<<s<<endl;
    this->mean = 0;
    this->pi = 3.1415;
    this->seed = 123456789;
    this->Lx = (Nx)*this->b;
    this->Ly = (Ny)*this->b;
    this->Lz = (Nz)*this->b;


    makeLattice();
}

void Lattice::makeLattice(){
    /*
      Private
      */
    vec posBase(3);
    for(int i=0;i<Nx;i++){
        for(int j=0;j<Ny;j++){
            for(int k=0;k<Nz; k++){
                posBase<<i*b<<j*b<<k*b;
                findPosAndMakeAtoms(posBase);
            }
        }
    }
}








void Lattice::writeVMDfile(string filename, string comment){
    /*
      This method writes a VMD-formatted file called filename.
      */
    ofstream writeToFile;
    numberOfAtoms = allAtoms.size();
    const char* filename_ch = filename.c_str();
    writeToFile.open(filename_ch);
    writeToFile<<numberOfAtoms<<endl;
    writeToFile<<comment<<endl;


    ///////////////////////////////////////////////////////////////////////////////////
    //runs through the atoms in allAtoms, and write the position and velocity to file//
    ///////////////////////////////////////////////////////////////////////////////////
    for(int i=0;i<numberOfAtoms;i++){
        Atom* a = allAtoms[i];
        writeToFile<<element<<" ";
        for(int j=0;j<3;j++){
            //AAngstrom
            writeToFile<<a->position[j]<<" ";
        }
        for(int k=0;k<3;k++){
             //aangstrom per femtosekund
            writeToFile<<a->velocity[k]<<" ";
        }
        writeToFile<<"\n";
    }


}

void Lattice::findPosAndMakeAtoms(vec posBase){
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
    allAtoms.push_back(xy);
    allAtoms.push_back(yz);
    allAtoms.push_back(zx);
    allAtoms.push_back(base);

}


double Lattice::gauss(double s, double mean){
    /*
      finds a random number from a gaussian distribution. the seed changes when the method is called
      */

    return r8_normal(mean,s,seed);
}



