#include "lattice.h"





Lattice::Lattice(int Nx, int Ny, int Nz, string element, double b)
{
    this->Nx = Nx;
    this->Ny = Ny;
    this->Nz = Nz;
    this->N = Nx*Ny*Nz;
    this->element = element;
    this->b = b;
    makeLattice();
}

void Lattice::makeLattice(){
    int i,j,k,p,counter;

    arma_lattice = zeros<mat>(4*N, 6);
    FaceCube *unit_cell = new FaceCube(b,element);
    mat unit_cell_table = unit_cell->mat_cube;
    counter = 0;
    for(i=0;i<Nx;i++){
        for(j=0;j<Ny;j++){
            for(k=0;k<Nz;k++){
                for(p=0;p<4;p++){
                    arma_lattice.row(counter++) = unit_cell_table.row(p) + mat("i j k 0 0 0");


                }
            }
        }
    }



}
void Lattice::get_vmd(string comment){
    cout<<4*N<<endl;cout<< comment <<endl;
    int i;
    for(i=0;i<4*N;i++){
        cout<<(arma_lattice.row(i))<<endl;
    }
}
