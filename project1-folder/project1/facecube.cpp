#include "facecube.h"

FaceCube::FaceCube(double b, string element)
{
    this->b = b;
    this->element = element;
    makeCube();
}


void FaceCube:: makeCube()
{

    vec pos00(3);
    vec posxy(3);
    vec posyz(3);
    vec posxz(3);

    vec vecxy(3);
    vec vecyz(3);
    vec vecxz(3);
    vec vec00(3);

    pos00 <<0<< 0<< 0;
    posxy <<0.5*b<< 0.5*b<<0;
    posyz <<0<< 0.5*b<< 0.5*b;
    posxz <<0.5*b<< 0<< 0.5*b;

    vecxy <<0<<0<<0;
    vecyz<<0<<0<<0;
    vecxz<<0<<0<<0;
    vec00<<0<<0<<0;
    //mat_cube = mat("0 0 0 0 0 0;0.5*b 0.5*b 0 0 0 0;0 0.5*b 0.5*b 0 0 0;0.5*b 0 0.5*b 0 0 0;");
    mat_cube = b*mat("0 0 0 0 0 0;0.5 0.5 0 0 0 0;0 0.5 0.5 0 0 0;0.5 0 0.5 0 0 0;");


    Atom* origo = new Atom(pos00, vec00, element);
    Atom* xy = new Atom(posxy, vecxy, element);
    Atom* yz = new Atom(posyz, vecyz, element);
    Atom* zx = new Atom(posxz, vecxz,element);



}
