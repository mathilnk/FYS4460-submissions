#include "atom.h"

Atom::Atom(){
}

Atom::Atom(vec position, vec velocity, string element)
{
    /*
      Constructor
      */
    this->position = position;
    this->real_position = position;
    this->initial_position = position;
    this->velocity = velocity;
    this->element = element;
    this->force = zeros(3,1);
    this->canMove = true;


}
