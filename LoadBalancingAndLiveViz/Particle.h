#ifndef PARTICLE_H
#define PARTICLE_H

/*
*Particle object with x&y coordinate components
*/


enum EColor
{
  e_RED,
  e_BLUE,
  e_GREEN
};

#include <vector>
class Particle  {
public:
  double x;
  double y;
  int color;
  Particle() { }
  Particle(double a, double b,EColor color_) { 
    x=a; y=b; color = color_; 
  }

  void pup(PUP::er &p){
    p|x;
    p|y;
    p|color; //need to test for correctness
  }

};

typedef std::vector<Particle>::iterator PartIter;

#endif
