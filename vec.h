#ifndef VEC_H
#define VEC_H

#include <cmath>

class vec // simple 3D vector class
{
  public:
    float x,y,z;

    vec();
    vec(float xx,float yy,float zz);

    vec   operator = (vec t); // assign
    vec   operator + (vec t); // addition
    vec   operator - (vec t); // subtraction
    vec   operator * (float t); // dot product
    float operator * (vec t); // scalar product
    vec   operator ^ (vec t); // cross product

    float length(); // pythagorean length
    vec   unit(); // normalized to a length of 1
    bool  zero(); // returns true if a zero vector
    bool  equals(vec t); // returns true if exactly equal
};

#endif // VEC_H
