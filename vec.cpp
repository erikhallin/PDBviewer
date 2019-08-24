#include "vec.h"

vec::vec()
{
    //ctor
}

vec::vec(float xx, float yy, float zz)
{
    x=xx; y=yy; z=zz;
}

vec vec::operator = (vec t) // assign
{
    x=t.x; y=t.y; z=t.z;
    return *this;
}

vec vec::operator + (vec t) // addition
{
    return vec(x+t.x,y+t.y,z+t.z);
}

vec vec::operator - (vec t) // subtraction
{
    return vec(x-t.x,y-t.y,z-t.z);
}

vec vec::operator * (float t) // dot product
{
    return vec(x*t,y*t,z*t);
}

float vec::operator * (vec t) // scalar product
{
    return x*t.x + y*t.y + z*t.z;
}

vec vec::operator ^ (vec t) // cross product
{
    return vec( y*t.z-z*t.y, t.x*z-x*t.z, x*t.y-y*t.x );
}

float vec::length(void) // pythagorean length
{
    return sqrt(x*x + y*y + z*z);
}

vec vec::unit(void) // normalized to a length of 1
{
    float l=length();
    if (l == 0.0) return vec(0.0,0.0,0.0);
    return vec(x/l,y/l,z/l);
}

bool vec::zero(void) // returns true if a zero vector
{
    return x==0 && y==0 && z==0;
}

bool vec::equals(vec t) // returns true if exactly equal
{
    return x==t.x && y==t.y && z==t.z;
}
