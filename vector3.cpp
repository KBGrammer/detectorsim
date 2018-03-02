#include "vector3.h"
#include <math.h>

Vector3::Vector3()
{
  
}

Vector3::Vector3(double x1, double y1, double z1)
{
  x = x1; y = y1; z = z1;
}

Vector3::~Vector3()
{

}

void Vector3::set_coords(double x1, double y1, double z1)
{
  x = x1; y = y1; z = z1;
}

void Vector3::return_coords(double & x1, double & y1, double & z1)
{
  x1 = x; y1 = y; z1 = z;
}

// void Vector3::make_unitvector()
// {
//   double mag = x*x+y*y+z*z;
//   x = x/mag;
//   y = y/mag;
//   z = z/mag;
// }

double Vector3::dot(const Vector3 & vec)
{
  return x*vec.x+y*vec.y+z*vec.z;
}

Vector3& Vector3::operator=(const Vector3& rhs)
{
  return *this;
}
