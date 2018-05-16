
#include "MVec.hpp"
#include <math.h>

MVec::MVec() {}
MVec::MVec( const MVec& v) : 
  d_x(v.d_x),
  d_y(v.d_y),
  d_z(v.d_z) {}

MVec::MVec( double x, double y, double z ) :
  d_x(x),
  d_y(y),
  d_z(z) {}

MVec::~MVec() {}


istream& operator>>(istream& stream, MVec& v) {
  stream >> v.d_x >> v.d_y >> v.d_z;
  return stream;
}


ostream& operator<<(ostream& stream, const MVec& v) {
  stream << ' ' << v.d_x << ' ' << v.d_y << ' ' << v.d_z;
  return stream;
}

MVec MVec::operator-(const MVec& ovec) const{

  MVec v;
  v.d_x = d_x - ovec.d_x;
  v.d_y = d_y - ovec.d_y;
  v.d_z = d_z - ovec.d_z;
  return v;
}

MVec MVec::operator+(const MVec& ovec) const {

  MVec v;
  v.d_x = d_x + ovec.d_x;
  v.d_y = d_y + ovec.d_y;
  v.d_z = d_z + ovec.d_z;
  return v;
}

MVec MVec::operator/(double val) const{

  MVec v;
  v.d_x = d_x/val;
  v.d_y = d_y/val;
  v.d_z = d_z/val;
  return v;
}

MVec MVec::operator*(double val) const{

  MVec v;
  v.d_x = d_x*val;
  v.d_y = d_y*val;
  v.d_z = d_z*val;
  return v;
}

// cross product
MVec MVec::operator*(const MVec& ovec) const {

  MVec v;
  v.d_x = d_y * ovec.d_z - d_z * ovec.d_y;
  v.d_y = d_z * ovec.d_x - d_x * ovec.d_z;
  v.d_z = d_x * ovec.d_y - d_y * ovec.d_x;
  return v;
}

boolean MVec::operator==( const MVec ovec) const {

  if( ovec.d_x > d_x + c_epsilon 
      || ovec.d_x < d_x - c_epsilon
      || ovec.d_y > d_y + c_epsilon 
      || ovec.d_y < d_y - c_epsilon
      || ovec.d_z > d_z + c_epsilon 
      || ovec.d_z < d_z - c_epsilon ) {
    return FALSE;
  }
  return ( (*this - ovec).length_squared() < c_epsilon2 ? TRUE : FALSE );
}


void MVec::move_to_line( const MVec& from, const MVec& to ){

  MVec base = to - from;
  MVec vec = *this - from;
  double b2 = base.length_squared();
  double dot = base%vec;
  if (dot <= 0.0) {
    *this = from;
    return;
  }
  double percent2 = ((dot*dot)/b2)/b2;
  if (percent2 >= 1.0){
    *this = to;
    return;
  }
  double delta = sqrt (percent2);
  *this = from * ( 1.0 - delta ) + to * delta;
}
