#ifndef MVEC_H
#define MVEC_H

#include <iostream.h>
#include <math.h>
#include "general.hpp"

class  MVec {
public:

  MVec();
  MVec( const MVec& v);
  MVec( double x, double y, double z );

  ~MVec();

  friend istream& operator>>(istream& stream, MVec& v);

  friend ostream& operator<<(ostream& stream, const MVec& v);

  MVec operator-(const MVec& ovec) const;
  MVec operator+(const MVec& ovec) const;
  MVec operator/(double val) const;
  MVec operator*(double val) const;

  const MVec& operator=( const MVec ovec);

  // compare - within epsilon range
  boolean operator==( const MVec ovec) const;

  // cross product
  MVec operator*(const MVec& ovec) const;

  // scalar product
  double operator%(const MVec& ovec) const;
 
  double length() const;
  double length_squared() const;
  double distance_between(MVec & end) const;

  void normalize();

  double x() const;
  double y() const;
  double z() const;

  void  x(double nx);
  void  y(double ny);
  void  z(double nz);

  double operator[](int i) const;

  void move_to_line( const MVec& from, const MVec& to );

private:

  double d_x;
  double d_y;
  double d_z;
  
};

const double c_epsilon = 1e-5;
const double c_epsilon2 = c_epsilon * c_epsilon;
const MVec c_vec_epsilon(c_epsilon,c_epsilon,c_epsilon);   

inline const MVec& MVec::operator=( const MVec ovec){

  d_x = ovec.d_x;
  d_y = ovec.d_y;
  d_z = ovec.d_z;
  return *this;
}

// scalar product
inline double MVec::operator%(const MVec& ovec) const {

  return d_x * ovec.d_x + d_y * ovec.d_y + d_z * ovec.d_z;
}


inline double MVec::length() const {
  return sqrt( d_x * d_x + d_y * d_y + d_z * d_z);
}

inline double MVec::length_squared() const {
  return d_x * d_x + d_y * d_y + d_z * d_z;
}

inline double MVec::distance_between(MVec & end) const {
  return sqrt( pow((d_x - end.x()),2) + 
               pow((d_y - end.y()),2) + 
	       pow((d_z - end.z()),2) ) ;
}

inline void MVec::normalize(){
  double len = length();

  d_x /= len;
  d_y /= len;
  d_z /= len;
}

inline double MVec::x() const { return d_x; }
inline double MVec::y() const { return d_y; }
inline double MVec::z() const { return d_z; }

inline void  MVec::x(double nx) { d_x = nx; }
inline void  MVec::y(double ny) { d_y = ny; }
inline void  MVec::z(double nz) { d_z = nz; }

inline double MVec::operator[](int i) const{
  return ( i == 0 ? d_x : ( i == 1 ? d_y : d_z ) );
}

#endif



