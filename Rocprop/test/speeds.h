/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: speeds.h,v 1.6 2008/12/06 08:45:28 mtcampbe Exp $

#ifndef _SPEEDS_H_
#define _SPEEDS_H_

#include "propbasic.h"

PROP_BEGIN_NAMESPACE

class Speed {
public:
  Speed() {}
  virtual ~Speed() {}

  // Obtain the velocity (vector) at a given point.
  // The appropriate functions must be defined by the sub-class, or a
  // runtime error (assertion failure) will occur.
  virtual Vector_3 get_velocity( const Point_3 &p, double t) const
  { assert(false); return Vector_3(0,0,0); }
};

// Define translation speed.
class Translate_speed : public Speed {
public:
  explicit Translate_speed( double s, bool b=false) : _spd(s) {}
  virtual ~Translate_speed() {}

  Vector_3 get_velocity( const Point_3 &, double t) const 
  { return Vector_3(_spd,0,0); }

protected:
  const double _spd;
};

// Define rotation speed with z=0.
class Rotate_speed : public Speed {
public:
  explicit Rotate_speed() {}
  virtual ~Rotate_speed() {}

  Vector_3 get_velocity( const Point_3 &pnt, double t) const { 

    return Vector_3(-pnt[1], pnt[0], 0);
  }
};

class Vortex_flow : public Speed {
public:
  // T is the reversal time
  explicit Vortex_flow( double T=2) : _t(T>0?T:2) {}
  virtual ~Vortex_flow() {}
  
  Vector_3 get_velocity( const Point_3 &p, double t) const {
    const double pi = 3.14159265358979;
    
    Vector_3 v;
    v[0] = square(std::sin(pi*p[0]))*(std::sin(2*pi*p[2])-std::sin(2*pi*p[1]));
    v[1] = square(std::sin(pi*p[1]))*(std::sin(2*pi*p[0])-std::sin(2*pi*p[2]));
    v[2] = square(std::sin(pi*p[2]))*(std::sin(2*pi*p[1])-std::sin(2*pi*p[0]));

    double s;
    if ( _t<0) s = std::cos(pi*(0.5-t/_t));
    else s = std::cos(pi*t/_t);

    return (v *= s);
  }

protected:
  double square( double x) const { return x*x; }
  double      _t;
};

class LeVeque_flow : public Speed {
public:
  // T is the reversal time
  explicit LeVeque_flow( double T=3) : _t(T) {}
  virtual ~LeVeque_flow() {}
  
  Vector_3 get_velocity( const Point_3 &p, double t) const {
    const double pi = 3.14159265358979;
    
    Vector_3 v;
    v[0] = 2*square(std::sin(pi*p[0]))*std::sin(2*pi*p[1])*std::sin(2*pi*p[2]);
    v[1] = -std::sin(2*pi*p[0])*square(std::sin(pi*p[1]))*std::sin(2*pi*p[2]);
    v[2] = -std::sin(2*pi*p[0])*std::sin(2*pi*p[1])*square(std::sin(pi*p[2]));

    double s;
    if ( _t<0) s = std::cos(pi*(0.5-t/_t));
    else s = std::cos(pi*t/_t);

    return (v *= s);
  }

protected:
  double square( double x) const { return x*x; }
  double      _t;
};


PROP_END_NAMESPACE

#endif // _SPPEDS_H_






