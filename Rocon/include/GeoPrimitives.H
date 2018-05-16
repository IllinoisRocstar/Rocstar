/// 
/// \file
/// \ingroup support
/// \brief Geometric Primitives interface
///
#ifndef _GeoPrimitives_H_
#define _GeoPrimitives_H_
#include <cmath>
#include <istream>
#include <ostream>
#include <vector>
#include <cassert>


namespace GeoPrim {

  const double TOL = 1e-16;
  const double PI  = 4.0*atan(1.0);

  class CLine;

  class CPoint {
    friend CPoint operator*(double,const CPoint &p);
    friend std::ostream &operator<<(std::ostream &oSt,const CPoint &p);
    friend std::istream &operator>>(std::istream &iSt,CPoint &p);
    friend double Distance(const CPoint &p,const GeoPrim::CLine &l);
  protected:
    double P[3];
  public:
    CPoint()
    {
      P[0] = P[1] = P[2] = 0.0;
    };
    CPoint(double ix,double iy,double iz = 0.0)
    {
      P[0] = ix;
      P[1] = iy;
      P[2] = iz;
    };
    CPoint(const double *p)
    {
      P[0] = p[0];
      P[1] = p[1];
      P[2] = p[2];
    };
    CPoint(const CPoint &in)
    {
      P[0] = in.P[0];
      P[1] = in.P[1];
      P[2] = in.P[2];
    };
    CPoint &init()
    {
      P[0] = 0.;
      P[1] = 0.;
      P[2] = 0;
      return(*this);
    };
    CPoint &init(double ix,double iy,double iz = 0.0)
    {
      P[0] = ix;
      P[1] = iy;
      P[2] = iz;
      return(*this);
    };
    CPoint &init(const double *p)
    {
      P[0] = p[0];
      P[1] = p[1];
      P[2] = p[2];
      return(*this);
    };
    CPoint &operator=(const CPoint &p)
    {
      P[0] = p.P[0];
      P[1] = p.P[1];
      P[2] = p.P[2];
      return(*this);
    };
    CPoint &copy(const CPoint &p)
    {
      *this = p;
      return(*this);
    };
    double norm() const
    {
      return(sqrt((x()*x()) + (y()*y()) + (z()*z())));
    };
    bool operator==(const CPoint &p) const
    {
      return( (fabs(P[0] - p.P[0]) < TOL) &&
	      (fabs(P[1] - p.P[1]) < TOL) &&
	      (fabs(P[2] - p.P[2]) < TOL) );
    };
    bool operator!=(const CPoint &p) const
    {
      return (!((*this)==p));
    };
    bool operator<(const CPoint &p) const
    {
      return (norm() < p.norm());
    };
    bool operator>(const CPoint &p) const
    {
      return (norm() > p.norm());
    };
    bool operator<=(const CPoint &p) const
    {
      return ((norm() < p.norm()) ||
	      (*this == p));
    };
    bool operator>=(const CPoint &p) const
    {
      return((norm() > p.norm()) ||
	     (*this == p));
    };
    CPoint &operator+=(const CPoint &p)
    {
      P[0] += p.P[0];
      P[1] += p.P[1];
      P[2] += p.P[2];
      return(*this);
    };
    CPoint &operator-=(const CPoint &p)
    {
      P[0] -= p.P[0];
      P[1] -= p.P[1];
      P[2] -= p.P[2];
      return(*this);
    };
    CPoint operator*(double num) const
    {
      CPoint rp(P[0]*num,P[1]*num,P[2]*num);
      return(rp);
    };
    CPoint &operator*=(double num)
    {
      P[0] *= num;
      P[1] *= num;
      P[2] *= num;
      return(*this);
    };
    double &x(){ return P[0]; };
    CPoint &x(double i) 
    {
      P[0] = i;
      return(*this);
    };
    const double &x() const { return P[0]; };
    double &y(){ return P[1]; };
    CPoint &y(double i)
    {
      P[1] = i;
      return(*this);
    };
    const double &y() const { return P[1]; };
    double &z(){ return P[2]; };
    CPoint &z(double i)
    {
      P[2] = i;
      return(*this);
    }
    const double &z() const { return P[2]; };
    const double &operator[](unsigned int i) const
    {
      if( i > 2 ) i = 0;
      return(P[i]);
    };
    double &operator[](unsigned int i)
    {
      if (i > 2) i = 0;
      return(P[i]);
    };
    CPoint operator+(const CPoint &p) const
    {
      CPoint rp(P[0]+p.P[0],P[1]+p.P[1],P[2]+p.P[2]);
      return(rp);
    };
    CPoint operator-(const CPoint &p) const
    {
      CPoint rp(P[0]-p.P[0],P[1]-p.P[1],P[2]-p.P[2]);
      return(rp);
    };
    double Distance(const CPoint &p) const
    {
      return((*this - p).norm());
    }
    CPoint &set(double ix,double iy,double iz)
    {
      P[0] = ix;
      P[1] = iy;
      P[2] = iz;
      return(*this);
    };
    CPoint &reset()
    {
      P[0] = 0.0;
      P[1] = 0.0;
      P[2] = 0.0;
      return(*this);
    };
  };
  
  class CVector {
    friend CVector operator*(double,const CVector &);
    friend std::ostream &operator<<(std::ostream &,const CVector &);
    friend std::istream &operator>>(std::istream &,CVector &); 
  protected:
    double V[3];
  public:
    CVector()
    {
      V[0] = V[1] = V[2] = 0.;
    };
    CVector(double ix,double iy,double iz = 0.0)
    {
      V[0] = ix;
      V[1] = iy;
      V[2] = iz;
    };
    CVector(const CPoint &p1,const CPoint &p2)
    {
      V[0] = (p2 - p1).x();
      V[1] = (p2 - p1).y();
      V[2] = (p2 - p1).z();
    };
    CVector(const CPoint &p)
    {
      V[0] = p.x();
      V[1] = p.y();
      V[2] = p.z();
    };
    CVector(const double *a)
    {
      V[0] = a[0];
      V[1] = a[1];
      V[2] = a[2];
    };
    CVector(const CVector &v)
    {
      V[0] = v.V[0];
      V[1] = v.V[1];
      V[2] = v.V[2];
    };
    CVector &init(const CPoint &p)
    {
      V[0] = p[0];
      V[1] = p[1];
      V[2] = p[2];
      return(*this);
    };
    CVector &init(const double *a)
    {
      V[0] = a[0];
      V[1] = a[1];
      V[2] = a[2];
      return(*this);
    };
    CVector &init(const CVector &v)
    {
      V[0] = v[0];
      V[1] = v[1];
      V[2] = v[2];
      return(*this);
    };
    CVector &init(double ix,double iy,double iz)
    {
      V[0] = ix;
      V[1] = iy;
      V[2] = iz;
      return(*this);
    };
    CVector &init(const CPoint &p1,const CPoint &p2)
    {
      V[0] = (p2 - p1).x();
      V[1] = (p2 - p1).y();
      V[2] = (p2 - p1).z();
      return(*this);
    };
    CVector &init()
    {
      V[0] = V[1] = V[2] = 0.0;
      return(*this);
    };
    double norm() const
    {
      return(sqrt((V[0] * V[0]) + (V[1] * V[1]) + (V[2] * V[2])));
    };
    double mag() const
    {
      return (norm());
    };
    double mag2() const
    {
      return ((V[0]*V[0])+(V[1]*V[1])+(V[2]*V[2]));
    };
    CVector &normalize()
    {
      double n = norm();
      V[0] /= n;
      V[1] /= n;
      V[2] /= n;
      return(*this);
    };
    CVector unit() const
    {
      CVector v(*this);
      double n = norm();
      v.V[0] /= n;
      v.V[1] /= n;
      v.V[2] /= n;
      return(v);
    };
    double operator*(const CVector &v2) const
    {
      return((V[0] * v2.V[0]) +
	     (V[1] * v2.V[1]) +
	     (V[2] * v2.V[2]));
    };
    CVector operator*(double num) const
    {
      CVector rv(num*V[0],num*V[1],num*V[2]);
      return(rv);
    };
    CVector &operator*=(double num) 
    {
      V[0] *= num;
      V[1] *= num;
      V[2] *= num;
      return(*this);
    };
    CVector operator%(const CVector &v2) const
    {
      CVector v;
      v.V[0] = ((v2.V[2] * V[1]) - (v2.V[1] * V[2]));
      v.V[1] = ((v2.V[0] * V[2]) - (v2.V[2] * V[0]));
      v.V[2] = ((v2.V[1] * V[0]) - (v2.V[0] * V[1]));
      return(v);
    };
    CVector &copy(const CVector &v2)
    {
      V[0] = v2.V[0];
      V[1] = v2.V[1];
      V[2] = v2.V[2];
      return(*this);
    };
    CVector &operator=(const CVector &v2)
    {
      return(copy(v2));
    };
    CVector &operator=(const CPoint &p)
    {
      V[0] = p.x();
      V[1] = p.y();
      V[2] = p.z();
      return(*this);
    };
    CVector &operator%=(const CVector &v2)
    {
      return(copy((*this)%v2));
    };
    CVector operator+(const CVector &v2) const
    {
      CVector v(V[0] + v2.V[0], V[1] + v2.V[1], V[2] + v2.V[2]);
      return(v);
    };
    CVector operator-(const CVector &v2) const
    {
      return((*this) + (-1.0 * v2));
    };
    CVector &operator+=(const CVector &v2)
    {
      V[0] += v2.V[0];
      V[1] += v2.V[1];
      V[2] += v2.V[2];
      return(*this);
    };
    CVector &operator+=(const CPoint &p)
    {
      V[0] += p.x();
      V[1] += p.y();
      V[2] += p.z();
      return(*this);
    };
    CVector &operator-=(const CVector &v2)
    {
      return(copy((*this) - v2));
    };
    double &x() 
    {
      return(V[0]);
    };
    double x() const
    {
      return(V[0]);
    };
    CVector &x(double i)
    {
      V[0] = i;
      return(*this);
    };
    double &y()
    {
      return(V[1]);
    };
    double y() const
    {
      return(V[1]);
    };
    CVector &y(double i)
    {
      V[1] = i;
      return(*this);
    };
    double &z()
    {
      return(V[2]);
    };
    double z() const
    {
      return(V[2]);
    };
    CVector &z(double i)
    {
      V[2] = i;
      return(*this);
    };
    const double &operator[](unsigned int i) const
    {
      if(i > 2) i = 0;
      return(V[i]);
    };
    double &operator[](unsigned int i)
    {
      if(i > 2) i = 0;
      return(V[i]);
    };
    bool operator==(const CVector &v) const
    {
      return( (fabs(V[0] - v.V[0]) < TOL) &&
	      (fabs(V[1] - v.V[1]) < TOL) &&
	      (fabs(V[2] - v.V[2]) < TOL) );
    };
    bool operator<(const CVector &v) const
    {
      return (norm() < v.norm());
    };
    bool operator>(const CVector &v) const
    {
      return (norm() > v.norm());
    };
    bool operator<=(const CVector &v) const
    {
      return ( (*this < v) ||
	       (*this == v));
    };
    bool operator>=(const CVector &v) const
    {
      return ( (*this > v) ||
	       (*this == v) );
    };
  };
  
  
  ///
  /// /brief Cartesian 3 Vector
  ///
  class C3Vector {
    friend C3Vector operator*(double,const C3Vector &);
    friend std::ostream &operator<<(std::ostream &,const C3Vector &);
    friend std::istream &operator>>(std::istream &,C3Vector &);
  private:
    bool _mydata;
  protected:
    double *V;
  public:
    C3Vector()
    {
      _create();
      V[0] = V[1] = V[2] = 0.;
    };
    C3Vector(const GeoPrim::CVector &inv)
    {
      _create();
      V[0] = inv.x();
      V[1] = inv.y();
      V[2] = inv.z();      
    };
    C3Vector(const CPoint &inp)
    {
      _create();
      V[0] = inp.x();
      V[1] = inp.y();
      V[2] = inp.z();
    };
    C3Vector(double ix,double iy,double iz = 0.0)
    {
      _create();
      V[0] = ix;
      V[1] = iy;
      V[2] = iz;
    };
    C3Vector(const C3Vector &p1,const C3Vector &p2)
    {
      _create();
      Copy(p2-p1);
    };
    C3Vector(double *a)
    {
      _mydata = false;
      V = a;
    };
    C3Vector(const double *a)
    {
      _create();
      V[0] = a[0];
      V[1] = a[1];
      V[2] = a[2];
    };
    C3Vector(const C3Vector &v)
    {
      _create();
      Copy(v);
    };
    ~C3Vector()
    {
      _destroy();
    };
    C3Vector &Copy(const C3Vector &v)
    {
      V[0] = v.V[0];
      V[1] = v.V[1];
      V[2] = v.V[2];
      return(*this);
    };
    C3Vector &Init(const double *a)
    {
      V[0] = a[0];
      V[1] = a[1];
      V[2] = a[2];
      return(*this);
    };
    C3Vector &Init(const C3Vector &v)
    {
      V[0] = v[0];
      V[1] = v[1];
      V[2] = v[2];
      return(*this);
    };
    C3Vector &Init(double ix,double iy,double iz)
    {
      V[0] = ix;
      V[1] = iy;
      V[2] = iz;
      return(*this);
    };
    C3Vector &Init(const C3Vector &p1,const C3Vector &p2)
    {
      Copy(p2 - p1);
      return(*this);
    };
    C3Vector &Init()
    {
      V[0] = V[1] = V[2] = 0.0;
      return(*this);
    };
    C3Vector &Clear()
    {
      V[0] = V[1] = V[2] = 0.0;
      return(*this);
    };
    double Norm() const
    {
      return(sqrt((V[0] * V[0]) + (V[1] * V[1]) + (V[2] * V[2])));
    };
    double Mag() const
    {
      return (Norm());
    };
    double Mag2() const
    {
      return ((V[0]*V[0])+(V[1]*V[1])+(V[2]*V[2]));
    };
    C3Vector &Normalize()
    {
      double n = Norm();
      V[0] /= n;
      V[1] /= n;
      V[2] /= n;
      return(*this);
    };
    C3Vector Unit() const
    {
      C3Vector v(*this);
      return(v.Normalize());
    };
    ///
    /// \brief Dot Product
    ///
    double operator*(const C3Vector &v2) const
    {
      return((V[0] * v2.V[0]) +
	     (V[1] * v2.V[1]) +
	     (V[2] * v2.V[2]));
    };
    ///
    /// \brief Scaling
    ///
    C3Vector operator*(double num) const
    {
      C3Vector rv(num*V[0],num*V[1],num*V[2]);
      return(rv);
    };
    C3Vector &operator*=(double num) 
    {
      V[0] *= num;
      V[1] *= num;
      V[2] *= num;
      return(*this);
    };
    ///
    /// \brief Cross Product
    ///
    C3Vector operator%(const C3Vector &v2) const
    {
      C3Vector v;
      v.V[0] = ((v2.V[2] * V[1]) - (v2.V[1] * V[2]));
      v.V[1] = ((v2.V[0] * V[2]) - (v2.V[2] * V[0]));
      v.V[2] = ((v2.V[1] * V[0]) - (v2.V[0] * V[1]));
      return(v);
    };
    C3Vector &operator=(const C3Vector &v2)
    {
      return(Copy(v2));
    };
    C3Vector &operator%=(const C3Vector &v2)
    {
      double t0 = V[0];
      double t1 = V[1];
      double t2 = V[2];
      V[0] = ((v2.V[2] * t1) - (v2.V[1] * t2));
      V[1] = ((v2.V[0] * t2) - (v2.V[2] * t0));
      V[2] = ((v2.V[1] * t0) - (v2.V[0] * t1));
      return(*this);
    };
    C3Vector operator+(const C3Vector &v2) const
    {
      C3Vector v(V[0] + v2.V[0], V[1] + v2.V[1], V[2] + v2.V[2]);
      return(v);
    };
    C3Vector operator-(const C3Vector &v2) const
    {
      return((*this) + (-1.0 * v2));
    };
    C3Vector &operator+=(const C3Vector &v2)
    {
      V[0] += v2.V[0];
      V[1] += v2.V[1];
      V[2] += v2.V[2];
      return(*this);
    };
    C3Vector &operator-=(const C3Vector &v2)
    {
      V[0] -= v2.V[0];
      V[1] -= v2.V[1];
      V[2] -= v2.V[2];
      return(*this);
    };
    double &x() 
    {
      return(V[0]);
    };
    double x() const
    {
      return(V[0]);
    };
    C3Vector &x(double i)
    {
      V[0] = i;
      return(*this);
    };
    double &y()
    {
      return(V[1]);
    };
    double y() const
    {
      return(V[1]);
    };
    C3Vector &y(double i)
    {
      V[1] = i;
      return(*this);
    };
    double &z()
    {
      return(V[2]);
    };
    double z() const
    {
      return(V[2]);
    };
    C3Vector &z(double i)
    {
      V[2] = i;
      return(*this);
    };
    const double &operator[](unsigned int i) const
    {
      assert(i>0 && i<4);
      return(V[i-1]);
    };
    double &operator[](unsigned int i)
    {
      assert(i>0 && i<4);
      return(V[i-1]);
    };
    bool operator==(const C3Vector &v) const
    {
      return( (fabs(V[0] - v.V[0]) < TOL) &&
	      (fabs(V[1] - v.V[1]) < TOL) &&
	      (fabs(V[2] - v.V[2]) < TOL) );
    };
    bool operator<(const C3Vector &v) const
    {
      return (Mag2() < v.Mag2());
    };
    bool operator>(const C3Vector &v) const
    {
      return (Mag2() > v.Mag2());
    };
    bool operator<=(const C3Vector &v) const
    {
      return ( (*this < v) ||
	       (*this == v));
    };
    bool operator>=(const C3Vector &v) const
    {
      return ( (*this > v) ||
	       (*this == v) );
    };
    void DestroyData() { _destroy(); };
    void SetData(double *d) { 
      V = d; 
      _mydata = false;
    };
  private:
    void _create()
    {
      V = new double [3];
      _mydata = true;
    };
    void _destroy()
    {
      if(_mydata){
	delete [] V;
	_mydata = false;
	V = NULL;
      }
    };
  };
  
  typedef C3Vector C3Point;

  
//   class C3Point {
//     friend C3Point operator*(double,const C3Point &p);
//     friend std::ostream &operator<<(std::ostream &oSt,const C3Point &p);
//     friend std::istream &operator>>(std::istream &iSt,C3Point &p);
//     friend double Distance(const C3Point &p,const CLine &l);
//   private:
//     bool mydata;
//   protected:
//     double *P;
//     //    double P[3];
//   public:
//     C3Point()
//     {
//       _create();
//       P[0] = P[1] = P[2] = 0.0;
//     };
//     C3Point(double ix,double iy,double iz = 0.0)
//     {
//       _create();
//       P[0] = ix;
//       P[1] = iy;
//       P[2] = iz;
//     };
//     C3Point(double *p)
//     {
//       mydata = false;
//       P = p;
//       P[0] = p[0];
//       P[1] = p[1];
//       P[2] = p[2];
//     };
//     C3Point(const double *p)
//     {
//       _create();
//       P[0] = p[0];
//       P[1] = p[1];
//       P[2] = p[2];
//     };
//     C3Point(C3Point &in)
//     {
//       mydata = false;
//       P = in.P;
//     };
//     C3Point(const C3Point &in)
//     {
//       Copy(in);
//     };
//     ~C3Point()
//     {
//       _destroy();
//     };
//     bool Good()
//     {
//       return(!(P==NULL));
//     };
//     void _destroy()
//     {
//       if(mydata){
// 	mydata = false;
// 	delete [] P;
// 	P = NULL;
//       }
//     };
//     void _create()
//     {
//       if(!mydata){
// 	mydata = true;
// 	P = new double [3];
//       }
//       Clear();
//     };
//     void Copy(const C3Point &in)
//     {
//       //      _create();
//       P[0] = in.P[0];
//       P[1] = in.P[1];
//       P[2] = in.P[2];
//     };
//     C3Point &Clear()
//     {
//       P[0] = 0.;
//       P[1] = 0.;
//       P[2] = 0;
//       return(*this);
//     };
//     C3Point &Create()
//     {
//       _create();
//       return(*this);
//     };
//     C3Point &Init(double ix,double iy,double iz = 0.0)
//     {
//       P[0] = ix;
//       P[1] = iy;
//       P[2] = iz;
//       return(*this);
//     };
//     C3Point &Init(const double *p)
//     {
//       P[0] = p[0];
//       P[1] = p[1];
//       P[2] = p[2];
//       return(*this);
//     };
//     C3Point &operator=(const C3Point &p)
//     {
//       P[0] = p.P[0];
//       P[1] = p.P[1];
//       P[2] = p.P[2];
//       return(*this);
//     };
//     //    CPoint &copy(const CPoint &p)
//     //    {
//     //      *this = p;
//     //      return(*this);
//     //    };
//     double Norm() const
//     {
//       return(sqrt((x()*x()) + (y()*y()) + (z()*z())));
//     };
//     bool operator==(const C3Point &p) const
//     {
//       return( (fabs(P[0] - p.P[0]) < TOL) &&
// 	      (fabs(P[1] - p.P[1]) < TOL) &&
// 	      (fabs(P[2] - p.P[2]) < TOL) );
//     };
//     bool operator!=(const C3Point &p) const
//     {
//       return (!((*this)==p));
//     };
//     bool operator<(const C3Point &p) const
//     {
//       return (Norm() < p.Norm());
//     };
//     bool operator>(const C3Point &p) const
//     {
//       return (Norm() > p.Norm());
//     };
//     bool operator<=(const C3Point &p) const
//     {
//       return ((Norm() < p.Norm()) ||
// 	      (*this == p));
//     };
//     bool operator>=(const C3Point &p) const
//     {
//       return((Norm() > p.Norm()) ||
// 	     (*this == p));
//     };
//     C3Point &operator+=(const C3Point &p)
//     {
//       P[0] += p.P[0];
//       P[1] += p.P[1];
//       P[2] += p.P[2];
//       return(*this);
//     };
//     C3Point &operator-=(const C3Point &p)
//     {
//       P[0] -= p.P[0];
//       P[1] -= p.P[1];
//       P[2] -= p.P[2];
//       return(*this);
//     };
//     C3Point operator*(double num) const
//     {
//       C3Point rp(P[0]*num,P[1]*num,P[2]*num);
//       return(rp);
//     };
//     C3Point &operator*=(double num)
//     {
//       P[0] *= num;
//       P[1] *= num;
//       P[2] *= num;
//       return(*this);
//     };
//     double &x(){ return P[0]; };
//     C3Point &x(double i) 
//     {
//       P[0] = i;
//       return(*this);
//     };
//     const double &x() const { return P[0]; };
//     double &y(){ return P[1]; };
//     C3Point &y(double i)
//     {
//       P[1] = i;
//       return(*this);
//     };
//     const double &y() const { return P[1]; };
//     double &z(){ return P[2]; };
//     C3Point &z(double i)
//     {
//       P[2] = i;
//       return(*this);
//     }
//     const double &z() const { return P[2]; };
//     const double &operator[](unsigned int i) const
//     {
//       assert(i > 0 && i < 4);
//       //      if( i > 2 ) i = 0;
//       return(P[i-1]);
//     };
//     double &operator[](unsigned int i)
//     {
//       assert(i > 0 && i < 4);
//       return(P[i-1]);
//     };
//     C3Point operator+(const C3Point &p) const
//     {
//       C3Point rp(P[0]+p.P[0],P[1]+p.P[1],P[2]+p.P[2]);
//       return(rp);
//     };
//     C3Point operator-(const C3Point &p) const
//     {
//       C3Point rp(P[0]-p.P[0],P[1]-p.P[1],P[2]-p.P[2]);
//       return(rp);
//     };
//     double Distance(const C3Point &p) const
//     {
//       return((*this - p).Norm());
//     }
//     C3Point &Set(double ix,double iy,double iz)
//     {
//       P[0] = ix;
//       P[1] = iy;
//       P[2] = iz;
//       return(*this);
//     };
//     C3Point &Reset()
//     {
//       P[0] = 0.0;
//       P[1] = 0.0;
//       P[2] = 0.0;
//       return(*this);
//     };
//   };
  
  class CLine {
    friend double Distance(const CPoint &p,const CLine &l);
    friend double Distance(const C3Point &p,const CLine &l);
  protected:
    CPoint  p;
    CVector v;
  public:
    CLine()
    {};
    CLine(const CPoint &p1,const CPoint &p2)
      : p(p1),v(p1,p2)
    {};
    CLine(const CLine &l)
      : p(l.p), v(l.v)
    {};
    CLine(const CPoint &p1,const CVector &v1)
      : p(p1),v(v1)
    {};
    CLine(const CVector &v1)
      : p(0,0,0),v(v1)
    {};
    CLine &move(const CPoint &p1)
    {
      p = p1;
      return(*this);
    };
    bool operator==(const CLine &l) const
    {
      return( (p == l.p)  &&
	      (v == l.v) );
    };
    bool has_point(const CPoint &p2) const
    {
      CVector v2(p2,p);
      v2.normalize();
      double r = v2 * v.unit();
      return( (fabs(r - 1.0) < TOL) || 
	      (fabs(r + 1.0) < TOL) );
    };
    CVector &slope()
    {
      return(v);
    };
    const CVector &slope() const
    {
      return(v);
    };
    CPoint &point()
    {
      return(p);
    };
    const CPoint &point() const
    {
      return(p);
    };
  };
  
  class CLineSegment {
  protected:
    CPoint p1;
    CPoint p2;
  public:
    CLineSegment()
      : p1(0,0,0),p2(0,0,0)
    {};
    CLineSegment(const CPoint &p)
      : p1(0,0,0),p2(p)
    {};
    CLineSegment(const CPoint &pb,const CPoint &pe)
      : p1(pb),p2(pe)
    {};
    double length() const
    {
      CVector v(p1,p2);
      return(v.norm());
    };
    CVector slope() const
    {
      CVector v(p1,p2);
    return(v);
    };
    bool has_point(const CPoint &p3) const
    {
      if(p3 == p1 || p3 == p2) return true;
      if(p1 == p2) return false;
      CLine l(p1,p2);
      if(!l.has_point(p3)) return false;
      double u;
      if(p1.x() != p2.x())
	u = (p3.x() - p1.x())/(p2.x() - p1.x());
      else if(p1.y() != p2.y())
	u = (p3.y() - p1.y())/(p2.y() - p1.y());
      else if(p1.z() != p2.z())
	u = (p3.z() - p1.z())/(p2.z() - p1.z());
      return ( u >= 0.0 && u <= 1.0 );
    };
    CPoint &point1()
    {
      return (p1);
    };
    CPoint &point2()
    {
      return (p2);
    };
    const CPoint &point1() const
    {
      return (p1);
    };
    const CPoint &point2() const
    {
      return (p2);
    };
    CLineSegment &point1(const CPoint &p)
    {
      p1 = p;
      return(*this);
    };
    CLineSegment &point2(const CPoint &p)
    {
      p2 = p;
      return(*this);
    };
  };
  
  class CBox {
    friend std::ostream &operator<<(std::ostream &,const CBox &);
  private:
    bool initd;
  protected:
    CPoint p1; // xmin,ymin,zmin
    CPoint p2; // xmax,ymax,zmax
  public:
    CBox()
    {
      p1.reset();
      p2.reset();
      initd=false;
    };
    CBox(const std::vector<CPoint> &pv)
    {
      std::vector<CPoint>::const_iterator pvi = pv.begin();
      p1.x(pvi->x());
      p2.x(p1.x());
      p1.y(pvi->y());
      p2.y(p1.y());
      p1.z(pvi->z());
      p2.z(p1.z());
      pvi++;
      while(pvi != pv.end()){
	if(pvi->x() < p1.x()) p1.x(pvi->x());
	if(pvi->x() > p2.x()) p2.x(pvi->x());
	if(pvi->y() < p1.y()) p1.y(pvi->y());
	if(pvi->y() > p2.y()) p2.y(pvi->y());
	if(pvi->z() < p1.z()) p1.z(pvi->z());
	if(pvi->z() > p2.z()) p2.z(pvi->z());
	pvi++;
      }
      initd = true;
    }
    CBox(const double *points,unsigned int n)
    {
      p1.x(points[0]);
      p1.y(points[1]);
      p1.z(points[2]);
      p2 = p1;
      for(unsigned int i = 1;i<n;i++){
	if(points[i*3] < p1.x()) p1.x(points[i*3]);
	if(points[i*3] > p2.x()) p2.x(points[i*3]);
	if(points[(i*3)+1] < p1.y()) p1.y(points[(i*3)+1]);
	if(points[(i*3)+1] > p2.y()) p2.y(points[(i*3)+1]);
	if(points[(i*3)+2] < p1.z()) p1.z(points[(i*3)+2]);
	if(points[(i*3)+2] > p2.z()) p2.z(points[(i*3)+2]);
      }
      initd = true;
    }; 
    CBox(const CPoint &ip1,const CPoint &ip2)
      : p1(ip1),p2(ip2)
    {initd = true;};
    CBox(const CBox &box)
      : p1(box.p1),p2(box.p2)
    {initd = true;};
    const CPoint &P1() const {return p1;};
    const CPoint &P2() const {return p2;};
    bool empty() const {return ((std::fabs(p2.x()-p1.x()) < TOL) ||
				(std::fabs(p2.y()-p1.y()) < TOL) ||
				(std::fabs(p2.z()-p1.z()) < TOL));};
    CPoint &P1() {return p1;};
    CPoint &P2() {return p2;};
    void init(const double *points,unsigned int n)
    {
      initd = true;
      p1.x(points[0]);
      p1.y(points[1]);
      p1.z(points[2]);
      p2 = p1;
      for(unsigned int i = 1;i<n;i++){
	if(points[i*3] < p1.x()) p1.x(points[i*3]);
	if(points[i*3] > p2.x()) p2.x(points[i*3]);
	if(points[(i*3)+1] < p1.y()) p1.y(points[(i*3)+1]);
	if(points[(i*3)+1] > p2.y()) p2.y(points[(i*3)+1]);
	if(points[(i*3)+2] < p1.z()) p1.z(points[(i*3)+2]);
	if(points[(i*3)+2] > p2.z()) p2.z(points[(i*3)+2]);
      }
    }; 
    CBox &operator=(const CBox &b)
    {
      p1 = b.p1;
      p2 = b.p2;
      initd = true;
      return(*this);
    };
    CBox around(const CPoint &p) const
    {
      CVector v(p1,p2);
      CPoint minp(p.x()-fabs(v.x()/2.0),
		p.y()-fabs(v.y()/2.0),
		  p.z()-fabs(v.z()/2.0));
      CPoint maxp(p.x()+fabs(v.x()/2.0),
		  p.y()+fabs(v.y()/2.0),
		  p.z()+fabs(v.z()/2.0));
      CBox box(minp,maxp);
      return(box);
    };
    bool operator==(const CBox &b) const
    {
      return(p1 == b.p1 && p2 == b.p2);
    };
    //    CBox operator&(const CBox &b) const
    //    {
    //      return(this->intersect(b));
    //    };
    CPoint &min()
    {
    return(p1);
    };
    const CPoint &min() const
    {
      return(p1);
    };
    CPoint &max()
    {
      return(p2);
    };
    const CPoint &max() const 
    {
      return(p2);
    };
    bool contains(const CPoint &p) const
    {
    return( p.x() <= p2.x() && p.x() >= p1.x() &&
	    p.y() <= p2.y() && p.y() >= p1.y() &&
	    p.z() <= p2.z() && p.z() >= p1.z() );
    };
    bool contains(const double *p) const
    {
      return( p[0] <= p2.x() && p[0] >= p1.x() &&
	      p[1] <= p2.y() && p[1] >= p1.y() &&
	      p[2] <= p2.z() && p[2] >= p1.z() );
    };
    // Volume of the box
    double size() const 
    {
      CVector v(p1,p2);
      return(fabs(v.x()*v.y()*v.z()));
    };
    bool operator<(const CBox &b) const
    {
      return(this->size() < b.size());
    };
    bool operator>(const CBox &b) const
    {
      return(this->size() > b.size());
    };    
    void AddPoint(const CPoint &p)
    {
      if(!initd){
	p1 = p;
	p2 = p;
	initd = true;
      }
      else{
	if(p.x() < p1.x()) p1.x() = p.x();
	if(p.x() > p2.x()) p2.x() = p.x();
	if(p.y() < p1.y()) p1.y() = p.y();
	if(p.y() > p2.y()) p2.y() = p.y();
	if(p.z() < p1.z()) p1.z() = p.z();
	if(p.z() > p2.z()) p2.z() = p.z();
      }
    };
    void merge(const CBox &b)
    {
      if(initd){
	CPoint p(p1);
	p1.x((b.p1.x() > p1.x() ? p.x() : b.p1.x()));
	p1.y((b.p1.y() > p1.y() ? p.y() : b.p1.y()));
	p1.z((b.p1.z() > p1.z() ? p.z() : b.p1.z()));
	p = p2;
	p2.x(b.p2.x() < p2.x() ? p.x() : b.p2.x());
	p2.y(b.p2.y() < p2.y() ? p.y() : b.p2.y());
	p2.z(b.p2.z() < p2.z() ? p.z() : b.p2.z());
      }
      else {
	p1 = b.p1;
	p2 = b.p2;
	initd = true;
      }
    };
    bool collide(const CBox &b) const
    {
      if(this->empty() || b.empty())
	return(false);
      return(!(b.p2.x() < p1.x() || b.p1.x() > p2.x()) &&
	     !(b.p2.y() < p1.y() || b.p1.y() > p2.y()) &&
	     !(b.p2.z() < p1.z() || b.p1.z() > p2.z()));
    };
    CBox intersect(const CBox &b) const
    {
      CBox retval;
      if(collide(b)){
	retval.p1.x(b.p1.x() > p1.x() ? b.p1.x() : p1.x());
	retval.p1.y(b.p1.y() > p1.y() ? b.p1.y() : p1.y());
	retval.p1.z(b.p1.z() > p1.z() ? b.p1.z() : p1.z());
	retval.p2.x(b.p2.x() < p2.x() ? b.p2.x() : p2.x());
	retval.p2.y(b.p2.y() < p2.y() ? b.p2.y() : p2.y());
	retval.p2.z(b.p2.z() < p2.z() ? b.p2.z() : p2.z());
      }
      return(retval);
    };
  };
  
  class CPlane {
  protected:
    CVector n;
    CPoint  p;
  public:
    CPlane(){};
    CPlane(const CPoint &p1,const CPoint &p2,const CPoint &p3)
      : p(p1)
    {
      double A = p1.y()*(p2.z() - p3.z()) + p2.y()*(p3.z() - p1.z()) +
	p3.y()*(p1.z() - p2.z());
      double B = p1.z()*(p2.x() - p3.x()) + p2.z()*(p3.x() - p1.x()) +
	p3.z()*(p1.x() - p2.x());
      double C = p1.x()*(p2.y() - p3.y()) + p2.x()*(p3.y() - p1.y()) +
	p3.x()*(p1.y() - p2.y());
      //      double D = -1.0 * ( p1.x()*(p2.y()*p3.z() - p3.y()*p2.z()) +
      //			  p2.x()*(p3.y()*p1.z() - p1.y()*p3.z()) +
      //			  p3.x()*(p1.y()*p2.z() - p2.y()*p1.z()) );
      
      n.x(A).y(B).z(C);
      n.normalize();
    };
    CPlane(const CVector &N,const CPoint &P)
      : n(N.unit()),p(P)
    {};
    CPlane(const CPlane &plane)
      : n(plane.n),p(plane.p)
    {};
    CPlane(const std::vector<CPoint> &vec)
      : p(vec[0])
    {
      double A = vec[0].y()*(vec[1].z() - vec[2].z()) + 
	vec[1].y()*(vec[2].z() - vec[0].z()) +
	vec[2].y()*(vec[0].z() - vec[1].z());
      double B = vec[0].z()*(vec[1].x() - vec[2].x()) + 
	vec[1].z()*(vec[2].x() - vec[0].x()) +
	vec[2].z()*(vec[0].x() - vec[1].x());
      double C = vec[0].x()*(vec[1].y() - vec[2].y()) + 
	vec[1].x()*(vec[2].y() - vec[0].y()) +
	vec[2].x()*(vec[0].y() - vec[1].y());
      //      double D = -1.0 * ( p1.x()*(p2.y()*p3.z() - p3.y()*p2.z()) +
      //			  p2.x()*(p3.y()*p1.z() - p1.y()*p3.z()) +
      //			  p3.x()*(p1.y()*p2.z() - p2.y()*p1.z()) );
      
      n.x(A).y(B).z(C);
      n.normalize();
    } ;     
    CPlane &operator=(const CPlane &plane)
    {
      n = plane.n;
      p = plane.p;
      return(*this);
    };
    double A() const { return(n.x()); };
    double B() const { return(n.y()); };
    double C() const { return(n.z()); };
    double D() const 
    { 
      return( -1.0 * ( n.x()*p.x() + n.y()*p.y() + n.z()*p.z())); 
    };
    CVector &normal() { return(n); };
    const CVector &normal() const { return(n); };
    CPoint &point() { return(p); };
    const CPoint &point() const { return(p); };
    bool contains_point(const CPoint &P) const
    {
      return(fabs(n.x()*P.x() + n.y()*P.y() + n.z()*P.z() + D()) <= TOL);
    };
    bool operator==(const CPlane &plane) const
    {
      return((plane.n.unit() == n.unit()) && contains_point(plane.p));
    };
  };
  
  class C3Plane {
  protected:
    C3Vector n;
    C3Point  p;
  public:
    C3Plane(){};
    C3Plane(const C3Point &p1,const C3Point &p2,const C3Point &p3)
      : p(p1)
    {
      double A = p1.y()*(p2.z() - p3.z()) + p2.y()*(p3.z() - p1.z()) +
	p3.y()*(p1.z() - p2.z());
      double B = p1.z()*(p2.x() - p3.x()) + p2.z()*(p3.x() - p1.x()) +
	p3.z()*(p1.x() - p2.x());
      double C = p1.x()*(p2.y() - p3.y()) + p2.x()*(p3.y() - p1.y()) +
	p3.x()*(p1.y() - p2.y());
      //      double D = -1.0 * ( p1.x()*(p2.y()*p3.z() - p3.y()*p2.z()) +
      //			  p2.x()*(p3.y()*p1.z() - p1.y()*p3.z()) +
      //			  p3.x()*(p1.y()*p2.z() - p2.y()*p1.z()) );
      
      n.x(A).y(B).z(C);
      n.Normalize();
    };
    C3Plane(const C3Vector &N,const C3Point &P)
      : n(N.Unit()),p(P)
    {};
    C3Plane(const C3Plane &plane)
      : n(plane.n),p(plane.p)
    {};
    C3Plane(const std::vector<C3Point> &vec)
      : p(vec[0])
    {
      double A = vec[0].y()*(vec[1].z() - vec[2].z()) + 
	vec[1].y()*(vec[2].z() - vec[0].z()) +
	vec[2].y()*(vec[0].z() - vec[1].z());
      double B = vec[0].z()*(vec[1].x() - vec[2].x()) + 
	vec[1].z()*(vec[2].x() - vec[0].x()) +
	vec[2].z()*(vec[0].x() - vec[1].x());
      double C = vec[0].x()*(vec[1].y() - vec[2].y()) + 
	vec[1].x()*(vec[2].y() - vec[0].y()) +
	vec[2].x()*(vec[0].y() - vec[1].y());
      //      double D = -1.0 * ( p1.x()*(p2.y()*p3.z() - p3.y()*p2.z()) +
      //			  p2.x()*(p3.y()*p1.z() - p1.y()*p3.z()) +
      //			  p3.x()*(p1.y()*p2.z() - p2.y()*p1.z()) );
      
      n.x(A).y(B).z(C);
      n.Normalize();
    } ;     
    C3Plane &operator=(const C3Plane &plane)
    {
      n = plane.n;
      p = plane.p;
      return(*this);
    };
    double A() const { return(n.x()); };
    double B() const { return(n.y()); };
    double C() const { return(n.z()); };
    double D() const 
    { 
      return( -1.0 * ( n.x()*p.x() + n.y()*p.y() + n.z()*p.z())); 
    };
    C3Vector &normal() { return(n); };
    const C3Vector &normal() const { return(n); };
    C3Point &point() { return(p); };
    const C3Point &point() const { return(p); };
    bool contains_point(const C3Point &P) const
    {
      return(fabs(n.x()*P.x() + n.y()*P.y() + n.z()*P.z() + D()) <= TOL);
    };
    bool operator==(const C3Plane &plane) const
    {
      return((plane.n.Unit() == n.Unit()) && contains_point(plane.p));
    };
  };
  
  template<typename PointContainer>
  C3Point C3Centroid(const PointContainer &pc)
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double np = 1.0/static_cast<double>(pc.size());
    typename PointContainer::const_iterator pi = pc.begin();
    while(pi != pc.end()){
      x += (*pi)->x();
      y += (*pi)->y();
      z += (*pi)->z();
      pi++;
    }
    x *= np;
    y *= np;
    z *= np;
    return(C3Point(x,y,z));
  };

  template<typename PointContainer>
  CPoint Centroid(const PointContainer &pc)
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double np = 1.0/static_cast<double>(pc.size());
    typename PointContainer::const_iterator pi = pc.begin();
    while(pi != pc.end()){
      x += pi->x();
      y += pi->y();
      z += pi->z();
      pi++;
    }
    x *= np;
    y *= np;
    z *= np;
    return(CPoint(x,y,z));
  }; 

  class CFacet : public std::vector<CPoint>
  {
  public:
    CFacet(){};
    CFacet(const CPoint &P1,const CPoint &P2,const CPoint &P3)
    {
      this->resize(3);
      (*this)[0] = P1;
      (*this)[1] = P2;
      (*this)[2] = P3;
    };
    CFacet(const CPoint &P1,const CPoint &P2,const CPoint &P3,
	   const CPoint &P4)
    {
      this->resize(4);
      (*this)[0] = P1;
      (*this)[1] = P2;
      (*this)[2] = P3;
      (*this)[3] = P4;
    };
    CFacet(const std::vector<CPoint> &v)
    {
      std::vector<CPoint>::const_iterator pi = v.begin();
      while(pi != v.end())
	this->push_back(*pi++);
    };
    CFacet(const CFacet &face)
    {
      std::vector<CPoint>::const_iterator fi = face.begin();
      while(fi != face.end())
	this->push_back(*fi++);
    };
    CVector normal() const
    {
      CVector v1((*this)[0],(*this)[1]);
      CVector v2((*this)[0],(*this)[2]);
      return(.5*(v1%v2));
    };
    CPoint centroid() const
    {
      return(Centroid(*this));
    };
    CFacet &operator=(const CFacet &face)
    {
      this->clear();
      std::vector<CPoint>::const_iterator fi = face.begin();
      while(fi != face.end())
	this->push_back(*fi++);
      return(*this);
    };
    bool operator==(const CFacet &face) const
    {
      std::vector<CPoint>::const_iterator li = this->begin();
      std::vector<CPoint>::const_iterator fi = face.begin();
      if(this->size() != face.size()) return false;
      while(li != this->end())
	if(!(*li++ == *fi++))
	  return(false);
      return(true);
    };
    bool contains_point(const CPoint &P)
    {
      CPlane plane((*this)[0],(*this)[1],(*this)[2]);
      if(!(plane.contains_point(P)))
	return false;
      std::vector<CPoint>::const_iterator fi = this->begin();
      CVector v1(P,*fi);
      CVector vp(v1);
      fi++;
      double angle = 0.0;
      while(fi != this->end()){
	CVector v(P,(*fi));
	angle += acos(v*vp);
	vp = v;
	fi++;
      }
      angle += acos(vp*v1);
      return((fabs(angle) - 2*PI) <= TOL);
    };
    
    
  };  
  

  class C3Facet : public std::vector<const C3Point *>
  {
  private:
    bool _mydata;
  public:
    C3Facet()
    {
      _mydata = false;
    };
    C3Facet(const C3Point *P1,const C3Point *P2,const C3Point *P3)
    {
      _mydata = false;
      this->resize(3);
      (*this)[0] = P1;
      (*this)[1] = P2;
      (*this)[2] = P3;
    };
    C3Facet(const C3Point *P1,const C3Point *P2,const C3Point *P3,const C3Point *P4)
    {
      _mydata = false;
      this->resize(4);
      (*this)[0] = P1;
      (*this)[1] = P2;
      (*this)[2] = P3;
      (*this)[3] = P4;
    };
    C3Facet(const C3Point &P1,const C3Point &P2,const C3Point &P3)
    {
      _mydata = true;
      this->resize(3);
      (*this)[0] = _copycreate_point(P1);
      (*this)[1] = _copycreate_point(P2);
      (*this)[2] = _copycreate_point(P3);
    };
    C3Facet(const C3Point &P1,const C3Point &P2,const C3Point &P3,
	    const C3Point &P4)
    {
      _mydata = true;
      this->resize(4);
      (*this)[0] = _copycreate_point(P1);
      (*this)[1] = _copycreate_point(P2);
      (*this)[2] = _copycreate_point(P3);
      (*this)[3] = _copycreate_point(P4);
    };
    C3Facet(const std::vector<const C3Point *> &v)
    {
      std::vector<const C3Point *>::const_iterator pi = v.begin();
      while(pi != v.end())
	this->push_back(*pi++);
    };
//     C3Facet(const std::vector<const C3Point> &v)
//     {
//       std::vector<const C3Point>::const_iterator pi = v.begin();
//       while(pi != v.end())
// 	this->push_back(_copycreate_point(*pi++));
//     };
    C3Facet(const C3Facet &face)
    {
      C3Facet::const_iterator fi = face.begin();
      while(fi != face.end())
	this->push_back(*fi++);
    };
    C3Vector Normal() const
    {
      C3Vector v1(*((*this)[0]),*((*this)[1]));
      C3Vector v2(*((*this)[0]),*((*this)[2]));
      return(.5*(v1%v2));
    };
    C3Point Centroid() const
    {
      return(C3Centroid(*this));
    };
    C3Facet &operator=(const C3Facet &face)
    {
      this->clear();
      C3Facet::const_iterator fi = face.begin();
      while(fi != face.end())
	this->push_back(*fi++);
      return(*this);
    };
    bool operator==(const C3Facet &face) const
    {
      C3Facet::const_iterator li = this->begin();
      C3Facet::const_iterator fi = face.begin();
      if(this->size() != face.size()) return false;
      while(li != this->end())
	if(!(*(*li++) == *(*fi++)))
	  return(false);
      return(true);
    };
//     bool contains_point(const C3Point &P)
//     {
//       CPlane plane(*(*this)[0],*(*this)[1],*(*this)[2]);
//       if(!(plane.contains_point(P)))
// 	return false;
//       std::vector<CPoint>::const_iterator fi = this->begin();
//       CVector v1(P,*fi);
//       CVector vp(v1);
//       fi++;
//       double angle = 0.0;
//       while(fi != this->end()){
// 	CVector v(P,(*fi));
// 	angle += acos(v*vp);
// 	vp = v;
// 	fi++;
//       }
//       angle += acos(vp*v1);
//       return((fabs(angle) - 2*PI) <= TOL);
//     };
  private:
    C3Point *_copycreate_point(const C3Point &p)
    {
      C3Point *rp = new C3Point;
      rp->Copy(p);
      return(rp);
    };
    C3Point *_create_point()
    {
      C3Point *rp = new C3Point;
      return(rp);
    };
    
  };  
  
  void 
  Transpose(CVector matrix[]);
  void 
  Transpose_2x3(CVector matrix[],double tpose[][2]);
  
  
  double
  Distance(const CPoint &p,const CLine &l);

}
#endif // GeoPrimitives
