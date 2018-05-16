///
/// \file
/// \ingroup support
/// \brief Geometric primitive implementation
///
#include "GeoPrimitives.H"


namespace GeoPrim {
  void 
  Transpose(CVector matrix[])
  {
    CVector tpose[3];
    
    tpose[0].init(matrix[0].x(),matrix[1].x(),matrix[2].x());
    tpose[1].init(matrix[0].y(),matrix[1].y(),matrix[2].y());
    tpose[2].init(matrix[0].z(),matrix[1].z(),matrix[2].z());
    matrix[0].init(tpose[0]);
    matrix[1].init(tpose[1]);
    matrix[2].init(tpose[2]);
  }
  
  void 
  Transpose_2x3(CVector matrix[],double tpose[][2])
  {
    tpose[0][0] = matrix[0].x();
    tpose[1][0] = matrix[0].y();
    tpose[2][0] = matrix[0].z();
    tpose[0][1] = matrix[1].x();
    tpose[1][1] = matrix[1].y();
    tpose[2][1] = matrix[1].z();
  }
  
  double
  Distance(const CPoint &p,const CLine &l)
  {
    double u = (((p.x() - l.p.x())*l.v.x() + (p.y() - l.p.y())*l.v.y() + 
		 (p.z() - l.p.z())*l.v.z())/l.v.mag2());
    CPoint p2(l.p.x()+(u*l.v.x()),l.p.y()+(u*l.v.y()),l.p.z()+(u*l.v.z()));
    CVector v(p,p2);
    return(v.mag());
  }

  double
  Distance(const C3Point &p,const CLine &l)
  {
    double u = (((p.x() - l.p.x())*l.v.x() + (p.y() - l.p.y())*l.v.y() + 
		 (p.z() - l.p.z())*l.v.z())/l.v.mag2());
    C3Point p2(l.p.x()+(u*l.v.x()),l.p.y()+(u*l.v.y()),l.p.z()+(u*l.v.z()));
    C3Vector v(p,p2);
    return(v.Mag());
  }

  std::ostream &
  operator<<(std::ostream &oS,const GeoPrim::CBox &b)
  {
    oS << b.P1() << "  " << b.P2();
    return(oS);
  }
  
  std::ostream &
  operator<<(std::ostream &oS,const GeoPrim::CPoint &p)
  {
    oS << p.x() << "  " << p.y() << "  " << p.z();
    return(oS);
  }
  
  std::istream &
  operator>>(std::istream &iS,GeoPrim::CPoint &p)
  {
    iS >> p[0] >> p[1] >> p[2];
    return(iS);
  }
  
//   std::ostream &
//   operator<<(std::ostream &oS,const GeoPrim::C3Point &p)
//   {
//     oS << p.x() << "  " << p.y() << "  " << p.z();
//     return(oS);
//   }
  
//   std::istream &
//   operator>>(std::istream &iS,GeoPrim::C3Point &p)
//   {
//     iS >> p[0] >> p[1] >> p[2];
//     return(iS);
//   }
  
  std::ostream &
  operator<<(std::ostream &oS,const GeoPrim::CVector &v)
  {
    oS << v.x() << "  " << v.y() << "  " << v.z();
    return(oS);
  }
  
  std::istream &
  operator>>(std::istream &iS,GeoPrim::CVector &v)
  {
    iS >> v[0] >> v[1] >> v[2];
    return(iS);
  }
  
  std::ostream &
  operator<<(std::ostream &oS,const GeoPrim::C3Vector &v)
  {
    oS << v.x() << "  " << v.y() << "  " << v.z();
    return(oS);
  }
  
  std::istream &
  operator>>(std::istream &iS,GeoPrim::C3Vector &v)
  {
    iS >> v[0] >> v[1] >> v[2];
    return(iS);
  }
  
  GeoPrim::CPoint operator*(double scalar,const GeoPrim::CPoint &p)
  {
    GeoPrim::CPoint rv(p);
    rv.x() *= scalar;
    rv.y() *= scalar;
    rv.z() *= scalar;
    return(rv);
  }
  
  GeoPrim::C3Vector operator*(double scalar,const GeoPrim::C3Vector &v)
  {
    GeoPrim::C3Vector rv(v);
    return(rv *= scalar);
  }
  
  GeoPrim::CVector operator*(double scalar,const GeoPrim::CVector &v)
  {
    GeoPrim::CVector rv(v);
    return(rv *= scalar);
  }
  
}
