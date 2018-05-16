/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-

/*! \file PowerQualityMetric.hpp
\brief
Header file for the Mesquite::PowerQualityMetric class

  \author Michael Brewer
  \date   2002-09-05
 */


#ifndef PowerQualityMetric_hpp
#define PowerQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"
#include "Exponent.hpp"

namespace Mesquite
{
  class MsqError;
  class Vector3D;

     /*! \class PowerQualityMetric
       \brief Raises a single quality metrics (qMetric1) to an arbitrary
       power (a double value, scaleAlpha) for two- and three-diminsional
       elements.  
     */
   class PowerQualityMetric : public QualityMetric
   {
  public:
       /*! Ensures that qm1 is not NULL.  If qm1 is only valid
         on a certain feasible, then the composite metric has the same
         constraint.  The composite metric also has the same negate flag
         as qm1.
       */
     PowerQualityMetric(QualityMetric* qm1, double pow_double,MsqError &err);
     
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~PowerQualityMetric()
        {  }
     
     bool evaluate_element(PatchData& pd, MsqMeshEntity *element,double &value,
                             MsqError &err);
     bool evaluate_vertex(PatchData& pd, MsqVertex *vertex, double &value,
                            MsqError &err);

  private:
  
    QualityMetric* qualMetric;
    Mesquite::Exponent mPower;
   };
   

} //namespace


#endif // PowerQualityMetric_hpp







