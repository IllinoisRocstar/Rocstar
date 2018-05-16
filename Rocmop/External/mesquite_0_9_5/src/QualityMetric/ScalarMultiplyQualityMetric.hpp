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

/*! \file ScalarMultiplyQualityMetric.hpp
\brief
Header file for the Mesquite::ScalarMultiplyQualityMetric class

  \author Todd Munson
  \date   2004-12-21
 */


#ifndef ScalarMultiplyQualityMetric_hpp
#define ScalarMultiplyQualityMetric_hpp

#include "Mesquite.hpp"
#include "QualityMetric.hpp"

namespace Mesquite
{
     /*! \class ScalarMultiplyQualityMetric
       \brief Multiplies quality metric value by a number (a double).
     */
   class ScalarMultiplyQualityMetric : public QualityMetric
   {
  public:
       /*! Ensures that qm1 is not NULL.  If qm1 is only valid
         on a certain feasible, then the composite metric has the same
         constraint.  The composite metric also has the same negate flag
         as qm1.
       */
     ScalarMultiplyQualityMetric(QualityMetric* qm1, double scalar_double,
                                 MsqError &err);
     
       // virtual destructor ensures use of polymorphism during destruction
     virtual ~ScalarMultiplyQualityMetric()
        {  }
     
     bool evaluate_element(PatchData& pd, MsqMeshEntity *element,double &value,
                           MsqError &err);
     bool evaluate_vertex(PatchData& pd, MsqVertex *vertex, double &value,
                          MsqError &err);

  private:
    
    QualityMetric* qualMetric;
    double scaleFactor;
          
   };
   

} //namespace


#endif // ScalarMultiplyQualityMetric_hpp


