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

/*! \file ConditionNumberQualityMetric.hpp

Header file for the Mesquite::ConditionNumberQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef ConditionNumberQualityMetric_hpp
#define ConditionNumberQualityMetric_hpp


#include "Mesquite.hpp"
#include "ShapeQualityMetric.hpp"


namespace Mesquite
{
     /*! \class ConditionNumberQualityMetric
       \brief Computes the condition number of given element.

       The metric does not use the sample point functionality or the
       compute_weighted_jacobian.  It evaluates the metric at
       the element vertices, and uses the isotropic ideal element.
       It does require a feasible region, and the metric needs
       to be minimized.
     */
   class ConditionNumberQualityMetric : public ShapeQualityMetric
   {
  public:
     ConditionNumberQualityMetric();
     
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~ConditionNumberQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  protected:
     
  private:
    
  };
    
   

} //namespace


#endif // ConditionNumberQualityMetric_hpp


