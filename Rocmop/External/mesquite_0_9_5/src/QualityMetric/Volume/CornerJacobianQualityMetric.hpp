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

/*! \file CornerJacobianQualityMetric.hpp

Header file for the Mesquite::CornerJacobianQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef CornerJacobianQualityMetric_hpp
#define CornerJacobianQualityMetric_hpp

#include "Mesquite.hpp"
#include "VolumeQualityMetric.hpp"


namespace Mesquite
{
     /*! \class CornerJacobianQualityMetric
       \brief Computes the volume or area of the element, as appropriate.
       This metric uses the average of the corner Jacobian determinants
       for the approximation to the volume of hex.

        The metric does not use the sample point functionality or the
        compute_weighted_jacobian (except for possibly indirectly
        when evaluating the metric value for a hex).  It evaluates
        the signed area of surface elements and the signed volume
        of volume elements.  It does require a feasible region,
        and (in general) the metric needs to be minimized.
     */
     /*!\todo MB:  CornerJacobianQualityMetric is currently not being
       evaluated using the corner jacobian method (instead it is using
       the area functions from MsqMeshEntity... needs to be modified
       to do so.*/
   class CornerJacobianQualityMetric : public VolumeQualityMetric
   {
  public:
 
     CornerJacobianQualityMetric();

       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~CornerJacobianQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  };

} //namespace


#endif // CornerJacobianQualityMetric_hpp


