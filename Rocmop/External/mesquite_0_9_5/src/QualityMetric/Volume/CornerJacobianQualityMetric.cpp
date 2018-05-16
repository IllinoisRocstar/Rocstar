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
/*!
  \file   CornerJacobianQualityMetric.cpp
  \brief  

  \author Michael Brewer
  \date   2002-06-9
*/
#include "CornerJacobianQualityMetric.hpp"
#include "Vector3D.hpp"
#include "ShapeQualityMetric.hpp"
#include "QualityMetric.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif

using namespace Mesquite;

CornerJacobianQualityMetric::CornerJacobianQualityMetric()
{
  MsqError err;
  set_metric_type(ELEMENT_BASED);
  set_element_evaluation_mode(ELEMENT_VERTICES, err); MSQ_CHKERR(err);
  avgMethod=QualityMetric::NONE;
  feasible=0;
  set_name("Corner Jacobian Volume (or Area)");
}

bool CornerJacobianQualityMetric::evaluate_element(PatchData &pd,
                                                    MsqMeshEntity *element,
                                                    double &fval,
                                                   MsqError &err)
{
  switch(element->get_element_type()){
    case TRIANGLE:
    case QUADRILATERAL:
      fval=element->compute_unsigned_area(pd,err);  MSQ_ERRZERO(err);
      break;
    case TETRAHEDRON:
    case HEXAHEDRON:
      fval=element->compute_unsigned_volume(pd,err);  MSQ_ERRZERO(err);
      break;
    default:
      fval=MSQ_MAX_CAP;
  }// end switch over element type
  return true;
}


