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
  \file   ScalarMultiplyQualityMetric.cpp
  \brief 

  \author Todd Munson
  \date   2004-12-21
*/

#include "ScalarMultiplyQualityMetric.hpp"
#include "QualityMetric.hpp"
#include "Vector3D.hpp"

using namespace Mesquite;

ScalarMultiplyQualityMetric::ScalarMultiplyQualityMetric(QualityMetric* qm1,
                                       			 double scalar_double,
                                       			 MsqError &err)
  : qualMetric(qm1), scaleFactor(scalar_double)
{
  if(qm1 == NULL){
    MSQ_SETERR(err)("NULL pointer", MsqError::INVALID_ARG);
    return;
  }
  feasible=qm1->get_feasible_constraint();
  int n_flag=qm1->get_negate_flag();
  set_negate_flag(n_flag);
  
  set_metric_type(qm1->get_metric_type());
  
  set_name("Scalar Multiply Quality Metric");
}

/*! Returns scalar times qMetric1->evaluate_element(element, err).
 */
bool ScalarMultiplyQualityMetric::evaluate_element(PatchData& pd,
                                                   MsqMeshEntity *element,
                                                   double &value,
                                                   MsqError &err)
{
  bool valid_flag;
  double metric1;
  valid_flag=qualMetric->evaluate_element(pd, element, metric1, err); MSQ_ERRZERO(err);
  if(!valid_flag)
    return false;
  value = scaleFactor*metric1;
  return valid_flag;
}

/*! Returns scalar times qMetric1->evaluate_vertex(...).
 */
bool ScalarMultiplyQualityMetric::evaluate_vertex(PatchData& pd,
                                                  MsqVertex* vert,
                                                  double &value,
                                                  MsqError& err)
{
  bool valid_flag;
  double metric1;
  valid_flag=qualMetric->evaluate_vertex(pd, vert, metric1, err); MSQ_ERRZERO(err);
  if(!valid_flag)
    return false;
  value = scaleFactor*metric1;
  return valid_flag;
}

