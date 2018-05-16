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

/*! \file ASMQualityMetric.hpp

Header file for the Mesquite::ASMQualityMetric class

  \author Michael Brewer
  \date   2002-06-19
 */


#ifndef ASMQualityMetric_hpp
#define ASMQualityMetric_hpp

#include "Mesquite.hpp"
#include "SmoothnessQualityMetric.hpp"


namespace Mesquite
{
     /*! \class ASMQualityMetric
       \brief Computes the ASM (Area Smoothness Metric) of the given element.
       
        The metric does not use the sample point functionality or the
        compute_weighted_jacobian.  It computes the unsigned area
        or volume (a_0)of the element and of neighboring elements (a_i).
        In 2-d, elements are considered neighbors if they share an
        edge.  In 3-d, elements are considered neighbors if they
        share a face.  The metric is then taken to be the average
        (default is MAXIMUM) of abs(a_i-a_0)/(a_i+a_0) , and the metric
        needs to be minimized.  The ideal metic value is zero which occurs
        when an element and its neighbors have the same unsigned area
        (or volume in three-dimensions).
     */
   class ASMQualityMetric : public SmoothnessQualityMetric
   {
  public:
     ASMQualityMetric();

     //! virtual destructor ensures use of polymorphism during destruction
     virtual ~ASMQualityMetric()
        {}
     
       //! evaluate using mesquite objects 
     bool evaluate_element(PatchData &pd, MsqMeshEntity *element,double &fval,
                           MsqError &err); 
          
  };

   
   

} //namespace


#endif // ASMQualityMetric_hpp


