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

/*! \file AspectRatioGammaQualityMetric.hpp
  \brief
  Header file for the Mesquite::AspectRatioGammaQualityMetric class

  \author Michael Brewer
  \date   2002-05-16
 */


#ifndef AspectRatioGammaQualityMetric_hpp
#define AspectRatioGammaQualityMetric_hpp

#include "Mesquite.hpp"
#include "MsqError.hpp"
#include "ShapeQualityMetric.hpp"
namespace Mesquite
{
     /*! \class AspectRatioGammaQualityMetric
       \brief Object for computing the aspect ratio gamma of
       simplicial elements.
     */
   class AspectRatioGammaQualityMetric : public ShapeQualityMetric
   {
   public:     
     AspectRatioGammaQualityMetric()
        {
          MsqError err;
          set_metric_type(ELEMENT_BASED);
          set_element_evaluation_mode(ELEMENT_VERTICES, err);
          fourDivRootThree=4.0/sqrt(3.0);
          twelveDivRootTwo=12.0/sqrt(2.0);
          feasible=0;
          set_name("Aspect Ratio Gamma");
        }
       //! virtual destructor ensures use of polymorphism during destruction
     virtual ~AspectRatioGammaQualityMetric()
        {}
     
   protected:
     
     
   private:
       //constants used in metric calculations
     double fourDivRootThree;
     double twelveDivRootTwo;
       //!Computes the aspect ratio gamma of element.  If element
       //!is not a tetrahedron or triangle, sets an error.
     bool evaluate_element(PatchData& pd,
                           MsqMeshEntity* element, double &fval,
                           MsqError &err);
   };
   
   
} //namespace


#endif // AspectRatioGammaQualityMetric_hpp


