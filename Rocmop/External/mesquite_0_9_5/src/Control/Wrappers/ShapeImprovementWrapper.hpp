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
//
//   SUMMARY: 
//     USAGE:
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 14-Nov-02 at 16:51:36
//  LAST-MOD: 22-May-03 at 09:04:07 by Michael Brewer


/*! \file ShapeImprovementWrapper.hpp
  ShapeImprovementWrapper header file.

*/
// DESCRIP-END.
//


#ifndef ShapeImprovementWrapper_hpp
#define ShapeImprovementWrapper_hpp

#include "IdealWeightInverseMeanRatio.hpp" 
#include "FeasibleNewton.hpp"
#include "LPtoPTemplate.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"

#include "UntangleBetaQualityMetric.hpp"
#include "ConjugateGradient.hpp"

namespace Mesquite { 
  /*! \class ShapeImprovementWrapper
       \brief Wrapper which performs a Feasible Newton solve using
       an \f$\ell_2^2 \f$ objective function template with inverse mean
       ratio.
       
     */
  class ShapeImprovementWrapper : public InstructionQueue {
     
  public:  
      //Constructor sets the instructions in the queue.
    ShapeImprovementWrapper(MsqError& err,
                            double cpu_time = 0.0, 
                            double grad_norm =1.e-6);
    
      //! Destructor must delete the objects inserted in the queue.
    virtual ~ShapeImprovementWrapper();
    
      //! run_instructions runs the wrapper on the given MeshSet.
    virtual void run_instructions(MeshSet &ms, MsqError &err);
    
    
  private:

    ShapeQualityMetric* inverseMeanRatio; 
    LPtoPTemplate* objFunc;
    FeasibleNewton* feasNewt;
    QualityAssessor* mQA;
    TerminationCriterion* termOuter; 
    TerminationCriterion* termInner;

    bool timerNeeded;
    double maxTime;
      //arbitraryily chosen variables
    double untBeta;
    double successiveEps;
  

  };
  
  
} // namespace

#endif // ShapeImprovementWrapper_hpp
