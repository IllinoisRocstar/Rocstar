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
//  LAST-MOD: 23-Jul-03 at 18:06:13 by Thomas Leurent


/*! \file LaplacianIQ.hpp

This is the second possibility for wrappers. It is based on the InctructionQueue concept. 

 */
// DESCRIP-END.
//


#ifndef LaplacianIQ_hpp
#define LaplacianIQ_hpp

#include "IdealWeightInverseMeanRatio.hpp" 
#include "LaplacianSmoother.hpp"
#include "QualityAssessor.hpp"
#include "InstructionQueue.hpp"
#include "TerminationCriterion.hpp"

namespace Mesquite { 

   class LaplacianIQ : public InstructionQueue {
   private:
      ShapeQualityMetric* inverseMeanRatio;
      LaplacianSmoother* lapl1;
      QualityAssessor* mQA;
      TerminationCriterion* mTerm;

   public:
      
      //! Constructor sets the instructions in the queue.  
      LaplacianIQ() {
         MsqError err;
         // creates a mean ratio quality metric ...
         inverseMeanRatio = new IdealWeightInverseMeanRatio(err);
     
         // creates the laplacian smoother  procedures
         lapl1 = new LaplacianSmoother(err);
         mQA = new QualityAssessor(inverseMeanRatio,QualityAssessor::MAXIMUM, err);
     
         //**************Set stopping criterion****************
         mTerm = new TerminationCriterion();
         mTerm->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,10,err);
 
            lapl1->set_outer_termination_criterion(mTerm);
            // sets a culling method on the first QualityImprover
            lapl1->add_culling_method(PatchData::NO_BOUNDARY_VTX);
      
            // adds 1 pass of pass1 
            this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
            this->set_master_quality_improver(lapl1, err); MSQ_CHKERR(err);
            this->add_quality_assessor(mQA,err); MSQ_CHKERR(err);
      }

      
      //! Destructor must delete the objects inserted in the queue.
      virtual ~LaplacianIQ()
      {
         delete inverseMeanRatio;
         delete lapl1;
         delete mQA;
         delete mTerm;
      }
  
   };


} // namespace

#endif // LaplacianIQ_hpp
