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

/*! \file ShapeImprovementWrapper.cpp

Member functions of the Mesquite::ShapeImprovementWrapper class

  \author Michael Brewer
  \date   June 6, 2003
 */

#include "InstructionQueue.hpp"
#include "ShapeImprovementWrapper.hpp"
#include "MeshSet.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {

/*! The consturctor allows for two values.  The first is a 
  time bound (in seconds) used as a termination criterion.  If
  this value is non-positive, no time bound will be set.
  By default, the value is set to zero and no time bound
  is used.  The second value is the tolerance for the gradient
  norm termination criteria.  The default value is 1.e-6.*/
ShapeImprovementWrapper::ShapeImprovementWrapper(MsqError& err,
                                                 double cpu_time,
                                                 double grad_norm) 
 : inverseMeanRatio(0),
   objFunc(0),
   feasNewt(0),
   mQA(0),
   termOuter(0),
   termInner(0)
{

    //arbitrarily chosen variables
  untBeta=1.e-8;
  successiveEps=1.e-4;  
  
  inverseMeanRatio = new IdealWeightInverseMeanRatio(err); MSQ_ERRRTN(err);
  inverseMeanRatio->set_gradient_type(QualityMetric::ANALYTICAL_GRADIENT);
  inverseMeanRatio->set_hessian_type(QualityMetric::ANALYTICAL_HESSIAN);
  inverseMeanRatio->set_averaging_method(QualityMetric::LINEAR,err);  MSQ_ERRRTN(err);
    // creates the l_2 squared objective function
  objFunc = new LPtoPTemplate(inverseMeanRatio, 2, err);  MSQ_ERRRTN(err);
  objFunc->set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
    //creates a FeasibleNewtone improver
  feasNewt = new FeasibleNewton(objFunc);
  feasNewt->set_patch_type(PatchData::GLOBAL_PATCH, err,1 ,1);  MSQ_ERRRTN(err);

  mQA = new QualityAssessor(inverseMeanRatio,QualityAssessor::MAXIMUM, err); MSQ_ERRRTN(err);
  mQA->disable_printing_results();
          //**************Set stopping criterion*e***************
  termInner = new TerminationCriterion();
  termOuter = new TerminationCriterion();
  //termInner->add_criterion_type_with_double(TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE,grad_norm,err);  MSQ_ERRRTN(err);
  //termInner->add_criterion_type_with_double(TerminationCriterion::SUCCESSIVE_IMPROVEMENTS_RELATIVE,successiveEps,err);  MSQ_ERRRTN(err);
  termInner->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err); MSQ_ERRRTN(err);
  termOuter->add_criterion_type_with_int(TerminationCriterion::NUMBER_OF_ITERATES,1,err); MSQ_ERRRTN(err);
  feasNewt->add_culling_method(PatchData::NO_BOUNDARY_VTX);
  feasNewt->set_inner_termination_criterion(termInner);
  feasNewt->set_outer_termination_criterion(termOuter);
      
}


ShapeImprovementWrapper::~ShapeImprovementWrapper()
{
  delete inverseMeanRatio;
  delete objFunc;
  delete feasNewt;
  delete mQA;
  delete termInner;
  delete termOuter;
}


/*!Run instructions first calls the global untangler.  If the
  resulting mesh is tangled after that pre-conditioning step,
  The mesh is iteratively smoothed with a local and then global
  untangler until the mesh is untangled or until a certain time
  constraint has been exceeded.  If the mesh was successfully
  untangled and there is still time remaining, an inverse mean ratio
  shape improvement is then performed.*/
void ShapeImprovementWrapper::run_instructions(MeshSet &ms, MsqError &err)
{

  mQA->loop_over_mesh(ms, err);  
  MSQ_ERRRTN(err);

  feasNewt->loop_over_mesh(ms, err);  
  MSQ_ERRRTN(err);

  mQA->loop_over_mesh(ms, err);  
  MSQ_ERRRTN(err);

}

} // namespace Mesquite

  
