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
  \file    CompositeOFScalarAdd.cpp
  \brief  

  This Objective Function combines two Objective Functions by addition
  \author Michael Brewer
  \date   2002-06-24
*/
#include <math.h>
#include "CompositeOFScalarAdd.hpp"
#include "PatchData.hpp"
#include "MsqTimer.hpp"
using namespace Mesquite;


/*!
Sets the QualityMetric pointer to the metric associated with Obj1.  
The new ObjectiveFunction's negateFlag is also the
same as that of Obj1.  This objective function defaults to the analytical
gradient which essentially just calls Obj1's gradient function.
  \param alp (double)
  \param Obj1 (ObjectiveFunction*)
 */
CompositeOFScalarAdd::CompositeOFScalarAdd(double alp, ObjectiveFunction* Obj1){
  set_quality_metric(Obj1->get_quality_metric());
  objFunc=Obj1;
  mAlpha=alp;
  set_negate_flag(1);
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}


//Michael:  need to clean up here
CompositeOFScalarAdd::~CompositeOFScalarAdd(){

}

/*!Computes fval= mAlpha+objFunc->evaluate(patch,err).  Note that since Obj's
  evaluate() function is called (as opposed to its concrete_evaluate) the
  returned value has been multiplied by objFunc's negateFlag (that is,
  if objFunc needed to be maximized then the value has been multiplied
  by negative one so that it may be minimized instead.)
  Functions returns `false' if and only if objFunc->evaluate() returns
  `false'.
*/
bool CompositeOFScalarAdd::concrete_evaluate(PatchData &patch, double &fval,
                                             MsqError &err){
    //if invalid return false without calculating fval.
  bool b = objFunc->evaluate(patch, fval, err);
  if( MSQ_CHKERR(err) || !b ){
    fval = 0.0;
    return false;
  }
  
  fval+=mAlpha;
  return true;
}

//!Returns the QualityMetric list assossiated with objFunc.
msq_std::list<QualityMetric*> CompositeOFScalarAdd::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}

/*! Analytically computes the composite objective function's gradient
  by scaling the gradient returned objFunc->compute_gradient().
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFScalarAdd::compute_analytical_gradient(PatchData &patch,
                                                       Vector3D *const &grad,
                                                       double &OF_val,
                                                       MsqError &err,
                                                       size_t array_size)
{
  MSQ_FUNCTION_TIMER( "CompositeOFScalarAdd::compute_analytical_gradient" );
  bool rval=objFunc->compute_gradient(patch, grad, OF_val,err, array_size); MSQ_ERRZERO(err);
  OF_val+=mAlpha;
  return rval;
}

