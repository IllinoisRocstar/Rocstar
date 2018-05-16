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
  \file    CompositeOFScalarMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "CompositeOFScalarMultiply.hpp"
#include "MsqTimer.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include "PatchData.hpp"

namespace Mesquite {


/*!
Sets the QualityMetric pointer to the metric associated with Obj. 
However, if alp is less than zero, the new
ObjectiveFunction's negateFlag is the opposite of Obj's.  This objective
function defaults to the analytical gradient.
  \param alp (double)
  \param Obj (ObjectiveFunction*)
 */
CompositeOFScalarMultiply::CompositeOFScalarMultiply(double alp, ObjectiveFunction* Obj){

  set_quality_metric(Obj->get_quality_metric());
  objFunc=Obj;
  mAlpha=alp;
  if(alp<0)
    set_negate_flag(-1);
  else if (alp>0)
    set_negate_flag(1);
  else
    MSQ_DBGOUT(1) << "ObjectiveFunction being scaled by zero.";
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}

//Michael:  need to clean up here
CompositeOFScalarMultiply::~CompositeOFScalarMultiply(){

}

/*!Computes fval= mAlpha*objFunc->evaluate(patch,err).  Note that since Obj's
  evaluate() function is called (as opposed to its concrete_evaluate) the
  returned value has been multiplied by objFunc's negateFlag (that is,
  if objFunc needed to be maximized then the value has been multiplied
  by negative one so that it may be minimized instead.)
  Function returns `false' if and only if objFunc->evaluate() returns `false'.
*/
bool CompositeOFScalarMultiply::concrete_evaluate(PatchData &patch,
                                                  double &fval, MsqError &err){
    //if invalid return false without calculating fval
  bool b = objFunc->evaluate(patch, fval, err);
  if (MSQ_CHKERR(err) || !b){
    fval = 0.0;
    return false;
    
  }
  fval*=mAlpha;
  return true;
}

//!Returns the QualityMetric list assossiated with objFunc.
msq_std::list<QualityMetric*> CompositeOFScalarMultiply::get_quality_metric_list(){
  return objFunc->get_quality_metric_list();
}
	

/*! Analytically computes the composite objective function's gradient
  by scaling the gradient returned objFunc->compute_gradient().  If
  mAlpha is less than zero, the gradient is scaled by negatvie mAlpha
  because the ObjectiveFunction is multiplied by negative one so that
  it may be minimized.
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFScalarMultiply::compute_analytical_gradient(PatchData &patch,
                                                            Vector3D *const
                                                            &grad,
                                                            double &OF_val,
                                                            MsqError &err,
                                                            size_t array_size)
{
  MSQ_FUNCTION_TIMER( "CompositeOFScalarMultiply::compute_analytical_gradient" );
  
  double scale_factor=(get_negate_flag()*mAlpha);
  bool rval=objFunc->compute_gradient(patch, grad, OF_val, err, array_size); MSQ_ERRZERO(err);
  int num_vert=patch.num_vertices();
  int i=0;
    //If the objFunc was successful in calculating the gradient
  if(rval){
      //scale the gradient by alpha
    for(i=0;i<num_vert;++i){
      grad[i]*=scale_factor;
    }
    //scale the OF_val
    OF_val*=mAlpha;
  }
  else{
    OF_val=0.0;
  }
  return rval;
}

} // namespace Mesquite
