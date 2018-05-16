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
  \file    CompositeOFMultiply.cpp
  \brief  

  This Objective Function combines two Objective Functions by mulitplication
  \author Michael Brewer
  \date   2002-01-23
*/
#include <math.h>
#include "CompositeOFMultiply.hpp"
#include "MsqTimer.hpp"
#include "PatchData.hpp"

namespace Mesquite {


/*!
Sets the QualityMetric pointer to the metric associated with Obj1 and Obj2
if Obj1 and Obj2 are associated with the same metric.  Otherwise, it sets
the QualityMetric pointer to NULL.  The new
ObjectiveFunction's negateFlag is set to negative one only if both Obj1 and
Obj2's negateFlag are negative one (because obj1 and obj2's evaluate function
multiply their return values by negative one if their respective
function needs to be maximized.  If both of these functions needed to
be maximized, then the negative ones will have cancelled out).  Otherwise,
the negateFlag is set to one.  Defaults to the analytical gradient.
  \param Obj1 (ObjectiveFunction*)
  \param Obj2 (ObjectiveFunction*)
 */
CompositeOFMultiply::CompositeOFMultiply(ObjectiveFunction* Obj1, ObjectiveFunction* Obj2){
  if(Obj1->get_quality_metric()==Obj2->get_quality_metric()){
    set_quality_metric(Obj1->get_quality_metric());
  }
  else
    set_quality_metric(NULL);
  objFunc1=Obj1;
  objFunc2=Obj2;
    //if both obj1 and ob2 have been negated
  if(Obj1->get_negate_flag()-Obj2->get_negate_flag()==-2)
    set_negate_flag(-1);
  else
    set_negate_flag(1);
  set_gradient_type(ObjectiveFunction::ANALYTICAL_GRADIENT);
}

//Michael:  need to clean up here
CompositeOFMultiply::~CompositeOFMultiply(){

}

/*! Returns the QualityMetric list associated with objFunc1 merged
  with the QualityMetric list associated with objFunc2.  The entries
  in this merged list may not be unique.
*/
msq_std::list<QualityMetric*> CompositeOFMultiply::get_quality_metric_list()
{
  msq_std::list<QualityMetric*> temp_list=objFunc1->get_quality_metric_list();
  msq_std::list<QualityMetric*> temp_list2=objFunc2->get_quality_metric_list();
  temp_list.merge(temp_list2);
  return temp_list;
}


/*!Computes fval= objFunc1->evaluate(patch,err)*objFunc2->evaluate(patch,err).
  Note that since objFunc1 and objFunc2's evaluate() functions are called
  (as opposed to their concrete_evaluates) the returned values have
  already been multiplied by the respective negateFlag (that is,
  if objFunc1 (or objFunc2) needed to be maximized then the value has
  been multiplied by negative one so that it may be minimized instead.)
  Function returns `false' if either objFunc1->evaluate() or
  objFunc2->evaluate() returns `false'; otherwise function returns `true'.
*/
bool CompositeOFMultiply::concrete_evaluate(PatchData &patch, double &fval,
                                            MsqError &err){
  double second_val;
    //if invalid, return false without calculating fval.
  bool b = objFunc1->evaluate(patch, fval, err); 
  if(MSQ_CHKERR(err) || !b){
    fval=0.0;
    return false;
  }
  
  b = objFunc2->evaluate(patch, second_val, err); 
  if(MSQ_CHKERR(err) || !b){
    fval=0.0;
    return false;
  }
  fval*=second_val;
  return true;
}

/*! Analytically computes the composite objective function's gradient,
  using the multiplication rule.
  by scaling the gradient returned objFunc->compute_gradient().
    \param patch The PatchData object for which the objective function
           gradient is computed.
    \param grad An array of Vector3D, at least the size of the number
           of vertices in the patch.
    \param array_size is the size of the grad Vector3D[] array and
    must correspond to the number of vertices in the patch.
*/
bool CompositeOFMultiply::compute_analytical_gradient(PatchData &patch,
						      Vector3D *const &grad,
						      double &OF_val,
						      MsqError &err,
						      size_t array_size)
{
  MSQ_FUNCTION_TIMER( "CompositeOFMultiply::compute_analytical_gradient" );
  
  double obj_2_val=0.0;
  //get the first gradient and objective function value
  bool rval=objFunc1->compute_gradient(patch, grad, OF_val,  err, array_size); MSQ_ERRZERO(err);
  //if the above is valid, get the second gradient
  if(rval){
    int num_vert=patch.num_vertices();
    Vector3D* second_grad = new Vector3D[num_vert];
    //get second objective function's gradient
    rval=objFunc2->compute_gradient(patch, second_grad,obj_2_val, 
				    err, num_vert);
    //if both objective functions gradients were successfully computed, 
    //use the multiplaction rule to get the complete gradient.
    if(rval){
      int i=0;
      for(i=0;i<num_vert;++i){
        grad[i]*=obj_2_val;
        grad[i]+=(second_grad[i]*OF_val);
      }
    }
    //delete the dynamically allocated space for the second gradient
    delete []second_grad;
  }
  
    //true if both gradient and both evaluate were successful.
  //compute the objective function value by mulitiplying
  //OF_val and obj_2_val
  if(rval)
    OF_val*=obj_2_val;
  else
    OF_val=0.0;
  return rval;
}
	
	
} // namespace Mesquite
