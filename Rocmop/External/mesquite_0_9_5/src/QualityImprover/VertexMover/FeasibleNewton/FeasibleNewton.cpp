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
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 15-Jan-03 at 08:05:56
//  LAST-MOD: 15-Jun-04 at 15:45:00 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*!
  \file   FeasibleNewton.cpp
  \brief  

  Implements the FeasibleNewton class member functions.
  
  \author Thomas Leurent
  \author Todd Munson
  \date   2003-01-15
*/
// DESCRIP-END.
//


#include "FeasibleNewton.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqDebug.hpp"

using namespace Mesquite;

FeasibleNewton::FeasibleNewton(ObjectiveFunction* of) :
  VertexMover()
{
  coordsMem=NULL;
  objFunc=of;
  MsqError err;
  convTol=1e-6;
  this->set_name("FeasibleNewton");
  set_patch_type(PatchData::GLOBAL_PATCH, err);
  TerminationCriterion* default_crit=get_inner_termination_criterion();
  default_crit->add_criterion_type_with_double(
                    TerminationCriterion::GRADIENT_L2_NORM_ABSOLUTE, 5e-5, err);
  MSQ_CHKERR(err);
}  
  

void FeasibleNewton::initialize(PatchData &pd, MsqError &err)
{
  // Cannot do anything.  Variable sizes with maximum size dependent
  // upon the entire MeshSet.
  coordsMem = pd.create_vertices_memento(err); MSQ_CHKERR(err);
}

void FeasibleNewton::initialize_mesh_iteration(PatchData &pd, MsqError &/*err*/)
{
  pd.reorder();
}

void FeasibleNewton::optimize_vertex_positions(PatchData &pd, 
                                               MsqError &err)
{
  MSQ_FUNCTION_TIMER( "FeasibleNewton::optimize_vertex_positions" );
  MSQ_DBGOUT(2) << "\no  Performing Feasible Newton optimization.\n";

  const double sigma   = 1e-4;
  const double beta0   = 0.25;
  const double beta1   = 0.80;
  const double tol1    = 1e-8;
  const double tol2    = 1e-12;
  const double epsilon = 1e-10;
  double original_value, new_value;
  double beta;
  
  int nv = pd.num_vertices();
  msq_std::vector<Vector3D> grad_vect(nv), d_vect(nv);
  Vector3D* grad = &grad_vect[0];
  Vector3D* d = &d_vect[0];
  bool fn_bool=true;// bool used for determining validity of patch

  int i;

  // TODD -- Don't blame the code for bad things happening when using a
  //         bad termination test or requesting more accuracy than is
  //	     possible.
  //
  //         Also, 

  // 1.  Allocate a hessian and calculate the sparsity pattern.
  mHessian.initialize(pd, err); MSQ_ERRRTN(err);

  // 3.  Calculate the norm of the gradient for the patch
  MSQ_DBGOUT(3) << "  o  gradient norm: " << length(grad, nv) << msq_stdio::endl;
  
  // does the Feasible Newton iteration until stopping is required.
  // Terminate when inner termination criterion signals.

  /* Computes the value of the stopping criterion*/
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  while ( !term_crit->terminate() ) {
    fn_bool = objFunc->compute_hessian( pd, mHessian, grad, original_value, err );
    MSQ_ERRRTN(err);
    if (!fn_bool) {
      MSQ_SETERR(err)("invalid patch for hessian calculation", MsqError::INTERNAL_ERROR);
      return; 
    }
  
    // Prints out free vertices coordinates. 
    if (MSQ_DBG(3)) {
      MSQ_DBGOUT(3) << "\n  o Free vertices ("<< pd.num_free_vertices(err)
                <<")original coordinates:\n ";
      MSQ_ERRRTN(err);
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_ERRRTN(err);
      MsqFreeVertexIndexIterator ind1(&pd, err); MSQ_ERRRTN(err);
      ind1.reset();
      while (ind1.next()) {
        MSQ_DBGOUT(3) << "\t\t\t" << toto1[ind1.value()];
      }
    }
      
    // 4. Calculate a direction using preconditionned conjugate gradients
    //    to find a zero of the Newton system of equations (H*d = -g)
    //    (a) stop if conjugate iteration limit reached
    //    (b) stop if relative residual is small
    //    (c) stop if direction of negative curvature is obtained

    mHessian.cg_solver(d, grad, err); MSQ_ERRRTN(err);

    // 5. Check for descent direction (inner produce of gradient and
    //    direction is negative.
    double alpha = inner(grad, d, nv);
    // TODD -- Add back in if you encounter problems -- do a gradient
    //         step if the direction from the conjugate gradient solver
    //         is not a descent direction for the objective function.  We
    //         SHOULD always get a descent direction from the conjugate
    //         method though, unless the preconditioner is not positive
    //         definite.
    // If direction is positive, does a gradient (steepest descent) step.

    if (alpha > -epsilon) {
      MSQ_PRINT(1)("Newton direction not guaranteed descent.  Ensure preconditioner is positive definite.\n");

      // TODD: removed performing gradient step here since we will use
      // gradient if step does not produce descent.  Instead we set
      // alpha to a small negative value.

      alpha = -epsilon;

      // alpha = inner(grad, grad, nv); // compute norm squared of gradient
      // if (alpha < 1) alpha = 1;	// take max with constant
      // for (i = 0; i < nv; ++i) {
      //   d[i] = -grad[i] / alpha; 	// compute scaled gradient
      // }
      // alpha = inner(grad, d, nv);  	// recompute alpha
      // 				// equal to one for large gradient
    }
    
    alpha *= sigma;
    beta = 1.0;
    pd.recreate_vertices_memento(coordsMem, err); MSQ_ERRRTN(err);
    
    // TODD: Unrolling the linesearch loop.  We do a function and
    // gradient evaluation when beta = 1.  Otherwise, we end up
    // in the linesearch regime.  We expect that several
    // evaluations will be made, so we only do a function evaluation
    // and finish with a gradient evaluation.  When beta = 1, we also
    // check the gradient for stability.

    // TODD -- the Armijo linesearch is based on the objective function,
    //         so theoretically we only need to evaluate the objective
    //         function.  However, near a very accurate solution, say with
    //         the two norm of the gradient of the objective function less
    //         than 1e-5, the numerical error in the objective function
    //         calculation is enough that the Armijo linesearch will
    //         fail.  To correct this situation, the iterate is accepted
    //         when the norm of the gradient is also small.  If you need
    //         high accuracy and have a large mesh, talk with Todd about
    //         the numerical issues so that we can fix it.

    // TODD -- the Armijo linesearch here is only good for differentiable
    //         functions.  When projections are involved, you should change
    //	       to a form of the linesearch meant for nondifferentiable
    //         functions.

    pd.move_free_vertices_constrained(d, nv, beta, err); MSQ_ERRRTN(err);
    fn_bool = objFunc->compute_gradient(pd, grad, new_value, err); MSQ_ERRRTN(err);
    if ((fn_bool && (original_value - new_value >= -alpha*beta - epsilon)) ||
        (fn_bool && (length(grad, nv) < 100*convTol))) {
      // Armijo linesearch rules passed.
    }
    else {
      if (!fn_bool) {
	// Function undefined.  Use the higher decrease rate.
        beta *= beta0;
      }
      else {
	// Function defined, but not sufficient decrease
        // Use the lower decrease rate.
        beta *= beta1;
      }
      pd.set_to_vertices_memento(coordsMem, err); MSQ_ERRRTN(err);
      
      // Standard Armijo linesearch rules
 
      while (beta >= tol1) {
        // 6. Search along the direction
        //    (a) trial = x + beta*d
        pd.move_free_vertices_constrained(d, nv, beta, err); MSQ_ERRRTN(err);
        //    (b) function evaluation
        fn_bool = objFunc->evaluate(pd, new_value, err);  MSQ_ERRRTN(err);
        //    (c) check for sufficient decrease and stop
        if (!fn_bool) { 
	  // function not defined at trial point
          beta *= beta0;
        }
        else if (original_value - new_value >= -alpha*beta - epsilon ) {
          // iterate is acceptable.
          break; 
        }
        else {
          // iterate is not acceptable -- shrink beta
          beta *= beta1;
        }
        pd.set_to_vertices_memento(coordsMem, err); MSQ_ERRRTN(err);
      } 

      if (beta < tol1) {
	// assert(pd.set_to_vertices_memento called last)

        // TODD -- Lower limit on steplength reached.  Direction does not 
	//         appear to make sufficient progress decreasing the 
	//         objective function.  This can happen when you are 
        //         very near a solution due to numerical errors in 
	//         computing the objective function.  It can also happen 
        //         when the direction is not a descent direction and when
	//         you are projecting the iterates onto a surface.
	//
	// 	   The latter cases require the use of a linesearch on
	//         a gradient step.  If this linesearch terminate with
	//         insufficient decrease, then you are at a critical 
	//         point and should stop!
	//
	//	   The numerical errors with the objective function cannot
	//	   be overcome.  Eventually, the gradient step will
	//	   fail to compute a new direction and you will stop.

        MSQ_PRINT(1)("Sufficient decrease not obtained in linesearch; switching to gradient.\n");

	alpha = inner(grad, grad, nv); 	// compute norm squared of gradient
	if (alpha < 1) alpha = 1;	// take max with constant
	for (i = 0; i < nv; ++i) {
	  d[i] = -grad[i] / alpha; 	// compute scaled gradient
	}
	alpha = inner(grad, d, nv);  	// recompute alpha
	alpha *= sigma;                 // equal to one for large gradient
	beta = 1.0;

	// Standard Armijo linesearch rules
	while (beta >= tol2) {
	  // 6. Search along the direction
	  //    (a) trial = x + beta*d
	  pd.move_free_vertices_constrained(d, nv, beta, err); MSQ_ERRRTN(err);
	  //    (b) function evaluation
	  fn_bool = objFunc->evaluate(pd, new_value, err);  MSQ_ERRRTN(err);
	  //    (c) check for sufficient decrease and stop
	  if (!fn_bool) { 
	    // function not defined at trial point
	    beta *= beta0;
	  }
	  else if (original_value - new_value >= -alpha*beta - epsilon ) {
	    // iterate is acceptable.
	    break; 
	  }
	  else {
	    // iterate is not acceptable -- shrink beta
	    beta *= beta1;
	  }
	  pd.set_to_vertices_memento(coordsMem, err); MSQ_ERRRTN(err);
	} 

	if (beta < tol2) {
	  // assert(pd.set_to_vertices_memento called last)
	  
	  // TODD -- Lower limit on steplength reached.  Gradient does not 
	  //         appear to make sufficient progress decreasing the 
	  //         objective function.  This can happen when you are 
	  //         very near a solution due to numerical errors in 
	  //         computing the objective function.  Most likely you
	  //         are at a critical point for the problem.

	  MSQ_PRINT(1)("Sufficient decrease not obtained with gradient; critical point likely found.\n");
	  break;
	}
      }

      // Compute the gradient at the new point -- needed by termination check
      fn_bool = objFunc->compute_gradient(pd, grad, new_value, err); MSQ_ERRRTN(err);
    }

    // Prints out free vertices coordinates. 
    if (MSQ_DBG(3)) {
      MSQ_DBGOUT(3) << "  o Free vertices new coordinates: \n";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_ERRRTN(err);
      MsqFreeVertexIndexIterator ind(&pd, err); MSQ_ERRRTN(err);
      ind.reset();
      while (ind.next()) {
        MSQ_DBGOUT(3) << "\t\t\t" << toto1[ind.value()];
      }
    }

    // checks stopping criterion 
    term_crit->accumulate_inner( pd, new_value, grad, err ); MSQ_ERRRTN(err);
    term_crit->accumulate_patch( pd, err ); MSQ_ERRRTN(err);
  }
  MSQ_PRINT(2)("FINISHED\n");
}


void FeasibleNewton::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{

    //Michael::  Should the vertices memento be delete here???
  //  cout << "- Executing FeasibleNewton::iteration_complete()\n";
}
  
void FeasibleNewton::cleanup()
{
  delete coordsMem; coordsMem = NULL;
}
  

