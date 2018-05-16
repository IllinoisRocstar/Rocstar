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
  \file   SteepestDescent.cpp
  \brief  

  Implements the SteepestDescent class member functions.
  
  \author Thomas Leurent
  \date   2002-06-13
*/

#include "SteepestDescent.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqTimer.hpp"
#include "MsqDebug.hpp"

namespace Mesquite {

SteepestDescent::SteepestDescent(ObjectiveFunction* of) :
  VertexMover()
{
  objFunc=of;
  MsqError err;
  gradientLessThan=.01;
  maxIteration=6;
  this->set_name("SteepestDescent");
  set_patch_type(PatchData::GLOBAL_PATCH, err);
}  
  

void SteepestDescent::initialize(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
}

void SteepestDescent::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
  MSQ_FUNCTION_TIMER( "SteepestDescent::optimize_vertex_positions" );
    //PRINT_INFO("\no  Performing Steepest Descent optimization.\n");
  // Get the array of vertices of the patch. Free vertices are first.
  int num_vertices = pd.num_vertices();
  msq_std::vector<Vector3D> gradient(num_vertices), dk(num_vertices);
  int nb_iterations = 0;
  double norm=10e6;
  bool sd_bool=true;//bool for OF values
  double smallest_edge = 0.4; // TODO -- update -- used for step_size
  TerminationCriterion* term_crit=get_inner_termination_criterion();
  
  // does the steepest descent iteration until stopping is required.
  while ( (nb_iterations<maxIteration &&
          norm>gradientLessThan ) && !term_crit->terminate()) {
    
    ++nb_iterations;
    double original_value = 0.0;
      //get intial objective function value, original_value, and gradient
    objFunc->compute_gradient(pd, &gradient[0], original_value, 
			      err, gradient.size()); MSQ_ERRRTN(err);
    
    // Prints out free vertices coordinates. 
    if (MSQ_DBG(3)) {
      int num_free_vertices = pd.num_free_vertices(err); MSQ_ERRRTN(err);
      MSQ_DBGOUT(3) << "\n  o Free vertices ("<< num_free_vertices <<")original coordinates:\n ";
      MsqVertex* toto1 = pd.get_vertex_array(err); MSQ_ERRRTN(err);
      MsqFreeVertexIndexIterator ind1(&pd, err); MSQ_ERRRTN(err);
      ind1.reset();
      while (ind1.next()) {
        MSQ_DBGOUT(3) << "\t\t\t" << toto1[ind1.value()];
      }
    }
      
      // computes the gradient norm
    norm=0;
    for (int n=0; n<num_vertices; ++n) 
      norm += gradient[n] % gradient[n]; // dot product
    norm = sqrt(norm);
    MSQ_DBGOUT(3) << "  o  gradient norm: " << norm << msq_stdio::endl;
  
    if (norm <= gradientLessThan) {
      break;
    }

    // ******** Chooses the search direction ********
    // i.e., -gradient for the steepest descent
    for (int i=0; i<num_vertices; ++i)
      for (int j=0; j<3; ++j)
        dk[i][j] = -gradient[i][j] / norm;

    // ******* Improve Quality *******
    
      //set an error if initial patch is invalid.
    if(!sd_bool){
      MSQ_SETERR(err)("SteepestDescent passed invalid initial patch.",
                      MsqError::INVALID_ARG);
      return;
    }
    MSQ_DBGOUT(3) << "  o  original_value: " << original_value << msq_stdio::endl;
    
    double new_value = original_value+1;
    // reduces the step size until we get an improvement
    int nb_iter = 0;

    // saves the PatchData coordinates in a memento
    PatchDataVerticesMemento* pd_previous_coords;
    pd_previous_coords = pd.create_vertices_memento(err); 
    if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }
    // Loop to find a step size that improves quality
    double step_size = smallest_edge;
    while (new_value > original_value
           && nb_iter < 10 ) {
      nb_iter++;
        // change vertices coordinates in PatchData according to descent
        //direction.
      pd.move_free_vertices_constrained(&dk[0], dk.size(), step_size, err);
      if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }      
      // and evaluate the objective function with the new node positions.
      sd_bool=objFunc->evaluate(pd, new_value, err);  
      if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }
      if(!sd_bool){
        MSQ_SETERR(err)("SteepestDescent created invalid patch.",
                        MsqError::INVALID_MESH);
        delete pd_previous_coords;
        return;
      }
      MSQ_DBGOUT(3) << "    o  step_size: " << step_size << msq_stdio::endl; 
      MSQ_DBGOUT(3) << "    o  new_value: " << new_value << msq_stdio::endl;
      
      // if no improvement
      if (new_value > original_value) {
        // undoes node movement
        pd.set_to_vertices_memento( pd_previous_coords, err );  
        if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }
        // and reduces step size to try again.
        step_size /= 2;
      }
    }

    // Prints out free vertices coordinates. 
    if (MSQ_DBG(3)) {
      MSQ_DBGOUT(3) << "  o Free vertices new coordinates: \n";
      MsqVertex* toto1 = pd.get_vertex_array(err);  
      if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }
      MsqFreeVertexIndexIterator ind(&pd, err);  
      if (MSQ_CHKERR(err)) { delete pd_previous_coords; return; }
      ind.reset();
      while (ind.next()) {
        MSQ_DBGOUT(3) << "\t\t\t" << toto1[ind.value()];
      }
    }
    
    delete pd_previous_coords; // user manages the memento.
    
    term_crit->accumulate_inner( pd, err ); MSQ_ERRRTN(err); 
    
  }
}


void SteepestDescent::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing SteepestDescent::iteration_complete()\n";
}
  
void SteepestDescent::cleanup()
{
  //  cout << "- Executing SteepestDescent::iteration_end()\n";
}
  
} // namespace Mesquite
