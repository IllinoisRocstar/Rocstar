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
  \file   SteepestDescent.hpp
  \brief  

  The SteepestDescent Class implements the steepest descent algorythm in
  order to move a free vertex to an optimal position given an
  ObjectiveFunction object and a QaulityMetric object.

  \author Thomas Leurent
  \date   2002-06-13
*/

#ifndef Mesquite_SteepestDescent_hpp 
#define Mesquite_SteepestDescent_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "ObjectiveFunction.hpp"


namespace Mesquite
{

  /*! \class SteepestDescent

      This is a very basic implementation of the steepest descent optimization algorythm.
      It works on patches of any size but the step size is hard-wired.
      Obvisouly, this is for testing purposed only. */ 
  class SteepestDescent : public VertexMover 
  {
  public:
    SteepestDescent(ObjectiveFunction* of);

    virtual ~SteepestDescent() { }

    /*! sets the maximum number of iteration of the steepest descent algorythm,
      i.e. the number of times we compute the gradient and try to move the nodes in the
      opposite direction. This is different from the number of passes over the mesh. */
    void set_maximum_iteration(int iter){
        maxIteration=iter;}

    /*! Sets a minimum value for the gradient. If the gradient is below that value,
      we stop iterating. */  
    void set_lower_gradient_bound(double gradc){
        gradientLessThan=gradc;}
    
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    double gradientLessThan;
    int maxIteration;
  };
  
}

#endif
