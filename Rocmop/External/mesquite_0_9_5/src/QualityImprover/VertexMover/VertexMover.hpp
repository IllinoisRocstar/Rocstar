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
  \file   VertexMover.hpp
  \brief  

  The VertexMover Class is the base class for all the smoothing algorythms 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_VertexMover_hpp 
#define Mesquite_VertexMover_hpp


#include "Mesquite.hpp"
#include "QualityImprover.hpp"
#include "PatchData.hpp"
#include "ObjectiveFunction.hpp"

namespace Mesquite
{

  /*! \class VertexMover
    Base class for all Vertex Movers.
   */  
  class VertexMover : public QualityImprover 
  {
  protected:
    VertexMover();
  public:
    // virtual destructor ensures use of polymorphism during destruction
    virtual ~VertexMover() { };
    
    virtual double loop_over_mesh(MeshSet &ms, MsqError &err);

  protected:

    virtual void initialize(PatchData &pd, MsqError &err) = 0;
    virtual void cleanup() = 0;
    virtual void optimize_vertex_positions(PatchData &pd, 
                                           MsqError &err) = 0; // modifies the PatchData object

    virtual void initialize_mesh_iteration(PatchData &pd, 
                                         MsqError &err) = 0;
    virtual void terminate_mesh_iteration(PatchData &, 
                                         MsqError &err) = 0;

      //!CHECK FEASIBLE IS NOT YET IMPLEMENTED.
    size_t check_feasible(PatchData &pd, MsqError &err);
    
    ObjectiveFunction* objFunc;
  };

  
/*!
  Takes a PatchData object (by reference) and returns whether the
  patch is within the feasible region, 0, or outside the region, 1.
*/
  inline size_t VertexMover::check_feasible(PatchData &pd, MsqError &err)
  {
    MsqMeshEntity* elems=pd.get_element_array(err);
    size_t num_elements=pd.num_elements();
    msq_std::vector<Vector3D> sample_points;
    Vector3D jacobian_vectors[3];
    short num_jacobian_vectors;
    size_t i =0;
    for(i=0;i<num_elements;++i){
      elems[i].get_sample_points(QualityMetric::ELEMENT_VERTICES,sample_points,err);
      msq_std::vector<Vector3D>::iterator iter=sample_points.begin();
      while(iter!=sample_points.end()){
        elems[i].compute_weighted_jacobian(pd, (*iter),
                                           jacobian_vectors,
                                           num_jacobian_vectors, err);
        if(num_jacobian_vectors==2){
            //2-d not yet implemented
        }
        else if(num_jacobian_vectors==3){
          if(jacobian_vectors[0]%(jacobian_vectors[1]*
                                   jacobian_vectors[2])<=0.0){
            return 1;
          }
        }
        ++iter;
      }
    }
    
    return 0;
  }
    

} // namespace
#endif // Mesquite_VertexMover_hpp
