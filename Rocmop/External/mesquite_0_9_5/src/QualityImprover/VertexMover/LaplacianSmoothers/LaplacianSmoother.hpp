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
  \file   LaplacianSmoother.hpp
  \brief  

  The LaplacianSmoother Class implements the Laplacian smoothing
  for a patch with one free vertex. 

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef Mesquite_LaplacianSmoother_hpp 
#define Mesquite_LaplacianSmoother_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqFreeVertexIndexIterator.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif

namespace Mesquite
{

  /*! \class LaplacianSmoother
    Moves free center vertex to the average of the neighboring vertices.
   */  
  class LaplacianSmoother : public VertexMover 
  {
  public:
    LaplacianSmoother(MsqError &err);
    ~LaplacianSmoother();
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();

  private:
    QualityMetric* edgeQM;

  };

  
  inline void centroid_smooth_mesh(PatchData &pd, size_t num_adj_vtx,
                                   msq_std::vector<size_t> adj_vtx_ind,
                                   size_t free_ind,
                                   size_t dimension, MsqError &err)
  {
    MsqVertex* verts=pd.get_vertex_array(err);
    msq_std::vector<size_t>::iterator iter;
    
    size_t j;
    double scale_val=1.0;
    if (num_adj_vtx==0) {
      MSQ_SETERR(err)("Number of incident vertices is zero",MsqError::INVALID_ARG);
      return;
    }
    else
      scale_val=1.0/((double) num_adj_vtx);
    double avg[3];
      //loop over the two or three dimensions
    for(j=0;j<dimension;++j) {
        //set the iterator to the beginning ob adj_vtx_ind
      iter=adj_vtx_ind.begin();
      avg[j] = 0.;
      while(iter != adj_vtx_ind.end()){
        avg[j]+=verts[*iter][j];
        ++iter;
      }
        //divide the average by the number of adj. verts
      verts[free_ind][j] = avg[j]*scale_val;
    }

    return;
  }

  
}

#endif
