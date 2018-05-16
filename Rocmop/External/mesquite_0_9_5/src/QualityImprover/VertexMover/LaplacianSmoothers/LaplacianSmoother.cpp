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
  \file   LaplacianSmoother.cpp
  \brief  

  The LaplacianSmoother Class is the concrete class
  that performs Laplacian Smoothing

  \author Thomas Leurent
  \date   2002-01-17
*/

#include "LaplacianSmoother.hpp"
#include "LPtoPTemplate.hpp"
#include "EdgeLengthQualityMetric.hpp"


namespace Mesquite {


LaplacianSmoother::LaplacianSmoother(MsqError &err) 
{
  this->set_name("LaplacianSmoother");
  
  set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err,1,1);MSQ_ERRRTN(err);

  edgeQM = new EdgeLengthQualityMetric;
  edgeQM->set_averaging_method(QualityMetric::RMS,err);MSQ_ERRRTN(err);
  objFunc = new LPtoPTemplate(edgeQM, 2, err);MSQ_ERRRTN(err);
  
}  

LaplacianSmoother::~LaplacianSmoother() 
{
  delete edgeQM;
  delete objFunc;
}    
  

void LaplacianSmoother::initialize(PatchData& /*pd*/, MsqError& /*err*/)
{
 
}


void LaplacianSmoother::initialize_mesh_iteration(PatchData &/*pd*/,
                                                  MsqError &/*err*/)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}


/*! \todo Michael:  optimize_vertex_position is probably not implemented
  in an optimal way.  We used to use all of the vertices in
  the patch as 'adjacent' vertices.  Now we call get_adjacent_vertex_indices.
  We could use a VERTICES_ON_VERTEX type of patch or a global patch?
*/
void LaplacianSmoother::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //default the laplacian smoother to 3 even for 2-d elements.
    //int dim = get_mesh_set()->space_dim();
  size_t dim = 3;
  
  
  // does the Laplacian smoothing
  MsqFreeVertexIndexIterator free_iter(&pd, err);  MSQ_ERRRTN(err);
  free_iter.reset();
  free_iter.next();
    //m is the free vertex.
  size_t m=free_iter.value();
  msq_std::vector<size_t> vert_indices;
  vert_indices.reserve(25);
    //get vertices adjacent to vertex m
  pd.get_adjacent_vertex_indices(m,vert_indices,err);  MSQ_ERRRTN(err);
    //move vertex m
  centroid_smooth_mesh(pd, vert_indices.size(), vert_indices,
                       m, dim, err); MSQ_ERRRTN(err);
    //snap vertex m to domain
  pd.snap_vertex_to_domain(m,err);
  
}
  
void LaplacianSmoother::terminate_mesh_iteration(PatchData &/*pd*/,
                                                 MsqError &/*err*/)
{
  //  cout << "- Executing LaplacianSmoother::iteration_complete()\n";
}
  
void LaplacianSmoother::cleanup()
{
  //  cout << "- Executing LaplacianSmoother::iteration_end()\n";
}
  
} // namespace Mesquite

