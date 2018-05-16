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
  \file   Randomize.cpp
  \brief  

  The Randomize Class is the concrete class that randomizes
  the vertex positions.          

  \author Michael Brewer
  \date   2002-10-27
*/


#include "Randomize.hpp"

using namespace Mesquite;


Randomize::Randomize() 
{
  this->set_name("Randomize");
  mPercent=.05;
}  

Randomize::Randomize(double percent) 
{
  this->set_name("Randomize");
  mPercent=percent;
}  
  
void Randomize::initialize(PatchData &/*pd*/, MsqError &err)
{
  this->set_patch_type(PatchData::ELEMENTS_ON_VERTEX_PATCH, err, 1);
}

void Randomize::initialize_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}

void Randomize::optimize_vertex_positions(PatchData &pd, 
                                                MsqError &err)
{
    //cout << "- Executing Randomize::optimize_vertex_position()\n";

  int num_local_vertices = pd.num_vertices();
  // gets the array of coordinates for the patch and print it 
  MsqVertex *patch_coords = pd.get_vertex_array(err); MSQ_ERRRTN(err);
  // does the randomize smooth
  MsqFreeVertexIndexIterator free_iter(&pd, err); MSQ_ERRRTN(err);
  free_iter.reset();
  free_iter.next();
    //find the free vertex.
  int m=free_iter.value();
  randomize_vertex(pd, num_local_vertices,
                   patch_coords[m], err);  MSQ_ERRRTN(err);
  pd.snap_vertex_to_domain(m,err); MSQ_ERRRTN(err);
}
  
void Randomize::terminate_mesh_iteration(PatchData &/*pd*/, MsqError &/*err*/)
{
  //  cout << "- Executing Randomize::iteration_complete()\n";
}
  
void Randomize::cleanup()
{
  //  cout << "- Executing Randomize::iteration_end()\n";
}
  

