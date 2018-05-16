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
  \file   Randomize.hpp
  \brief  

  The Randomize Class implements the Randomize Vertex Mover
  for a patch with one free vertex. 

  \author Michael Brewer      
  \date   2002-10-27
*/

#ifndef Mesquite_Randomize_hpp 
#define Mesquite_Randomize_hpp

#include "Mesquite.hpp"
#include "VertexMover.hpp"
#include "MsqFreeVertexIndexIterator.hpp"
#include "MsqDebug.hpp"
#include <math.h>
namespace Mesquite
{

  /*! \class Randomize
   \brief Randomly perftubs the (un-culled) vertices.
  */ 
  class Randomize : public VertexMover 
  {
  public:
      //!Constructor defaulting mPercent to .05.
    Randomize();
      //!Constructor allowing user to set mPercent
    Randomize(double percent);

    virtual ~Randomize() { }
    
  protected:
    virtual void initialize(PatchData &pd, MsqError &err);
    virtual void optimize_vertex_positions(PatchData &pd,
                                         MsqError &err);
    virtual void initialize_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void terminate_mesh_iteration(PatchData &pd, MsqError &err);
    virtual void cleanup();
  private:
      //! \param The percentage of the scale factor each vertex will be moved.
    double mPercent;
      /*!Function calculates a scale factor for the patch, then moves
        the incident vertex randomly in each of the three coordinate
        directions (relative to the scale factor multiplied by mPercent).
      */
    inline void randomize_vertex(PatchData &pd,size_t num_vtx,
                                 MsqVertex &free_vtx,
                                 MsqError &err);
  };

  
    //!Perturbs the free vertex randomly.
  inline void Randomize::randomize_vertex(PatchData &pd, size_t num_vtx,
                                          MsqVertex &free_vtx,
                                          MsqError &err)
  {
    size_t i;
    short j;
    MsqVertex* verts = pd.get_vertex_array(err);  MSQ_ERRRTN(err);
      //a scale w.r.t. the patch size
    double scale_factor=0.0;
      //a "random" number between -1 and 1
    double rand_double=0.0;
      //a "random" int
    int rand_int=0;
    if (num_vtx<=1){
      MSQ_PRINT(1)("WARNING: Number of incident vertex is zero.  Returning.\n");
      return;
    }

    size_t free_ind = pd.get_vertex_index(&(free_vtx));
    
    
    for (i=0;i<num_vtx;++i){
      if(i != free_ind)
        scale_factor+=(verts[i]-free_vtx).length();
    }
    scale_factor/=( (double) num_vtx - 1.0 );    
    for (j=0;j<3;++j){
      rand_int = rand();
        //number between 0 and 1000
      rand_int = rand_int%1000;
        //number between -1 and 1
      rand_double = (((double) rand_int)/500.0)-1.0;
      free_vtx[j] += scale_factor*rand_double*mPercent;
    }
    
    return;
  }

  
}

#endif
