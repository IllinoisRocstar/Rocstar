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
  \file   MsqFreeVertexIndexIterator.hpp
  \brief    This file contains the MsqFreeVertexIndexIterator class

  \author Thomas Leurent
  \date   2002-01-17
*/

#ifndef MsqFreeVertexIndexIterator_hpp
#define MsqFreeVertexIndexIterator_hpp

#ifdef MSQ_USE_OLD_C_HEADERS
#  include <stdlib.h>
#  include <stddef.h>
#else
#  include <cstddef>
#  include <cstdlib>
#endif

#include "Mesquite.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"

namespace Mesquite
{
  class MsqError;

  /*! \class MsqFreeVertexIndexIterator
    \brief iterates over indexes of free vetices in a PatchData.

    A free vertex is defined as not having the MSQ_SOFT_FIXED and MSQ_HARD_FIXED
    flags activated.
    
    Use the iterator as follow:
    MsqFreeVertexIndexIterator ind(&patch_data,err);
    ind.reset();
    while (ind.next()) {
      cout << ind.value();
    }  .*/
  class MsqFreeVertexIndexIterator {
  public:
    MsqFreeVertexIndexIterator(PatchData *pd, MsqError &err) :
      iterOriginator(pd), iterCurrentIndex(0), initialState(true)
    { iterVertexArray = pd->get_vertex_array(err); }
    //! Resets the iterator. 
    //! The next call to next() will set the iterator on the first free vertex. 
    void reset() { initialState=true; iterCurrentIndex=0; }
    //! Increments the iterator. returns false if there is no more free vertex.
    inline bool next();
    //! Returns an index corresponding to a free vertex.
    size_t value() {return iterCurrentIndex;}
  private:
    PatchData* iterOriginator;
    size_t iterCurrentIndex;
    MsqVertex* iterVertexArray;
    bool initialState;
  };
  

  /*! \fn inline bool MsqFreeVertexIndexIterator::next() */
  inline bool MsqFreeVertexIndexIterator::next()
  {
    bool fixed=true;
    while ( fixed ) 
      {
        if ( initialState==true )  initialState=false;
        else  ++iterCurrentIndex;
        
        if ( iterCurrentIndex == iterOriginator->num_vertices() ) {
          return false; 
        }
        fixed = iterVertexArray[iterCurrentIndex].is_flag_set(MsqVertex::MSQ_SOFT_FIXED) ||
        iterVertexArray[iterCurrentIndex].is_flag_set(MsqVertex::MSQ_HARD_FIXED) ;
          //fixed = !(iterVertexArray[iterCurrentIndex].is_free_vertex());
      }
    return true;
  }
  


} // namespace

#endif //  MsqFreeVertexIndexIterator_hpp
