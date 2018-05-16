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

/*! \file MsqVertex.hpp
  \brief Mesquite's vertex object.

  \author Darryl Melander
  \author Thomas Leurent
*/
#ifndef MSQVERTEX_HPP
#define MSQVERTEX_HPP

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{
  class MeshSet;
  
    /*!
      \class MsqVertex
      \brief MsqVertex is the Mesquite object that stores information about
      the vertices in the mesh.

      This class has no virtual destructor for performance reasons.
      !!! Make sure NOT to delete a MsqVertex object from a pointer
          to Vector3D !!!
    */
  class MsqVertex : public Vector3D
   {
   public:
       //!Construct vertex using three doubles.
     MsqVertex(double x, double y, double z) 
       : Vector3D(x, y, z), vertexBitFlags(0)
       {}
     
       //!Construct vertex using Vector3D.
     MsqVertex(const Vector3D &vec) 
       : Vector3D(vec), vertexBitFlags(0)
       {}
     
       //!Construct default vertex with coordinates (0.0,0.0,0.0)
     MsqVertex() 
       : Vector3D(0,0,0), vertexBitFlags(0)
       {}

       //!Construct default vertex with coordinates (0.0,0.0,0.0)
     MsqVertex(const MsqVertex& rhs) 
       : Vector3D(rhs), vertexBitFlags(rhs.vertexBitFlags)
       {}

       //! Initializes with coordinates. Sets tag data/pointer to 0.
     MsqVertex& operator=(const Vector3D& rhs)
       { Vector3D::operator=(rhs);
         vertexBitFlags = 0;
         return *this; }
     
       // This allows for 8 flag bits.
       // I don't think we'll want more than that (yet).
     typedef unsigned char FlagMask;
     
       //! \enum FlagMaskID
       //!   Those are the available flags... currently only return
       //!   is_free.
       //!   Developers: The values used in that enum are used by a bitset,
       //!               so they have to be 2-based (2,4,8,16,32, ...)
     enum FlagMaskID
     {
       MSQ_NO_VTX_FLAG = 0,
       MSQ_ALGO_FLAG0 = 1<<0, //!< vertex is "free"
       MSQ_SOFT_FIXED = 1<<1, //!< vertex is fixed. This flag can be set on and off. 
       MSQ_HARD_FIXED = 1<<2, //!< vertex is always fixed. This can only be set on and never off.
       MSQ_COORDS_CHANGED = 1<<3,
       MSQ_FLAG_3 =     1<<4,
       MSQ_FLAG_4 =     1<<5,
       MSQ_ALGO_FLAG1 = 1<<6, //!< free bit, to be used by algorithm if needed.
       MSQ_ALGO_FLAG2 = 1<<7 //!< free bit, to be used by algorithm if needed. 
     };
       //!Returns true if vertex is ``free''.
     bool is_free_vertex() const
       { return (vertexBitFlags & (MSQ_SOFT_FIXED|MSQ_HARD_FIXED)) == 0; }
     
     void set_soft_fixed_flag()
       { vertexBitFlags|=MSQ_SOFT_FIXED; }
     
     void remove_soft_fixed_flag()
       { vertexBitFlags &= (~MSQ_SOFT_FIXED); }
     
     void set_hard_fixed_flag()
       { vertexBitFlags|=MSQ_HARD_FIXED; }
     
     void set_vertex_flag(FlagMaskID flag)
       { vertexBitFlags|=flag; }
     
     void remove_vertex_flag(FlagMaskID flag)
       { vertexBitFlags &= (~flag); }
     
     bool is_flag_set(FlagMaskID flag) const
       { return (vertexBitFlags & flag) != 0; }
     
   private:
     FlagMask vertexBitFlags;

     friend class Mesquite::MeshSet;
   };

} //namespace
  

#endif // MsqVertex_hpp
