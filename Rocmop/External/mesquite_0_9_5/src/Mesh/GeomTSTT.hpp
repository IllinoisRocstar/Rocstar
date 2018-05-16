/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
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
 
    kraftche@cae.wisc.edu    
   
  ***************************************************************** */
/*!
  \file   GeomTSTT.cpp
  \brief  Mesquite::MeshDomain implemented on TSTT interfaces
  \author Jason Kraftcheck
  \date   2005-01-13
*/

#ifndef GEOM_TSTT_HPP
#define GEOM_TSTT_HPP

#include "Mesquite.hpp"
#include "MeshInterface.hpp"

namespace TSTTG {
  class Geometry;
}
namespace TSTTR {
  class Relate;
}
namespace TSTTM {
  class Mesh;
}

namespace Mesquite
{

    /**\brief A base class describing a Mesquite::MeshDomain implemented
     *        on top of the TSTT geometry and classification interfaces.
     *
     * A base class for Mesquite::MeshDomain implementations on top
     * of the TSTT geometry and classification interfaces and static
     * methods for constructing concrete implementations.  The concrete
     * implementations are not themselves in a header to avoid the
     * need to include all TSTT headers if this header is included, for
     * example indirectly through Mesquite_all_headers.hpp
     */
  class GeomTSTT : public Mesquite::MeshDomain
  {
    public:
    
      /**\brief Create MeshDommain from TSTT classification and geometry interfaces
        *
        * Create an instance of an implementation of MeshDomain that uses the 
        * TSTTR/Lasso interface to get a handle for the geometric entity 
        * associated with a mesh entity and the TSTT geometry interface to
        * evaluate the geometry.
        */
    static GeomTSTT* create( TSTTG::Geometry& geom, 
                             TSTTM::Mesh& mesh,
                             TSTTR::Relate& assoc,
                             MsqError& err );
    
      /**\brief Create a MeshDomain for a single geometric entity that uses
       *        the TSTT geometry interface for geometric evaluation.
       *
       * Create a TSTT geometry MeshDomain for a single geometric entity.
       * This implementation will be faster than the one that uses the
       * classification interface because it assumes all entities in the
       * mesh are in the interior of the single geometric entity specified.
       * This implemenation in intended to be used only in the case where 
       * the mesh of a single surface is being smoothed and the mesh 
       * vertices on the boundary of the surface are fixed.
       */
    static GeomTSTT* create( TSTTG::Geometry& geom,
                             void* geom_ent_handle,
                             MsqError& err );
  
    virtual ~GeomTSTT();
  };
}

#endif
