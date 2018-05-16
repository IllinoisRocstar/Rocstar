/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
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

#include "CornerTag.hpp"
#include "MeshSet.hpp"

namespace Mesquite {

size_t CornerTagHandles::size( Mesh::TagType type )
{
  switch (type) {
    case Mesh::BYTE:   return 1;
    case Mesh::BOOL:   return sizeof(bool);
    case Mesh::INT:    return sizeof(int);
    case Mesh::DOUBLE: return sizeof(double);
    case Mesh::HANDLE: return sizeof(Mesh::EntityHandle);
    default:           return 0;
  }
}
  
int CornerTagHandles::num_corners( PatchData* pd, int elem_index )
{
  return pd->element_by_index(elem_index).vertex_count();
}


Mesh* CornerTagHandles::get_current_mesh( PatchData* pd )
{
  return pd->get_mesh_set() ? pd->get_mesh_set()->get_current_mesh() : 0;
}

TagHandle CornerTagHandles::get_handle( Mesh* mesh, unsigned corners, MsqError& err )
{
    // Resize vector as necessary
  if (corners >= cornerHandles.size())
  {
    cornerHandles.resize( corners+1, 0 );
  }
  
    // Get handle if we don't already have it
  if (!cornerHandles[corners])
  {
      // Construct tag name
    char numbuf[16];
    sprintf(numbuf, "%d", corners );
    msq_std::string name = tagName + numbuf;

      // Try to get existing handle
    cornerHandles[corners] = mesh->tag_get( name, err );  MSQ_ERRZERO(err);
      // If got handle, make sure type is correct
    if (cornerHandles[corners])
    {
      Mesh::TagType type;
      unsigned size;
      mesh->tag_properties( cornerHandles[corners], name, type, size, err );
      MSQ_ERRZERO(err); 

      if (type != tagType || size != corners*tagLen)
      {
        MSQ_SETERR(err)("Invalid tag type", MsqError::INVALID_ARG );
        return 0;
      }
    }
      // If didn't get handle, try to create it
    else 
    {
      cornerHandles[corners] = mesh->tag_create( name, tagType, tagLen * corners, 0, err );
      MSQ_ERRZERO(err);
    }
  }
  
  return cornerHandles[corners];
}

void CornerTagHandles::save_load_tags( bool load, PatchData* pd, 
                                       size_t elem_index, void* data, 
                                       size_t bytes, MsqError& err )
{
  Mesh* mesh = get_current_mesh( pd );
  if (!mesh) // this happens for some tests in the test suite.
    return;
    
  MsqMeshEntity* elem_array = pd->get_element_array( err ); 
  unsigned num_corners = elem_array[elem_index].vertex_count();
  const Mesh::ElementHandle handle = pd->get_element_handles_array()[elem_index];
  TagHandle tag = get_handle( mesh, num_corners, err ); MSQ_ERRRTN(err);
  if (load)
    mesh->tag_get_element_data( tag, 1, &handle, data, err ); 
  else
    mesh->tag_set_element_data( tag, 1, &handle, data, err );
  MSQ_CHKERR(err);
}

} // namespace Mesquite
