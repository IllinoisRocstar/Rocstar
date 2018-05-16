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

#ifndef MSQ_CORNER_TAG_HPP
#define MSQ_CORNER_TAG_HPP

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
#endif
#include <string>
#include <assert.h>

#include "MeshInterface.hpp"
#include "MsqError.hpp"

namespace Mesquite {

class PatchData;

/** \brief Utility class to manage tag handles for corner tags 
 *
 * Tags on corners are saved as an array of values on the corresponding
 * element.  As tags must be of a constant size, this necessitates using
 * different tags for elements with different numbers of corners.  This
 * class manages the set of all tag handles for a given corner tag, where
 * the specific tag handle can be retrieved by specifying the number of
 * corners in the element.
 */ 
class CornerTagHandles
{
public:
  CornerTagHandles( const char* tag_name, 
                    Mesh::TagType type, 
                    unsigned tag_len ) 
    : tagName( tag_name ), 
      tagType(type), 
      tagLen(tag_len) 
    {}
  
  /** Get the tag handle for storing this tag type with the specified 
   *  number of values (corners).
   */  
  TagHandle get_handle( Mesh* mesh, unsigned num_corners, MsqError& err );
  
  void save_load_tags( bool load, PatchData* pd, size_t elem_index, void* data, size_t tag_byes, MsqError& err );

  static Mesh* get_current_mesh( PatchData* pd );
  
  static size_t size( Mesh::TagType type );
  
  static int num_corners( PatchData* pd, int elem_index );

private:

  const msq_std::string tagName;
  const Mesh::TagType tagType;
  const unsigned tagLen;

  msq_std::vector<TagHandle> cornerHandles;
};


/** \brief A class for caching and managing Tags on element corners.
 *
 * This class provides:
 * - An abstraction of the mechanism for storing tags on element corners.
 * - Management of accessing tag data in the Mesh
 * - Caching of tag data for elements in a patch
 * - Converting between tag data types and Mesquite objects
 *
 * Due to limitations in accessing tag data (and for more efficient
 * access to tag data) this class assumes that either all or none
 * of the corner tags have been defined for the elements in a patch.
 *
 * When creating new tags, first call \ref allocate_new_tags to allocate
 * a local cache of the tag data for all elements in the patch.  Then 
 * use \ref get_element_corner_tags to retrieve the allocated space for
 * the corners of each element.  Calling \ref get_element_corner_tags
 * without first calling \ref allocate_new_tags will result in an 
 * attempt to read the tag data for all elements from the Mesh instance.
 * 
 * Tag data is not saved to the mesh unless \ref save_tag_data is
 * explicitly called.
 */
template <typename T> class CornerTag 
{
  public:
  
      /**\brief Initialize
       *\param mesh  A pointer to the Mesh instance.
       *\param name  The tag name.
       *\param type  The native tag type.  If the type "T" can be
       *             cast to one or an array of some native type,
       *             specify that type here.  Otherwise use \ref Mesh::BYTE
       */
    CornerTag( const char* name, Mesh::TagType type = Mesh::BYTE);
    
    ~CornerTag();
    
      /** Clear cached data.  Any changes will be lost if
       *  \ref save_tag_data has not been called.
       */
    inline void clear( );

      /** Get a pointer to the array of all corner tag values for
       *  a given element.
       *\param pd       The PatchData
       *\param elem_idx The element, specified as it's index in the PatchData.
       */
    inline const T* get_element_corner_tags( PatchData* pd, int elem_idx, MsqError& err );
    
    inline void set_element_corner_tags( PatchData* pd, int elem_idx, const T* data, MsqError& err );
      
  private:

    CornerTagHandles tagHandles;
    msq_std::vector<T*> tagData; //< Cached tag data for all elems of patch
};


template <typename T>
CornerTag<T>::CornerTag( const char* tag_name, Mesh::TagType type )
 : tagHandles( tag_name, type, sizeof(T) / CornerTagHandles::size(type) )
{
    // We're going to assume this for all tag read/writes,
    // so make sure its true now.
  assert( (sizeof(T) % CornerTagHandles::size(type)) == 0 );
}

template <typename T>
CornerTag<T>::~CornerTag()
{
  clear();
}

template <typename T>
void CornerTag<T>::clear()
{
  for (typename msq_std::vector<T*>::iterator i = tagData.begin(); i != tagData.end(); ++i)
    delete [] *i;
  tagData.clear();
}

template <typename T>
const T* CornerTag<T>::get_element_corner_tags( PatchData* pd, 
                                                int elem_index, 
                                                MsqError& err )
{
  if (tagData.size() <= (unsigned)elem_index)
    tagData.resize( elem_index+1, 0 );
  
  int num_corners = tagHandles.num_corners( pd, elem_index );
  if (!tagData[elem_index]) {
    tagData[elem_index] = new T[num_corners];
    tagHandles.save_load_tags( true, pd, elem_index, tagData[elem_index], sizeof(T), err );
    MSQ_ERRZERO(err);
  }
  
  return tagData[elem_index];
}

template <typename T>
void CornerTag<T>::set_element_corner_tags( PatchData* pd, 
                                            int elem_index, 
                                            const T* data,
                                            MsqError& err )
{
  if (tagData.size() <= (unsigned)elem_index)
    tagData.resize( elem_index+1, 0 );
  
  int num_corners = tagHandles.num_corners( pd, elem_index );
  if (!tagData[elem_index]) 
    tagData[elem_index] = new T[num_corners];
  memcpy( tagData[elem_index], data, num_corners * sizeof(T) );
    
  tagHandles.save_load_tags( false, pd, elem_index, (void*)data, sizeof(T), err );
  MSQ_CHKERR(err);
}
  

} // namespace Mesquite

#endif
