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

#ifndef MESQUITE_PATCH_DATA_MEM_HPP
#define MESQUITE_PATCH_DATA_MEM_HPP

#include <stdlib.h>
#include <assert.h>

namespace Mesquite {

template <typename X> class PatchDataMem
{
  enum {
    SHRINK_PERCENT = 10  /**< Shink memory when used portion is less than this percent */
  };
  
  private:
  
    size_t arrayLength;
    size_t activeSize;
    X* arrayData;
    
    inline void resize_storage( size_t size );
  
  public:
    
    inline PatchDataMem();
    inline PatchDataMem( size_t size );
    inline PatchDataMem( const PatchDataMem<X>& );
    
    inline ~PatchDataMem();
  
    inline size_t size() const;
    inline bool empty() const;

    inline const X& operator[]( size_t index ) const;
    inline X& operator[]( size_t index ) ;
    
    inline const X* to_array() const;
    inline X* to_array();
    
    inline void resize( size_t new_size );
    inline void clear();
    
    inline PatchDataMem<X>& operator=( const PatchDataMem<X>& );
    
    typedef X* iterator;
    typedef const X* const_iterator;
    inline iterator begin();
    inline const_iterator begin() const;
    inline iterator end();
    inline const_iterator end() const;
};

template <typename X> 
PatchDataMem<X>::PatchDataMem()
  : arrayLength(0),
    activeSize(0),
    arrayData(0)
  {}

template <typename X> 
PatchDataMem<X>::PatchDataMem( size_t size )
  : arrayLength(0),
    activeSize(0),
    arrayData(0)
  { resize( size ); }
  
template <typename X> 
PatchDataMem<X>::PatchDataMem( const PatchDataMem<X>& other )
  : arrayLength(0),
    activeSize(0),
    arrayData(0)
  { operator=(other); }

template <typename X>
PatchDataMem<X>::~PatchDataMem()
  { clear(); free( arrayData ); }

template <typename X> 
size_t PatchDataMem<X>::size() const
  { return activeSize; }

template <typename X> 
const X* PatchDataMem<X>::to_array() const
  { return arrayData; }

template <typename X> 
X* PatchDataMem<X>::to_array() 
  { return arrayData; }

template <typename X> 
const X& PatchDataMem<X>::operator[]( size_t index ) const 
  { return arrayData[index]; }

template <typename X> 
X& PatchDataMem<X>::operator[]( size_t index ) 
  { return arrayData[index]; }

template <typename X> 
bool PatchDataMem<X>::empty() const
  { return !activeSize; }

template <typename X> 
typename PatchDataMem<X>::iterator PatchDataMem<X>::begin()
  { return arrayData; }

template <typename X> 
typename PatchDataMem<X>::const_iterator PatchDataMem<X>::begin() const
  { return arrayData; }

template <typename X> 
typename PatchDataMem<X>::iterator PatchDataMem<X>::end()
  { return arrayData + activeSize; }

template <typename X> 
typename PatchDataMem<X>::const_iterator PatchDataMem<X>::end() const
  { return arrayData + activeSize; }

template <typename X> 
void PatchDataMem<X>::resize_storage( size_t new_size )
{
  assert( new_size >= activeSize );
  X* new_array = (X*)malloc( sizeof(X) * new_size );
  for (size_t i = 0; i < activeSize; ++i)
  {
    //new_array[i].X( arrayData[i] );
    new (new_array + i) X( arrayData[i] );
    arrayData[i].~X();
  }
  free( arrayData );
  arrayData = new_array;
  arrayLength = new_size;
}

template <typename X> 
void PatchDataMem<X>::resize( size_t new_size )
{
  if (new_size > activeSize)
  {
    if (new_size > arrayLength)
      resize_storage( new_size );
    
    for (size_t i = activeSize; i < new_size; ++i)
    {
      //arrayData[i].X();
      new (arrayData + i) X();
    }
    activeSize = new_size;
  }
  else if (new_size < activeSize)
  {
    for (size_t i = new_size; i < activeSize; ++i)
      arrayData[i].~X();
    activeSize = new_size;
    
    if (activeSize && ((SHRINK_PERCENT*activeSize) < arrayLength))
      resize_storage( new_size );
  }
}

template <typename X> 
void PatchDataMem<X>::clear()
{
  for (size_t i = 0; i < activeSize; ++i)
    arrayData[i].~X();
  activeSize = 0;
}

template <typename X> 
PatchDataMem<X>& PatchDataMem<X>::operator=( const PatchDataMem<X>& other )
{
  clear();
  
  if (other.activeSize)
  {
    if (arrayLength < other.activeSize ||
        (SHRINK_PERCENT * other.activeSize) < arrayLength)
      resize_storage( other.activeSize );
  
    activeSize = other.activeSize;
    for (size_t i = 0; i < activeSize; ++i)
    {
      //arrayData[i].X( other.arrayData[i] );
      new (arrayData + i) X( other.arrayData[i] );
    }
  }
  
  return *this;
}  

} // namespace Mesquite

#endif
    
  
