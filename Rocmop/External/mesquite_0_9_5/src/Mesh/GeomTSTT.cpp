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
  \brief  


  \author Jason Kraftcheck
  \date   2005-01-13
*/

#include <sidl_cxx.hh>
#include "TSTTB.hh"
#include "TSTTM.hh"
#include "TSTTR.hh"
#include "TSTTG.hh"
#include "GeomTSTT.hpp"
#include "MsqDebug.hpp"
#include "MsqVertex.hpp"
#include "MsqError.hpp"
#include "TSTTUtil.hpp"

namespace Mesquite
{

/**\brief Common code for specific implementations of MeshDomain on TSTT interfaces.
 *
 * This class contains the common functionality used by concrete implementations
 * of MeshDomain on the TSTT geometry interface.
 */
class GeomTSTTCommon
{
public:

    /**\param geom The TSTT geometry interface implementation to query */
  GeomTSTTCommon( TSTTG::Geometry& geom ) throw (TSTTB::Error );
  
  virtual ~GeomTSTTCommon();

    /** Evaluate the closest point to the input position on the specified
     *  geometric entity and return the result in the passed position 
     *  argument (move the passed position onto the geometry.)
     */
  void move_to( void* geom_handle, Vector3D& coord ) const throw( TSTTB::Error );
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  void normal ( void* geom_handle, Vector3D& coord ) const throw( TSTTB::Error );
  
    /** Given a geometric entity and a position, evaluate the normal 
     *  on the geometric entity at the closest point on that entity
     *  to the input position, and pass back the result in the input
     *  coord vector.
     */
  void normal ( void* geom_handle, Vector3D coords[], unsigned count ) const 
    throw( TSTTB::Error );
    
    /** Given a geometric entity and a position, get point on 
     *  the geometric entity closest to the input position, and
     *  the surface normal at that position.
     */
  void closest_and_normal( void* geom_handle,
                           const Vector3D& position,
                           Vector3D& closest, 
                           Vector3D& normal ) const 
                         throw (TSTTB::Error);

    /** TSTT geometry interface implementation to query */
  mutable TSTTG::Shape  geomIface;
  
    /** Temporary storage for geometry entity handles */
  mutable sidl::array<void*> geomHandles;
    /** Temporary storate for input and output vectors */
  mutable sidl::array<double> positionsIn, positionsOut, normalsOut;
};


/* General MeshDomain on TSTT implementation */

class DomainTSTT : public GeomTSTT, protected GeomTSTTCommon
{
public:

  DomainTSTT( TSTTG::Geometry& geom,
              TSTTM::Mesh& mesh, 
              TSTTR::Relate& relate ) throw( TSTTB::Error );

  virtual ~DomainTSTT();
  
  DomainHint hint() const;

  void snap_to( Mesh::EntityHandle entity_handle,
                Vector3D& coordinat ) const;

  void normal_at( Mesh::EntityHandle entity_handle,
                  Vector3D& coordinate ) const;
  
  void normal_at( Mesh::EntityHandle handle,
                  Vector3D coordinates[],
                  unsigned count,
                  MsqError& err ) const;

  void closest_point( Mesh::EntityHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;
                      
protected:

    /** Get geometric entity owning a mesh entity */
  void* geom_from_mesh( void* handle ) const throw( TSTTB::Error );

private:

    /** TSTT interface implementation used to get dimension of geometric entities */
  mutable TSTTG::Topology topoIface;
    /** TSTT interface used to query mesh entity properties */
  mutable TSTTM::Mesh   meshIface;
    /** TSTT interface implementation for mesh->geometry association */
  mutable TSTTR::Relate relateIface;
  
    /** Temporary storage for mesh entity handles */
  mutable sidl::array<void*> oneMeshHandle;
    /** Temporary storage for geometry entity type */
  mutable sidl::array<TSTTG::GentityType> oneTypeOut;
};


/* Specialized one-geometry-entity MeshDomain on TSTT implementation */

class GeomEntTSTT : public GeomTSTT, protected GeomTSTTCommon
{
public:

  GeomEntTSTT( TSTTG::Geometry& geom,
               void* geom_ent_handle ) throw( TSTTB::Error );

  virtual ~GeomEntTSTT();
  
  DomainHint hint() const;

  void snap_to( Mesh::EntityHandle entity_handle,
                Vector3D& coordinat ) const;

  void normal_at( Mesh::EntityHandle entity_handle,
                  Vector3D& coordinate ) const;

  
  void normal_at( Mesh::EntityHandle handle,
                  Vector3D coordinates[],
                  unsigned count,
                  MsqError& err ) const;

  void closest_point( Mesh::EntityHandle handle,
                      const Vector3D& position,
                      Vector3D& closest,
                      Vector3D& normal,
                      MsqError& err ) const;
private:
  
    /** A handle for the geometry entity to evaluate */
  void* geomEntHandle;
};




/*****************  GeomTSTT base class methods *********************/

GeomTSTT* GeomTSTT::create( TSTTG::Geometry& geom,
                            TSTTM::Mesh& mesh,
                            TSTTR::Relate& relate,
                            MsqError& err )
{
  try {
    return new DomainTSTT( geom, mesh, relate );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
  return 0;
}

GeomTSTT* GeomTSTT::create( TSTTG::Geometry& geom,
                            void* geom_ent_handle,
                            MsqError& err )
{
  try {
    return new GeomEntTSTT( geom, geom_ent_handle );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
  return 0;
}
  

GeomTSTT::~GeomTSTT() {}




/***************** DomainTSTT class methods *********************/

DomainTSTT::DomainTSTT( TSTTG::Geometry& geom,
                            TSTTM::Mesh& mesh,
                            TSTTR::Relate& relate ) 
                            throw ( TSTTB::Error )
  : GeomTSTTCommon( geom ), 
    topoIface( geom ),
    meshIface(mesh), 
    relateIface( relate ),
    oneMeshHandle( alloc_sidl_vector<void*>(1) ),
    oneTypeOut( alloc_sidl_vector<TSTTG::GentityType>(1) )
{
}

DomainTSTT::~DomainTSTT() {}


DomainHint DomainTSTT::hint() const 
  { return NO_DOMAIN_HINT; }

void DomainTSTT::snap_to( Mesh::EntityHandle handle,
                            Vector3D& coordinate ) const
{
  try {
    void* geom = geom_from_mesh( (void*)handle ); 
    if (geom)
      move_to( geom, coordinate );
  }
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
  }
}

void DomainTSTT::normal_at( Mesh::EntityHandle handle,
                              Vector3D& coordinate ) const
{
  try {
    void* geom = geom_from_mesh( (void*)handle );
    if (geom)
      normal( geom, coordinate );
  }
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
  }
}

void DomainTSTT::normal_at( Mesh::EntityHandle handle,
                            Vector3D coordinates[],
                            unsigned count,
                            MsqError& err ) const
{
  try {
    void* geom = geom_from_mesh( (void*)handle );
    if (!geom) {
      MSQ_SETERR(err)(MsqError::INVALID_ARG);
      return;
    }
    
    normal( geom, coordinates, count );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
}

void DomainTSTT::closest_point( Mesh::EntityHandle handle,
                                const Vector3D& position,
                                Vector3D& closest,
                                Vector3D& normal,
                                MsqError& err ) const
{
  try {
    void* geom = geom_from_mesh( (void*)handle );
    if (!geom) {
      MSQ_SETERR(err)(MsqError::INVALID_ARG);
      return;
    }
    
    closest_and_normal( geom, position, closest, normal );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
}


void* DomainTSTT::geom_from_mesh( void* mesh_ent_handle ) const
                                    throw ( TSTTB::Error )
{
  int junk;
  oneMeshHandle.set( 0, mesh_ent_handle );
  relateIface.getMeshRelatedEntities( &geomIface,
                                      &meshIface,
                                      oneMeshHandle, 1,
                                      geomHandles, junk );
  // get dimension
  junk = 1;
  topoIface.gentityGetType( geomHandles, 1, oneTypeOut, junk );
                                      
  return oneTypeOut[0] == TSTTG::GentityType_GREGION ? 0 : geomHandles.get(0);
}




/***************** GeomEntTSTT class methods *********************/

GeomEntTSTT::GeomEntTSTT( TSTTG::Geometry& geom, void* geom_ent_handle ) 
                            throw ( TSTTB::Error )
  : GeomTSTTCommon( geom ), 
    geomEntHandle( geom_ent_handle )
{
}

GeomEntTSTT::~GeomEntTSTT() {}

DomainHint GeomEntTSTT::hint() const 
  { return SMOOTH_DOMAIN; }


void GeomEntTSTT::snap_to( Mesh::EntityHandle handle,
                            Vector3D& coordinate ) const
{
  try {
    move_to( geomEntHandle, coordinate );
  }
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
  }
}

void GeomEntTSTT::normal_at( Mesh::EntityHandle handle,
                              Vector3D& coordinate ) const
{
  try {
    normal( geomEntHandle, coordinate );
  }
  catch (TSTTB::Error& tstt_err ) {
    process_tstt_error(tstt_err);
  }
}


void GeomEntTSTT::normal_at( Mesh::EntityHandle handle,
                             Vector3D coordinates[],
                             unsigned count,
                             MsqError& err ) const
{
  try {
    normal( geomEntHandle, coordinates, count );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
}

void GeomEntTSTT::closest_point( Mesh::EntityHandle handle,
                                 const Vector3D& position,
                                 Vector3D& closest,
                                 Vector3D& normal,
                                 MsqError& err ) const
{
  try {
    closest_and_normal( geomEntHandle, position, closest, normal );
  }
  catch (TSTTB::Error& tstt_err ) {
    MSQ_SETERR(err)(process_tstt_error(tstt_err),MsqError::INTERNAL_ERROR);
  }
}




/***************** GeomTSTTCommon class methods *********************/

GeomTSTTCommon::GeomTSTTCommon( TSTTG::Geometry& geom ) throw ( TSTTB::Error )
  : geomIface( geom ),
    geomHandles( alloc_sidl_vector<void*>(1) ),
    positionsIn( alloc_sidl_vector<double>(3) )
{
}

GeomTSTTCommon::~GeomTSTTCommon() {}



void GeomTSTTCommon::move_to( void* geom, Vector3D& coord ) const
                            throw ( TSTTB::Error )
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
  memcpy( convert_from_sidl_vector(positionsIn), &coord, 3*sizeof(double) );
  positionsOut = convert_to_sidl_vector( reinterpret_cast<double*>(&coord), 3 );
  *convert_from_sidl_vector( geomHandles ) = geom;

  int junk = 3;
  geomIface.gentityClosestPoint( geomHandles,  1,
                                 positionsIn,  3,
                                 positionsOut, junk );
}

 
 
void GeomTSTTCommon::normal( void* geom, Vector3D& coord ) const
                            throw ( TSTTB::Error )
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
  memcpy( convert_from_sidl_vector(positionsIn), &coord, 3*sizeof(double) );
  normalsOut = convert_to_sidl_vector( reinterpret_cast<double*>(&coord), 3 );
  *convert_from_sidl_vector( geomHandles ) = geom;

  int junk;
  geomIface.gentityNormal( geomHandles, 1,
                           positionsIn, 3,
                           normalsOut,  junk );
}
 
void GeomTSTTCommon::normal( void* geom, Vector3D coords[], unsigned count ) const
                            throw ( TSTTB::Error )
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
  if (positionsIn.upper(0)+1 < 3*(int)count)
    positionsIn = positionsIn.create1d( 3*count );
  memcpy( convert_from_sidl_vector(positionsIn), coords, 3*count*sizeof(double) );
  
  normalsOut = convert_to_sidl_vector( reinterpret_cast<double*>(coords), 3*count );
  
  if (geomHandles.upper(0)+1 < (int)count)
    geomHandles = geomHandles.create1d( count );
  void** ptr = convert_from_sidl_vector( geomHandles );
  for (void** end = ptr + count; ptr != end; ++ptr)
    *ptr = geom;

  int junk;
  geomIface.gentityNormal( geomHandles, count,
                           positionsIn, 3*count,
                           normalsOut,  junk );
}

void GeomTSTTCommon::closest_and_normal( void* geom, 
                                         const Vector3D& position,
                                         Vector3D& closest,
                                         Vector3D& normal ) const
                                       throw (TSTTB::Error)
{
    // going to assume this in the following reinterpret_cast, so
    // check to make sure it is true
  assert( sizeof(Vector3D) == 3*sizeof(double) );
  
  memcpy( convert_from_sidl_vector(positionsIn), &position, 3*sizeof(double) );
  positionsOut = convert_to_sidl_vector( reinterpret_cast<double*>(&closest), 3 );
  normalsOut = convert_to_sidl_vector( reinterpret_cast<double*>(&normal), 3 );
  *convert_from_sidl_vector( geomHandles ) = geom;
  
  int junk1, junk2;
  geomIface.gentityClosestPointAndNormal( geomHandles,  1,
                                          positionsIn,  3,
                                          positionsOut, junk1,
                                          normalsOut,   junk2 );
}

} // namespace Mesquite

