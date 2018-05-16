;/* ***************************************************************** 
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
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 18-Jun-04 at 11:36:07 by Thomas Leurent
//
/*! \file MsqMeshEntity.cpp

\brief All elements in Mesquite are of type MsqMeshEntity. Their associated
functionality is implemented in this file. 
  
    \author Thomas Leurent
    \author Michael Brewer
    \author Darryl Melander
    \date 2002-05-16  
 */

#include "Mesquite.hpp"
#include "MsqMeshEntity.hpp"
#include "MsqVertex.hpp"
#include "PatchData.hpp"
#include "MeshSet.hpp"

#ifdef MSQ_USE_OLD_STD_HEADERS
#  include <vector.h>
#else
#  include <vector>
   using std::vector;
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  include <ostream.h>
#else
#  include <ostream>
   using std::ostream;
   using std::endl;
#endif

namespace Mesquite {

//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
void Mesquite::MsqMeshEntity::get_vertex_indices(vector<size_t> &vertices) const
{
  vertices.resize( vertex_count() );
  memcpy( &vertices[0], vertexIndices, vertex_count() * sizeof(size_t) );
}

//! Gets the indices of the vertices of this element.
//! The indices are only valid in the PatchData from which
//! this element was retrieved.
//! The order of the vertices is the canonical order for this
//! element's type.
//! The indices are placed appended to the end of the list.
//! The list is not cleared before appending this entity's vertices.
void Mesquite::MsqMeshEntity::append_vertex_indices(vector<size_t> &vertex_list) const
{
  vertex_list.insert(vertex_list.end(),
                     vertexIndices,
                     vertexIndices + vertex_count());
}

void MsqMeshEntity::get_node_indices( vector<size_t>& indices ) const
{
  indices.resize( node_count() );
  memcpy( &indices[0], vertexIndices, node_count() * sizeof( size_t ) );
}

void MsqMeshEntity::append_node_indices( vector<size_t>& indices ) const
{
  indices.insert( indices.end(), vertexIndices, vertexIndices + node_count() );
}


/*! The centroid of an element containing n vertices with equal masses is located at
  \f[ \b{x} = \frac{ \sum_{i=1}^{n} \b{x}_i }{ n }  \f]
  where \f$ \b{x}_i  ,\, i=1,...,n\f$ are the vertices coordinates.
*/
void MsqMeshEntity::get_centroid(Vector3D &centroid, const PatchData &pd, MsqError &err) const
{
  const MsqVertex* vtces = pd.get_vertex_array(err); MSQ_ERRRTN(err);
  size_t nve = vertex_count();
  for (size_t i=0; i<nve; ++i)
    centroid += vtces[vertexIndices[i]];
  centroid /= nve;
}
  

/*!fills array of Vector3D's with the jacobian vectors and the 
  number of jacobian vecotors.*/
void Mesquite::MsqMeshEntity::compute_weighted_jacobian(PatchData &pd,
                                                        Vector3D &sample_point,
                                                        Vector3D jacobian_vectors[],
                                                        short &num_jacobian_vectors,
                                                        MsqError &err )
{
    // v_v is just an alias for vertexIndices
  const size_t* v_v = &vertexIndices[0];
  
    //   vector<size_t> v_v;
    //   get_vertex_indices(v_v);
  MsqVertex *vertices=pd.get_vertex_array(err); MSQ_ERRRTN(err);
  switch (mType)
  {
      //Note:: For the linear tri case we do not use sample pt.
    case TRIANGLE:
      jacobian_vectors[0].set(vertices[v_v[1]]-vertices[v_v[0]]);
      jacobian_vectors[1].set((2.0*vertices[v_v[2]]-vertices[v_v[0]]-
                               vertices[v_v[1]])*MSQ_SQRT_THREE_INV);
      num_jacobian_vectors=2;
      break;
      
    case QUADRILATERAL:
      jacobian_vectors[0]=(vertices[v_v[1]]-vertices[v_v[0]]+sample_point[1]*
                           (vertices[v_v[2]]+vertices[v_v[0]]-vertices[v_v[3]]-
                            vertices[v_v[1]]));
      jacobian_vectors[1]=(vertices[v_v[3]]-vertices[v_v[0]]+sample_point[0]*
                           (vertices[v_v[2]]+vertices[v_v[0]]-vertices[v_v[3]]-
                            vertices[v_v[1]]));
      num_jacobian_vectors=2;
      break;
      
    case TETRAHEDRON:
      jacobian_vectors[0]=vertices[v_v[1]]-vertices[v_v[0]];
      jacobian_vectors[1]=(2.0*vertices[v_v[2]]-vertices[v_v[0]]-
                           vertices[v_v[1]])*MSQ_SQRT_THREE_INV;
      jacobian_vectors[2]=(3.0*vertices[v_v[3]]-vertices[v_v[2]]-
                           vertices[v_v[1]]-vertices[v_v[0]])*
        MSQ_SQRT_TWO_INV*MSQ_SQRT_THREE_INV;
      num_jacobian_vectors=3;
      break;
      
    case HEXAHEDRON:
      
      jacobian_vectors[0]=vertices[v_v[1]]-vertices[v_v[0]]+
                           (sample_point[1]*(vertices[v_v[2]]+vertices[v_v[0]]-
                                             vertices[v_v[3]]-vertices[v_v[1]]))+
                           (sample_point[2]*(vertices[v_v[5]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[1]]))+
                           (sample_point[1]*sample_point[2]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]])),
      jacobian_vectors[1]=vertices[v_v[3]]-vertices[v_v[0]]+
                           (sample_point[0]*(vertices[v_v[2]]+vertices[v_v[0]]-
                                             vertices[v_v[3]]-vertices[v_v[1]]))+
                           (sample_point[2]*(vertices[v_v[7]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[3]]))+
                           (sample_point[0]*sample_point[2]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]]));
                           
      jacobian_vectors[2]=vertices[v_v[4]]-vertices[v_v[0]]+
                           (sample_point[0]*(vertices[v_v[5]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[1]]))+
                           (sample_point[1]*(vertices[v_v[7]]+vertices[v_v[0]]-
                                             vertices[v_v[4]]-vertices[v_v[3]]))+
                           (sample_point[0]*sample_point[1]*(vertices[v_v[6]]+
                                                             vertices[v_v[4]]+
                                                             vertices[v_v[3]]+
                                                             vertices[v_v[1]]-
                                                             vertices[v_v[7]]-
                                                             vertices[v_v[5]]-
                                                             vertices[v_v[2]]-
                                                             vertices[v_v[0]]));

      num_jacobian_vectors=3;
      break;
      
    default:
      MSQ_SETERR(err)(
        "Compute_weighted_jacobian not yet defined for this entity.",
        MsqError::NOT_IMPLEMENTED);
  } 
  
}

//! Appends the coordinates of the sample point to 'coords'.
/*! Places Vector3Ds holding the sample point for a given EvaluationMode
  and element type combination into a given vector of Vector3D.
  \param QualityMetric::EvaluationMode mode Specifies the type of sample
  points being used.
  \param vector<Vecotr3D> &coords A vector of Vector3D passed by
  reference which is used to store the sample points.
*/
void Mesquite::MsqMeshEntity::get_sample_points(QualityMetric::ElementEvaluationMode mode,
                                      vector<Vector3D> &coords,
                                      MsqError &err){
  switch (mType)
  {
    case TRIANGLE:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(3);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO, 0.0));
          break;
          
            //The following need to be verified
        case (QualityMetric::LINEAR_GAUSS_POINTS):
          coords.reserve(1);
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          break;
          
        case (QualityMetric::QUADRATIC_GAUSS_POINTS):
          coords.reserve(3);
          coords.push_back(Vector3D(0.5, 0.0, 0.0));	
          coords.push_back(Vector3D(0.75, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          coords.push_back(Vector3D(0.25, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          break;
          
        case (QualityMetric::CUBIC_GAUSS_POINTS):
          coords.reserve(4);
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/2.0, 0.0));
          coords.push_back(Vector3D(0.2, 2.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          coords.push_back(Vector3D(0.8, 2.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          coords.push_back(Vector3D(0.5, 11.0* MSQ_SQRT_THREE_DIV_TWO/15.0,
                                    0.0));
          break;
        default:
            //return error saying sample points for mode not implem.    
          MSQ_SETERR(err)("Requested Sample Point Mode not implemented",
                          MsqError::NOT_IMPLEMENTED);
          return;
      }
      break;
    
    case QUADRILATERAL:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(4);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(1.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 1.0, 0.0));
          break;
            //THESE NEED TO BE VERIFIED
        case (QualityMetric::LINEAR_GAUSS_POINTS):
          coords.push_back(Vector3D(0.5, 0.5, 0.0));
          break;
        case (QualityMetric::QUADRATIC_GAUSS_POINTS):
        case (QualityMetric::CUBIC_GAUSS_POINTS):
        default:
            //return error saying sample points for mode not implem.
          MSQ_SETERR(err)("Requested Sample Point Mode not implemented",
                          MsqError::NOT_IMPLEMENTED);
          return;
      }
      break;

    case TETRAHEDRON:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(4);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO, 0.0));
          coords.push_back(Vector3D(0.5, MSQ_SQRT_THREE_DIV_TWO/3.0,
                                    MSQ_SQRT_TWO_DIV_SQRT_THREE));
          break;
        case (QualityMetric::LINEAR_GAUSS_POINTS):

        case (QualityMetric::QUADRATIC_GAUSS_POINTS):

        case (QualityMetric::CUBIC_GAUSS_POINTS):

        default:
            //return error saying sample points for mode not implem.    
          MSQ_SETERR(err)("Requested Sample Point Mode not implemented",
                          MsqError::NOT_IMPLEMENTED);
          return;
      }   
      break;

      case HEXAHEDRON:
      switch (mode)
      {
        case (QualityMetric::ELEMENT_VERTICES):
          coords.reserve(8);
          coords.push_back(Vector3D(0.0, 0.0, 0.0));	
          coords.push_back(Vector3D(1.0, 0.0, 0.0));
          coords.push_back(Vector3D(1.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 1.0, 0.0));
          coords.push_back(Vector3D(0.0, 0.0, 1.0));	
          coords.push_back(Vector3D(1.0, 0.0, 1.0));
          coords.push_back(Vector3D(1.0, 1.0, 1.0));
          coords.push_back(Vector3D(0.0, 1.0, 1.0));
          break;
        case (QualityMetric::LINEAR_GAUSS_POINTS):

        case (QualityMetric::QUADRATIC_GAUSS_POINTS):

        case (QualityMetric::CUBIC_GAUSS_POINTS):

        default:
            //return error saying sample points for mode not implem.    
          MSQ_SETERR(err)("Requested Sample Point Mode not implemented",
                          MsqError::NOT_IMPLEMENTED);
          return;
      }   
      break;
      
    default:
        //return error saying sample points for mode not implem.
          MSQ_SETERR(err)("Requested Sample Point Mode not implemented",
                          MsqError::NOT_IMPLEMENTED);
          return;
  }
}

/*!
  \brief Computes the area of the given element.  Returned value is
  always non-negative.  If the entity passed is not a two-dimensional
  element, an error is set.*/
double MsqMeshEntity::compute_unsigned_area(PatchData &pd, MsqError &err) {
  MsqVertex* verts=pd.get_vertex_array(err);MSQ_ERRZERO(err);
  double tem=0.0;
  switch (mType)
  {
   
    case TRIANGLE:
      tem =  ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
              (verts[vertexIndices[2]]-verts[vertexIndices[0]])).length()/2.0;
      return tem;
      
    case QUADRILATERAL:
      tem = ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
             (verts[vertexIndices[3]]-verts[vertexIndices[0]])).length();
      tem += ((verts[vertexIndices[3]]-verts[vertexIndices[2]])*
              (verts[vertexIndices[1]]-verts[vertexIndices[2]])).length();
      return (tem/2.0);
      
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute unsigned area.",
                      MsqError::INVALID_ARG);
      return 0;
  }
  return 0;
}
                                            
/*!
  \brief Computes the volume of the given element.  Returned value is
  always non-negative.  If the entity passed is not a three-dimensional
  element, an error is set.*/
double MsqMeshEntity::compute_unsigned_volume(PatchData &pd, MsqError &err) {
  Vector3D sample_point(.5,.5,.5);
  Vector3D jac_vecs[3];
  short num_jacobian_vectors=-1;
  double tem=0;
  MsqVertex *verts = pd.get_vertex_array(err);MSQ_ERRZERO(err);
  switch (mType)
  {
    case TETRAHEDRON:
      tem = (verts[vertexIndices[3]]-verts[vertexIndices[0]])%
        ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
         (verts[vertexIndices[1]]-verts[vertexIndices[0]]))/6.0;
      return fabs(tem);
      
    case HEXAHEDRON:
      compute_weighted_jacobian(pd,sample_point,jac_vecs,
                                num_jacobian_vectors, err ); MSQ_ERRZERO(err);
      return fabs(jac_vecs[2]%(jac_vecs[0]*jac_vecs[1]));
      
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute unsigned volume.",
                       MsqError::INVALID_ARG);
     return 0;
  }
  return 0;
}


/*!
  \brief Computes the area of the given element.  Returned value can be
  negative.  If the entity passed is not a two-dimensional element, an
  error is set.*/
double MsqMeshEntity::compute_signed_area(PatchData &pd, MsqError &err) {
  MsqVertex* verts=pd.get_vertex_array(err);MSQ_ERRZERO(err);
  double tem=0.0;
  double tem2=0.0;
  Vector3D surface_normal;
  Vector3D cross_vec;
  size_t element_index=pd.get_element_index(this);
  
  switch (mType)
  {
    
    case TRIANGLE:
      cross_vec=((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
      (verts[vertexIndices[2]]-verts[vertexIndices[0]]));
      pd.get_domain_normal_at_element(element_index,surface_normal,err);
      MSQ_ERRZERO(err);
      tem =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem *= -1;
      }
      
      return tem;
      
    case QUADRILATERAL:
      cross_vec=((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
                 (verts[vertexIndices[3]]-verts[vertexIndices[0]]));
      pd.get_domain_normal_at_element(element_index,surface_normal,err);
      MSQ_ERRZERO(err);
      tem =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem *= -1;
      }
      cross_vec=((verts[vertexIndices[3]]-verts[vertexIndices[2]])*
                 (verts[vertexIndices[1]]-verts[vertexIndices[2]]));
      tem2 =  (cross_vec.length()/2.0);
        //if normals do not point in same general direction, negate area
      if(cross_vec%surface_normal<0){ 
        tem2 *= -1;
          //test to make sure surface normal existed
          //if(surface_normal.length_squared()<.5){
          //err.set_msg("compute_signed_area called without surface_normal available.");
          //}  
      }
      return (tem + tem2);
      
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute unsigned area.",
                       MsqError::INVALID_ARG);
      return 0;
  };
  return 0.0;
}
    
/*!
  \brief Computes the volume of the given element.  Returned value can be
  negative.  If the entity passed is not a three-dimensional element,
  an error is set.*/
double MsqMeshEntity::compute_signed_volume(PatchData &pd, MsqError &err) {
  Vector3D sample_point(.5,.5,.5);
  Vector3D jac_vecs[3];
  short num_jacobian_vectors=-1;
  double tem=0;
  MsqVertex *verts = pd.get_vertex_array(err);MSQ_ERRZERO(err);
  switch (mType)
  {
    case TETRAHEDRON:
      tem = (verts[vertexIndices[3]]-verts[vertexIndices[0]])%
        ((verts[vertexIndices[1]]-verts[vertexIndices[0]])*
         (verts[vertexIndices[2]]-verts[vertexIndices[0]]))/6.0;
      return tem;
      
    case HEXAHEDRON:
      compute_weighted_jacobian(pd,sample_point,jac_vecs,
                                num_jacobian_vectors, err );
      MSQ_ERRZERO(err);
      return (jac_vecs[2]%(jac_vecs[0]*jac_vecs[1]));
      
    default:
      MSQ_SETERR(err)("Invalid type of element passed to compute signed volume.",
                       MsqError::INVALID_ARG);
      return 0;
  };
  return 0.0;      
}


/*!Appends the indices (in the vertex array) of the vertices to connected
  to vertex_array[vertex_index] to the end of the vector vert_indices.
  The connected vertices are right-hand ordered as defined by the
  entity.
  
*/
void Mesquite::MsqMeshEntity::get_connected_vertices(size_t vertex_index,
                                                     vector<size_t> &vert_indices,
                                                     MsqError &err)
{
    //i iterates through elem's vertices
  int i=0;
    //index is set to the index in the vertexIndices corresponding
    //to vertex_index
  int index=-1;
  
  switch (mType)
  {
    case TRIANGLE:
      while(i<3)
      {
        if(vertexIndices[i]==vertex_index)
        {
          index=i;
          break;
        }
        ++i;
      }
      if(index>=0)
      {
        vert_indices.push_back(vertexIndices[(index+1)%3]);
        vert_indices.push_back(vertexIndices[(index+2)%3]);
      }
      
      break;
      
    case QUADRILATERAL:
      while(i<4)
      {
        if(vertexIndices[i]==vertex_index)
        {
          index=i;
          break;
        }
        ++i;
      }
      if(index>=0)
      {
        vert_indices.push_back(vertexIndices[(index+1)%4]);
        vert_indices.push_back(vertexIndices[(index+3)%4]);
      }
          
      break;
      
    case TETRAHEDRON:
      while(i<4)
      {
        if(vertexIndices[i]==vertex_index)
        {
          index=i;
          break;
        }
        ++i;
      }
      if(index>=0)
      {
        vert_indices.push_back(vertexIndices[(index+1)%4]);
        vert_indices.push_back(vertexIndices[(index+2)%4]);
        vert_indices.push_back(vertexIndices[(index+3)%4]);
      }
      
      break;
      
    case HEXAHEDRON:
      while(i<8)
      {
        if(vertexIndices[i]==vertex_index)
        {
          index=i;
          break;
        }
        ++i;
      }
      
      if(index>=0)
      {
        if (index<4)
        {
          vert_indices.push_back(vertexIndices[(index+1)%4]);
          vert_indices.push_back(vertexIndices[(index+3)%4]);
          vert_indices.push_back(vertexIndices[(index)+4]);
        }
        else
        {
          vert_indices.push_back(vertexIndices[(index+3)%4+4]);
          vert_indices.push_back(vertexIndices[(index+1)%4+4]);
          vert_indices.push_back(vertexIndices[(index)-4]);
        }
      }
      
      break;
      
  default:
    MSQ_SETERR(err)("Element type not available", MsqError::INVALID_ARG);
    break;
  }
}

/*! Gives the normal at the surface point corner_pt ... but if not available,
    gives the normalized cross product of corner_vec1 and corner_vec2. 
  */
/*
void MsqMeshEntity::compute_corner_normal(size_t corner,
                                          Vector3D &normal,
                                          PatchData &pd, 
                                          MsqError &err)
{
  if ( get_element_type()==TRIANGLE || get_element_type()==QUADRILATERAL ) 
  {
      // There are two cases where we cannot get a normal from the 
      // geometry that are not errors:
      // 1) There is no domain set
      // 2) The vertex is at a degenerate point on the geometry (e.g. 
      //     tip of a cone.)
    
    bool have_normal = false;

      // Get normal from domain
    if (pd.domain_set()) 
    {
      size_t index = pd.get_element_index(this);
      pd.get_domain_normal_at_corner( index, corner, normal, err );
      MSQ_ERRRTN(err);
      
      double length = normal.length();
      if (length > DBL_EPSILON) 
      {
        have_normal = true;
        normal /= length;
      }
    }

      // If failed to get normal from domain, calculate
      // from adjacent vertices.
    if ( !have_normal )
    {
      const int num_corner = this->vertex_count();
      const int prev_corner = (corner + num_corner - 1) % num_corner;
      const int next_corner = (corner + 1 ) % num_corner;
      const size_t this_idx = vertexIndices[corner];
      const size_t prev_idx = vertexIndices[prev_corner];
      const size_t next_idx = vertexIndices[next_corner];
      normal = (pd.vertex_by_index(next_idx) - pd.vertex_by_index(this_idx))
             * (pd.vertex_by_index(prev_idx) - pd.vertex_by_index(this_idx));
      normal.normalize();
    }
  }
  else
    MSQ_SETERR(err)("Should only be used for faces (tri, quads, ...).",
                    MsqError::INVALID_ARG);
}
*/

void MsqMeshEntity::compute_corner_normals( Vector3D normals[],
                                            PatchData &pd, 
                                            MsqError &err)
{
  EntityTopology type = get_element_type();
  if (type != TRIANGLE && type != QUADRILATERAL && type != POLYGON)
  {
      MSQ_SETERR(err)("Should only be used for faces (tri, quads, ...).",
                    MsqError::INVALID_ARG);
      return;
  }
  
  
    // There are two cases where we cannot get a normal from the 
    // geometry that are not errors:
    // 1) There is no domain set
    // 2) The vertex is at a degenerate point on the geometry (e.g. 
    //     tip of a cone.)

    // Get normal from domain
  if (pd.domain_set()) 
  {
    size_t index = pd.get_element_index(this);
    pd.get_domain_normals_at_corners( index, normals, err );
    MSQ_ERRRTN(err);
  }

    // Check if normals are valid (none are valid if !pd.domain_set())
  const unsigned count = vertex_count();
  for (unsigned i = 0; i < count; ++i)
  {
      // If got valid normal from domain, 
      // make it a unit vector and continue.
    if (pd.domain_set()) 
    {
      double length = normals[i].length();
      if (length > DBL_EPSILON)
      {
        normals[i] /= length;
        continue;
      }
    }
    
    const size_t prev_idx = vertexIndices[(i + count - 1) % count];
    const size_t this_idx = vertexIndices[i];
    const size_t next_idx = vertexIndices[(i + 1) % count];

      // Calculate normal using edges adjacent to corner
    normals[i] = (pd.vertex_by_index(next_idx) - pd.vertex_by_index(this_idx))
               * (pd.vertex_by_index(prev_idx) - pd.vertex_by_index(this_idx));
    normals[i].normalize();
  }
}

/*!  \param pd  The PatchData the element belongs to. It contains the vertices coords.
     \param c_m3d An array of Matrix3D objects. There should be one matrix per element corner
            (4 for a tet, 8 for an hex). Each column of the matrix will contain a vector
            corresponding to a corner edge. 
     \param num_m3d The number of matrices in the c_m3d array. The function will check this number
            corresponds to the number of corner in the element. If not, an error is set. 
  */
void MsqMeshEntity::compute_corner_matrices(PatchData &pd, Matrix3D A[], int num_m3d, MsqError &err )
{
  MsqVertex* vertices = pd.get_vertex_array(err); MSQ_ERRRTN(err); 
  const size_t* v_i = &vertexIndices[0];

  // If 2D element, we will get the surface normal 
  Vector3D normals[4], vec1, vec2, vec3, vec4;

  
  switch(get_element_type()){
    
  case TRIANGLE:
    if (num_m3d != 3) {
      MSQ_SETERR(err)("num_m3d incompatible with element type.", MsqError::INVALID_ARG); 
      return;
    }
    
    compute_corner_normals( normals, pd, err ); MSQ_ERRRTN(err);

    vec1 = vertices[v_i[1]]-vertices[v_i[0]];
    vec2 = vertices[v_i[2]]-vertices[v_i[0]];
    vec3 = vertices[v_i[2]]-vertices[v_i[1]];
    
    A[0].set_column(0, vec1);
    A[0].set_column(1, vec2);
    A[0].set_column(2, normals[0]*MSQ_3RT_2_OVER_6RT_3);
    
    A[1].set_column(0, vec3);
    A[1].set_column(1, -vec1);
    A[1].set_column(2, normals[1]*MSQ_3RT_2_OVER_6RT_3);

    A[2].set_column(0, -vec2);
    A[2].set_column(1, -vec3);
    A[2].set_column(2, normals[2]*MSQ_3RT_2_OVER_6RT_3);

    break;
    
  case QUADRILATERAL:
    if (num_m3d != 4)  {
      MSQ_SETERR(err)("num_m3d incompatible with element type.", MsqError::INVALID_ARG); 
      return;
    }
    vec1 = vertices[v_i[1]]-vertices[v_i[0]];
    vec2 = vertices[v_i[3]]-vertices[v_i[0]];
    vec3 = vertices[v_i[2]]-vertices[v_i[1]];
    vec4 = vertices[v_i[3]]-vertices[v_i[2]];
    
    compute_corner_normals( normals, pd, err ); MSQ_ERRRTN(err);

    A[0].set_column(0, vec1);
    A[0].set_column(1, vec2);
    A[0].set_column(2, normals[0]);
    
    A[1].set_column(0, vec3);
    A[1].set_column(1, -vec1);
    A[1].set_column(2, normals[1]);

    A[2].set_column(0, vec4);
    A[2].set_column(1, -vec3);
    A[2].set_column(2, normals[2]);

    A[3].set_column(0, -vec2);
    A[3].set_column(1, -vec4);
    A[3].set_column(2, normals[3]);

    break;
    
  case TETRAHEDRON:
    if (num_m3d != 4) {
      MSQ_SETERR(err)("num_m3d incompatible with element type.", MsqError::INVALID_ARG); 
      return;
    }
    A[0].set_column(0, vertices[v_i[1]]-vertices[v_i[0]]);
    A[0].set_column(1, vertices[v_i[2]]-vertices[v_i[0]]);
    A[0].set_column(2, vertices[v_i[3]]-vertices[v_i[0]]);
    
    A[1].set_column(0, vertices[v_i[0]]-vertices[v_i[1]]);
    A[1].set_column(1, vertices[v_i[3]]-vertices[v_i[1]]);
    A[1].set_column(2, vertices[v_i[2]]-vertices[v_i[1]]);

    A[2].set_column(0, vertices[v_i[3]]-vertices[v_i[2]]);
    A[2].set_column(1, vertices[v_i[0]]-vertices[v_i[2]]);
    A[2].set_column(2, vertices[v_i[1]]-vertices[v_i[2]]);

    A[3].set_column(0, vertices[v_i[2]]-vertices[v_i[3]]);
    A[3].set_column(1, vertices[v_i[1]]-vertices[v_i[3]]);
    A[3].set_column(2, vertices[v_i[0]]-vertices[v_i[3]]);

    break;
    
  /*
   case PYRAMID:
      //We compute the pyramid's "condition number" by averaging
      //the 4 tet's condition numbers, where the tets are created
      //by removing one of the four base vertices from the pyramid.
      //transform to origina v_i[0]
      temp_vec[3]=vertices[v_i[1]]-vertices[v_i[0]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[0]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[0]];
      //find AW_inverse
      temp_vec[0]=temp_vec[3];
      temp_vec[1]=temp_vec[4]-temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,pd,met_vals[0],err);
      if(!return_flag)
      return return_flag;
      //transform to origina v_i[1]
      temp_vec[3]=vertices[v_i[2]]-vertices[v_i[1]];
      temp_vec[4]=vertices[v_i[3]]-vertices[v_i[1]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[1]];
      //find AW_inverse
      temp_vec[0]=temp_vec[3]-temp_vec[4];
      temp_vec[1]=temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,pd,met_vals[1],err);
      if(!return_flag)
      return return_flag;
      //transform to origina v_i[1]     
      temp_vec[3]=vertices[v_i[3]]-vertices[v_i[2]];
      temp_vec[4]=vertices[v_i[0]]-vertices[v_i[2]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[2]];
      //find AW_inverse
      temp_vec[0]=-temp_vec[3];
      temp_vec[1]=temp_vec[3]-temp_vec[4];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,pd,met_vals[2],err);
      if(!return_flag)
      return return_flag;
      //transform to origina v_i[1]     
      temp_vec[3]=vertices[v_i[0]]-vertices[v_i[3]];
      temp_vec[4]=vertices[v_i[1]]-vertices[v_i[3]];
      temp_vec[5]=vertices[v_i[4]]-vertices[v_i[3]];
      //find AW_inverse
      temp_vec[0]=temp_vec[4]-temp_vec[3];
      temp_vec[1]=-temp_vec[3];
      temp_vec[2]=MSQ_SQRT_TWO*(temp_vec[5]-(temp_vec[4]/2.0));
      return_flag=condition_number_3d(temp_vec,pd,met_vals[3],err);
      fval=average_metrics(met_vals, 4, err);
      if(!return_flag)
      return return_flag;
      break;
    */
    
  case HEXAHEDRON:
    if (num_m3d != 8) {
      MSQ_SETERR(err)("num_m3d incompatible with element type.", MsqError::INVALID_ARG); 
      return;
    }
    A[0].set_column(0, vertices[v_i[1]]-vertices[v_i[0]]);
    A[0].set_column(1, vertices[v_i[3]]-vertices[v_i[0]]);
    A[0].set_column(2, vertices[v_i[4]]-vertices[v_i[0]]);

    A[1].set_column(0, vertices[v_i[2]]-vertices[v_i[1]]);
    A[1].set_column(1, vertices[v_i[0]]-vertices[v_i[1]]);
    A[1].set_column(2, vertices[v_i[5]]-vertices[v_i[1]]);

    A[2].set_column(0, vertices[v_i[3]]-vertices[v_i[2]]);
    A[2].set_column(1, vertices[v_i[1]]-vertices[v_i[2]]);
    A[2].set_column(2, vertices[v_i[6]]-vertices[v_i[2]]);
    
    A[3].set_column(0, vertices[v_i[0]]-vertices[v_i[3]]);
    A[3].set_column(1, vertices[v_i[2]]-vertices[v_i[3]]);
    A[3].set_column(2, vertices[v_i[7]]-vertices[v_i[3]]);

    A[4].set_column(0, vertices[v_i[7]]-vertices[v_i[4]]);
    A[4].set_column(1, vertices[v_i[5]]-vertices[v_i[4]]);
    A[4].set_column(2, vertices[v_i[0]]-vertices[v_i[4]]);
    
    A[5].set_column(0, vertices[v_i[4]]-vertices[v_i[5]]);
    A[5].set_column(1, vertices[v_i[6]]-vertices[v_i[5]]);
    A[5].set_column(2, vertices[v_i[1]]-vertices[v_i[5]]);
    
    A[6].set_column(0, vertices[v_i[5]]-vertices[v_i[6]]);
    A[6].set_column(1, vertices[v_i[7]]-vertices[v_i[6]]);
    A[6].set_column(2, vertices[v_i[2]]-vertices[v_i[6]]);
    
    A[7].set_column(0, vertices[v_i[6]]-vertices[v_i[7]]);
    A[7].set_column(1, vertices[v_i[4]]-vertices[v_i[7]]);
    A[7].set_column(2, vertices[v_i[3]]-vertices[v_i[7]]);

    break;
    
  default:
    MSQ_SETERR(err)("element type not implemented.", MsqError::NOT_IMPLEMENTED);
    return;
  }// end switch over element type
}

ostream& operator<<( ostream& stream, const MsqMeshEntity& entity )
{
  size_t num_vert = entity.vertex_count();
  stream << "MsqMeshEntity " << &entity << " with vertices ";
  for (size_t i = 0; i < num_vert; ++i)
    stream << entity.vertexIndices[i] << " ";
  stream << endl;
  return stream;
}

} // namespace Mesquite
