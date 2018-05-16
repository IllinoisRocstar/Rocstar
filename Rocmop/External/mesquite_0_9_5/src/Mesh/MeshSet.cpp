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
// -*- Mode : c++; tab-width: 2; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 2 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 16-May-02 at 10:26:21
//  LAST-MOD: 22-Jun-04 at 14:31:34 by Thomas Leurent
//
/*! \file MeshSet.cpp

\brief This files implements all the memory management issues related
to the copy of the original TSTT (or other maybe) mesh entity handles
into Mesquite.
That copy is of course encapsulated in the MeshSet class.
  
    \author Thomas Leurent
    \date 2002-05-16  
 */
#include "Mesquite.hpp"
#include "MeshSet.hpp"
#include "QualityImprover.hpp"
#include "MsqError.hpp"
#include "MsqDebug.hpp"
#include "TopologyInfo.hpp"

#ifdef MSQ_USE_OLD_IO_HEADERS
#include <fstream.h>
#include <string.h>
#include <iomanip.h>
#else
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
#endif

namespace Mesquite {

MeshSet::MeshSet() :
  vertexIterator(NULL),
  spaceDim(0),
  csrOffsets(0),
//  csrData(0),
//  vertArray(NULL),
//  elemArray(NULL),
  elemTopologies(0),
  vertexOnBoundary(0),
  csrOffsetsSize(0),
//  csrDataSize(0),
  vertArraySize(0),
  elemArraySize(0),
  mDomain(NULL)
{
  cullFlag=MsqVertex::MSQ_SOFT_FIXED;
}

MeshSet::~MeshSet()
{
    // Delete the vertex iterator
  delete vertexIterator;
    // Release all of our meshes
  list<Mesquite::Mesh*>::iterator it = meshSet.begin();
  while(!(it == meshSet.end()))
    (*it++)->release();
  
    // Delete cache arrays
  delete [] csrOffsets;
//  delete [] csrData;
//  delete [] vertArray;
//  delete [] elemArray;
  delete [] elemTopologies;
  delete [] vertexOnBoundary;
}


/*! \fn  MeshSet::add_mesh(Mesquite::Mesh* mh, MsqError &err)

    adds a Mesquite::Mesh to the MeshSet. If used several times,
    it concatenates the vertices from several meshes into
    one Mesquite::MeshSet.
  */
void MeshSet::add_mesh(Mesquite::Mesh* mesh, MsqError &err)
{
    // sets MeshSet::SpaceDim
  int dim = mesh->get_geometric_dimension(err); MSQ_ERRRTN(err);
  if (spaceDim == 0) 
    spaceDim = dim;
  else if (dim != spaceDim)
  {
    MSQ_SETERR(err)( "Meshes of different dimensions added to the same MeshSet.",
                     MsqError::INVALID_STATE );
    return;
  }
  
    // adds the Mesh* to the MeshSet.
  meshSet.push_front(mesh);
}


/*! \fn MeshSet::reset(MsqError &err)

    Resets the MeshSet object.
    The current vertex is set back to the first vertex in the first mesh handle.
  */
void MeshSet::reset(MsqError& err)
{
    // If we have at least one mesh...
  if (meshSet.size())
  {
      // If we aren't already on the first mesh...
    if (!vertexIterator || !(currentMesh == meshSet.begin()))
    {
      currentMesh = meshSet.begin();
      delete vertexIterator;
      vertexIterator = (*currentMesh)->vertex_iterator(err); MSQ_ERRRTN(err);
    }
    else // We ARE on the first mesh...
    {
        // ...so we can re-use the iterator.
      vertexIterator->restart();
    }
  }
  else
  {
      // This is probably redundant...
    vertexIterator = NULL;
  }
}


void MeshSet::set_domain_constraint(MeshDomain* domain, MsqError &/*err*/)
{
    mDomain = domain;
}



/*! \fn MeshSet::get_next_patch(PatchData &pd, MsqError &err )
  \brief This function fills up a PatchData object with the mesh information
  necessary for optimization algorithms.
  The return value is true as long there exists a next patch, false otherwise.
  The list culling is performed in this function.
  
  This function is a friend of the PatchData class. Therefore the PatchData
  object passed as an argument is filled up directly.
  
  \param PatchData &pd: this is the PatchData object that will be filled up.
*/
bool MeshSet::get_next_patch(PatchData &pd,
                             PatchDataParameters &pd_params,
                             MsqError &err )
{
  MSQ_FUNCTION_TIMER( "MeshSet::get_next_patch" );

    // get rid of previous Patch information (but keep memory allocated).
  pd.clear();

    // Mark this MeshSet as the originator
  pd.meshSet = this;
  pd.domainSet = (mDomain != NULL);

  bool result = false;
  switch (pd_params.get_patch_type())
  {
    case PatchData::ELEMENTS_ON_VERTEX_PATCH:
      result = get_next_elem_on_vert_patch( pd, pd_params, err );
      break;
    case PatchData::GLOBAL_PATCH:
      result = get_next_global_patch( pd, pd_params, err );
      break;
    default:
      MSQ_SETERR(err)( "no implementation for specified patch type.",
                       MsqError::NOT_IMPLEMENTED );
      result = false;
      break;
  }

  return !MSQ_CHKERR(err) && result;
}


bool MeshSet::get_next_elem_on_vert_patch( PatchData& pd,
                                           PatchDataParameters &pd_params,
                                           MsqError& err )
{
  size_t i;

    // Get the patch parameters.
  long unsigned int culling_method_bits = pd_params.get_culling_method_bits();

    //variable to store the center vertex's fixed flag
  MsqVertex::FlagMask center_fixed_byte;
    // Make sure we're only getting a patch depth of 1
  int num_layers = pd_params.get_nb_layers(err); 
  if (MSQ_CHKERR(err)) return false;

  if ((unsigned)num_layers > (unsigned)1)
  {
    MSQ_SETERR(err)( "no implementation for patch depth > 1.",
                     MsqError::NOT_IMPLEMENTED ); 
    return false;
  }

    // Set the patch type
  pd.mType = PatchData::ELEMENTS_ON_VERTEX_PATCH;
  pd.domainHint = NO_DOMAIN_HINT;
  if (mDomain)
    pd.domainHint = mDomain->hint();

    // If this is our first time through the mesh,
    // initialize everything.
  if (!vertexIterator)
  {
    reset(err);
    if (MSQ_CHKERR(err)) return false;
  }

    // currentVertex is pointing at next potential center vertex.
    // Move forward in the list of vertices if necessary.
  bool next_vertex_identified = false;
  while (!next_vertex_identified)
  {
      // Move to next mesh if necessary
    while (vertexIterator->is_at_end())
    {
      delete vertexIterator;
      ++currentMesh;
      if (currentMesh == meshSet.end())
      {
        vertexIterator = NULL;
        return false;
      }
      vertexIterator = (*currentMesh)->vertex_iterator(err); MSQ_CHKERR(err);
    }

    bool on_bnd = false;
    bool is_mid = false;
    Mesquite::Mesh::VertexHandle vtx = **vertexIterator;
    (*currentMesh)->vertices_are_on_boundary(&vtx, &on_bnd, 1, err);
    if (MSQ_CHKERR(err)) return false;
    //(*currentMesh)->vertices_are_midnodes(&vtx, &is_mid, 1, err);
    //if (MSQ_CHKERR(err)) return false;

      //always skip midnodes
    if (is_mid)
    {
       vertexIterator->operator++();
    }
      //if this is a 'boundary' fixed flag, skip it now
    else if ((culling_method_bits & PatchData::NO_BOUNDARY_VTX)
             && (on_bnd==true))
    {
        ++(*vertexIterator);
    }
    else if ((culling_method_bits & PatchData::NO_INTERIOR_VTX)
             && (on_bnd==false))
    {
        ++(*vertexIterator);
    }
      //otherwise we check to see if this vertex has been culled
    else{
        //get the fixed_bit_flag for the center vertex
      (*currentMesh)->vertex_get_byte(**vertexIterator,&center_fixed_byte, err);
      if (MSQ_CHKERR(err)) return false;

        //remove the hard fixed flag if it has been set
      center_fixed_byte &= ~(MsqVertex::MSQ_HARD_FIXED);
        //if it is culled, skip it
      if(center_fixed_byte & cullFlag)
      {
        ++(*vertexIterator);
      }
      else
      {
          // We found the right one
        next_vertex_identified = true;
      }//end else (vertex was not fixed [boundary] or culled)
    }//end else (iterator was not at the end and vertex was not boundary)  
  }//end while (!next_vertex_identified)
  Mesh::VertexHandle vertex = **vertexIterator;
  vertexIterator->operator++();
  if(num_layers == 0 ){
    pd.vertexArray.resize( 1 );
    MsqVertex* pd_vert_array = pd.get_vertex_array(err);
    (*currentMesh)->vertices_get_coordinates(&vertex,
                                             pd_vert_array,
                                             1,
                                             err); MSQ_ERRZERO(err);
    pd_vert_array[0].vertexBitFlags=center_fixed_byte;
    
    pd.vertexHandlesArray.resize(1);
    pd.vertexHandlesArray[0]=vertex;
    
    pd.initialize_data( NULL, err ); MSQ_ERRZERO(err);
    return true;
  }
      // Get the number of elements in this vertex
  size_t num_elems =
    (*currentMesh)->vertex_get_attached_element_count(vertex, err);
  if (MSQ_CHKERR(err)) return false;
  pd.elementHandlesArray.resize( num_elems );
 
    // Get the elements attached to this vertex
  if (elemArraySize < num_elems)
  {
    delete [] elemTopologies;
    elemTopologies = new EntityTopology[num_elems];
    elemArraySize = num_elems;
  }

  (*currentMesh)->vertex_get_attached_elements(vertex,
                                               &pd.elementHandlesArray[0],
                                               num_elems, err);
  if (MSQ_CHKERR(err)) return false;

    // Get the topologies of those elements
  (*currentMesh)->elements_get_topologies(&pd.elementHandlesArray[0],
                                          elemTopologies,
                                          num_elems, err);
  if (MSQ_CHKERR(err)) return false;

    // Figure out how many vertices we need to allocate
  //size_t num_vert_uses = 1;
  //size_t i;
  //for (i = 0; i < num_elems; ++i)
  //  num_vert_uses += vertices_in_topology(elemTopologies[i]);
  size_t num_vert_uses = (*currentMesh)->
    get_vertex_use_count( &pd.elementHandlesArray[0], num_elems, err ); 
  MSQ_ERRZERO(err);

    // All elems share at least 1 vertex (the center vertex).  The
    // center vertex is used 1 time, but it was counted num_elems times.
  size_t num_verts = num_vert_uses - num_elems + 1;
  pd.vertexHandlesArray.resize( num_verts );
  pd.elementArray.resize( num_elems );
  pd.elemConnectivityArray.resize( num_vert_uses );

    // Get the vertices attached to those elements
  
  if (csrOffsetsSize < num_elems + 1)
  {
    delete [] csrOffsets;
    csrOffsets = new size_t[num_elems + 1];
    csrOffsetsSize = num_elems + 1;
  }
  (*currentMesh)->elements_get_attached_vertices(&pd.elementHandlesArray[0],
                                                 num_elems,
                                                 &pd.vertexHandlesArray[0],
                                                 num_verts,
                                                 &pd.elemConnectivityArray[0],
                                                 num_vert_uses,
                                                 csrOffsets,
                                                 err); 
  if (MSQ_CHKERR(err)) return false;
  
    // Update with actual vertex count
  pd.vertexHandlesArray.resize( num_verts );

    // Put the elements into the PatchData
  MsqMeshEntity* pd_elem_array = pd.get_element_array(err);
  for (i = 0; i < num_elems; ++i)
    pd_elem_array[i].set_element_type( elemTopologies[i] );
  pd.initialize_data( csrOffsets, err ); MSQ_ERRZERO(err);

    // Get the coordinates of the vertices and its flags.
  pd.vertexArray.resize( num_verts );
  MsqVertex* pd_vert_array = pd.get_vertex_array(err);
    //get the coordinates
  (*currentMesh)->vertices_get_coordinates(&pd.vertexHandlesArray[0],
                                           pd_vert_array,
                                           num_verts,
                                           err);
  if (MSQ_CHKERR(err)) return false;
  for (i = 0; i < num_verts; i++)
  {

      // If it's not the center vertex, mark it as hard fixed
    if (pd.vertexHandlesArray[i] != vertex)
    {
        // Get its flags
      (*currentMesh)->vertex_get_byte(pd.vertexHandlesArray[i],
                                      &(pd_vert_array[i].vertexBitFlags),
                                      err);
      if (MSQ_CHKERR(err)) return false;
      pd_vert_array[i].vertexBitFlags |= MsqVertex::MSQ_HARD_FIXED;
    }
      //else it is the center vertex.  We therefore already have
      //the fixed flag stored center_fixed_byte.  The hard fixed
      //flag has already been removed (when flag was retreived).
    else{
      pd_vert_array[i].vertexBitFlags = (center_fixed_byte);
    }
  }

  return true;
}


bool MeshSet::get_next_global_patch( PatchData& pd,
                                     PatchDataParameters& pd_params,
                                     MsqError& err )
{
    // We only support global patches for a single Mesh
  if (meshSet.size() != 1)
  {
    MSQ_SETERR(err)( 
       "Global patches only supported for single-Mesh MeshSets.",
       MsqError::NOT_IMPLEMENTED );
    return false;
  }

  pd.mType = PatchData::GLOBAL_PATCH;
  pd.domainHint = NO_DOMAIN_HINT;
  if (mDomain)
    pd.domainHint = mDomain->hint();

    // for a global patch, we always reset to start of the mesh.
  reset(err);
  if (MSQ_CHKERR(err)) return false;

  size_t i;

    // Get sizes for mesh data
  size_t num_verts, num_elems, num_uses;
  (*currentMesh)->get_all_sizes( num_verts, num_elems, num_uses, err ); MSQ_ERRZERO(err);

    // Get handles and connectivity
  pd.vertexHandlesArray.resize( num_verts );
  pd.elementHandlesArray.resize( num_elems );
  pd.elemConnectivityArray.resize( num_uses );
  msq_std::vector<size_t> offsets(num_elems+1);
  (*currentMesh)->get_all_mesh( &pd.vertexHandlesArray[0], num_verts,
                                &pd.elementHandlesArray[0], num_elems,
                                &offsets[0], offsets.size(),
                                &pd.elemConnectivityArray[0], 
                                pd.elemConnectivityArray.size(),
                                err ); MSQ_ERRZERO(err);
  
    // Get element topologies
  pd.elementArray.resize( num_elems );
  msq_std::vector<EntityTopology> elem_topologies(num_elems);
  (*currentMesh)->elements_get_topologies( &pd.elementHandlesArray[0], 
                                           &elem_topologies[0], 
                                           num_elems, err );MSQ_ERRZERO(err);

    // Put them into the patch
  MsqMeshEntity* pd_elem_array = pd.get_element_array(err);MSQ_ERRZERO(err);
  for (i = 0; i < num_elems; ++i)
    pd_elem_array[i].set_element_type( elem_topologies[i] );
    
    // Complete connectivity data in patch
  pd.initialize_data( &offsets[0], err ); MSQ_ERRZERO(err);


    // Get vertex coordinates
  pd.vertexArray.resize( num_verts );
  MsqVertex* pd_vert_array = pd.get_vertex_array(err);MSQ_ERRZERO(err);
  (*currentMesh)->vertices_get_coordinates(&pd.vertexHandlesArray[0],
                                           pd_vert_array,
                                           num_verts,
                                           err); MSQ_ERRZERO(err);

    // Get vertex boundary flag
  if (vertArraySize < num_verts)
  {
    delete [] vertexOnBoundary;
    vertArraySize = num_verts;
    vertexOnBoundary = new bool[vertArraySize];
  }
  (*currentMesh)->vertices_are_on_boundary( &pd.vertexHandlesArray[0],
                                            vertexOnBoundary,
                                            num_verts,
                                            err );MSQ_ERRZERO(err);

  for (i = 0; i < num_verts; i++)
  {
      // Get its flags
    /*(*currentMesh)->vertex_get_byte(vertArray[i],
                                    &(pd_vert_array[i].vertexBitFlags),
                                    err); MSQ_CHKERR(err);*/
      // Set its hard-fixed flag
    if (/*(*currentMesh)->vertex_is_fixed(vertArray[i], err) ||*/
       vertexOnBoundary[i])
    {
      pd_vert_array[i].vertexBitFlags |= MsqVertex::MSQ_HARD_FIXED;
    }
    else
    {
      pd_vert_array[i].vertexBitFlags &= ~(MsqVertex::MSQ_HARD_FIXED);
    }
  }

  return true;
}

// Currently, the only thing supported is updating each vertices
// coordinates and flags.  Connectivity changes aren't supported yet.
void Mesquite::MeshSet::update_mesh(const PatchData &pd, MsqError &err)
{
  MSQ_FUNCTION_TIMER( "MeshSet::update_mesh" );
  if (pd.num_nodes() == 0)
    return;
  
  size_t i;
  
  switch (pd.type())
  {
    // If the patch type is marked as local,
    // all handles belong to the currentMesh.
  case PatchData::ELEMENTS_ON_VERTEX_PATCH:
      // For each vertex, update the coordinates
        // and the "mesquite byte".
      for (i = 0; i < pd.num_nodes(); i++)
      {
        if(!pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_HARD_FIXED))
        {
          (*currentMesh)->vertex_set_coordinates(pd.vertexHandlesArray[i],
                                                 pd.vertexArray[i],
                                                 err); MSQ_ERRRTN(err);
        }
        
        (*currentMesh)->vertex_set_byte(pd.vertexHandlesArray[i],
                                        pd.vertexArray[i].vertexBitFlags,
                                        err); MSQ_ERRRTN(err);
      }
      break;
      
    // If the patch type is marked as global,
    // the handles may belong to more than
    // one Mesh.
  case PatchData::GLOBAL_PATCH:
    {
      list<Mesquite::Mesh*>::iterator mesh_itr = meshSet.begin();
      assert( mesh_itr != meshSet.end() );
      Mesquite::Mesh* cur_mesh = *mesh_itr;
      Mesquite::VertexIterator *vert_itr = cur_mesh->vertex_iterator(err);
      MSQ_ERRRTN(err);
      for (i = 0; i < pd.num_nodes(); i++)
      {
        if (vert_itr->is_at_end())
        {
          mesh_itr++;
          if ( mesh_itr==meshSet.end() )
            return;
          cur_mesh = *mesh_itr;
          delete vert_itr;
          vert_itr = cur_mesh->vertex_iterator(err); MSQ_ERRRTN(err);
        }
        if(!pd.vertexArray[i].is_flag_set( MsqVertex::MSQ_HARD_FIXED))
        {
          cur_mesh->vertex_set_coordinates(pd.vertexHandlesArray[i],
                                           pd.vertexArray[i],
                                           err); MSQ_ERRRTN(err);
        }
        cur_mesh->vertex_set_byte(pd.vertexHandlesArray[i],
                                  pd.vertexArray[i].vertexBitFlags,
                                  err); MSQ_ERRRTN(err);
      }
      delete vert_itr;
    }
    break;
  default:
    {
      MSQ_SETERR(err)("PatchData Type not accepted yet.", MsqError::NOT_IMPLEMENTED);
      break;
    }
  }
}

bool MeshSet::clear_all_soft_fixed_flags(MsqError &err)
{
    //variable to store the center vertex's fixed flag
  MsqVertex::FlagMask fixed_byte;
  bool finished_with_vertices=false;
    // initialize everything.
  if (!vertexIterator)
  {  reset(err); if(MSQ_CHKERR(err)) return false; }
    // currentVertex is pointing at next potential center vertex.
    
  while(!finished_with_vertices){
      // Move to next mesh if necessary
    if (vertexIterator->is_at_end())
    {
      delete vertexIterator;
      ++currentMesh;
      if (currentMesh == meshSet.end())
      {
        vertexIterator = NULL;
        finished_with_vertices=true;
      }
      if(!finished_with_vertices){
        vertexIterator = (*currentMesh)->vertex_iterator(err); 
        if (MSQ_CHKERR(err)) return false;
      }
    }
      //otherwise we check to see if this vertex has been culled
    else{
        //get the fixed_bit_flag 
      (*currentMesh)->vertex_get_byte(**vertexIterator,&fixed_byte, err);
      if (MSQ_CHKERR(err)) return false;
      fixed_byte &= (~MsqVertex::MSQ_SOFT_FIXED);
      (*currentMesh)->vertex_set_byte(**vertexIterator,fixed_byte, err);
      if (MSQ_CHKERR(err)) return false;
      ++(*vertexIterator);
    }
  }
  return true;
}



/* ************************************************************************* */
/* ************* Mesh Files can be written directly from the MeshSet ******* */
/* *************      Various formats are available below            ******* */
/* ************************************************************************* */



/*! Writes a VTK file directly from the MeshSet.
    This means that any mesh imported successfully into Mesquite
    can be outputed in VTK format.
    This is not geared for performance, since it has to load a global Patch from
    the mesh to write a mesh file. 
*/
void MeshSet::write_vtk(const char* out_filename,
                   Mesquite::MsqError &err)
{
    // Open the file
  msq_stdio::ofstream file(out_filename);
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // loads a global patch
  PatchData pd;
  PatchDataParameters pd_params;
  pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err); MSQ_ERRRTN(err);
  pd_params.no_culling_method();
  get_next_patch(pd, pd_params, err); MSQ_ERRRTN(err);
    
    // Write a header
  file << "# vtk DataFile Version 2.0\n";
  file << "Mesquite Mesh " << out_filename << " .\n";
  file << "ASCII\n";
  file << "DATASET UNSTRUCTURED_GRID\n";
  
    // Write vertex coordinates
  file << "POINTS " << pd.num_nodes() << " float\n";
  size_t i;
  for (i = 0; i < pd.num_nodes(); i++)
  {
    file << pd.vertexArray[i][0] << ' '
         << pd.vertexArray[i][1] << ' '
         << pd.vertexArray[i][2] << '\n';
  }
  
    // Write out the connectivity table
  size_t connectivity_size = 0;
  for (i = 0; i < pd.num_elements(); ++i)
    connectivity_size += pd.elementArray[i].node_count()+1;
    
  file << "CELLS " << pd.num_elements() << ' ' << connectivity_size << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    std::vector<size_t> vtx_indices;
    pd.elementArray[i].get_node_indices(vtx_indices);
    file << vtx_indices.size();
    for (msq_stdc::size_t j = 0; j < vtx_indices.size(); ++j)
    {
      file << ' ' << vtx_indices[j];
    }
    file << '\n';
  }
  
    // Write out the element types
  file << "CELL_TYPES " << pd.num_elements() << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    unsigned char type_id = 0;
    switch (pd.elementArray[i].get_element_type())
    {
      case Mesquite::TRIANGLE:
        type_id = 5;
        break;
      case Mesquite::QUADRILATERAL:
        type_id = 9;
        break;
      case Mesquite::TETRAHEDRON:
        type_id = 10;
        break;
      case Mesquite::HEXAHEDRON:
        type_id = 12;
        break;
    default:
      MSQ_SETERR(err)("element type not implemented",MsqError::NOT_IMPLEMENTED);
      break;
    }
    file << (int)type_id << '\n';
  }
  
    // Write out which points are fixed.
  file << "POINT_DATA " << pd.num_nodes()
       << "\nSCALARS fixed float\nLOOKUP_TABLE default\n";
  for (i = 0; i < pd.num_nodes(); ++i)
  {
    if (pd.vertexArray[i].is_free_vertex())
      file << "0\n";
    else
      file << "1\n";
  }
  
    // Close the file
  file.close();
}



/*! Writes a gnuplot file directly from the MeshSet.
    This means that any mesh imported successfully into Mesquite
    can be outputed in gnuplot format.

    Within gnuplot, use \b plot 'file1.gpt' w l, 'file2.gpt' w l  
    
    This is not geared for performance, since it has to load a global Patch from
    the mesh to write a mesh file. 
*/
void MeshSet::write_gnuplot(const char* out_filebase,
                   Mesquite::MsqError &err)
{
    // Open the file
  string out_filename = out_filebase;
  out_filename += ".gpt";
  ofstream file(out_filename.c_str());
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // loads a global patch
  PatchData pd;
  PatchDataParameters pd_params;
  pd_params.set_patch_type(PatchData::GLOBAL_PATCH, err); MSQ_ERRRTN(err);
  pd_params.no_culling_method();
  get_next_patch(pd, pd_params, err); MSQ_ERRRTN(err);
    
    // Write a header
  file << "\n";
  
  for (size_t i=0; i<pd.num_elements(); ++i)
  {
    std::vector<size_t> vtx_indices;
    pd.elementArray[i].get_node_indices(vtx_indices);
    for (size_t j = 0; j < vtx_indices.size(); ++j)
    {
      file << pd.vertexArray[vtx_indices[j]][0] << ' '
           << pd.vertexArray[vtx_indices[j]][1] << ' '
           << pd.vertexArray[vtx_indices[j]][2] << '\n';
    }
      file << pd.vertexArray[vtx_indices[0]][0] << ' '
           << pd.vertexArray[vtx_indices[0]][1] << ' '
           << pd.vertexArray[vtx_indices[0]][2] << '\n';
    file << '\n';
  }
  
    // Close the file
  file.close();
}

} // namespace Mesquite
