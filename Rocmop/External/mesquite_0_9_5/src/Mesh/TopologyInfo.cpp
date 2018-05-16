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

#include "MsqError.hpp"
#include "TopologyInfo.hpp"

#include <string.h>
#include <assert.h>

namespace Mesquite {

TopologyInfo TopologyInfo::instance;


TopologyInfo::TopologyInfo()
{
  memset( dimMap, 0, sizeof(dimMap) );
  memset( adjMap, 0, sizeof(adjMap) );
  memset( edgeMap, 0, sizeof(edgeMap) );
  memset( faceMap, 0, sizeof(faceMap) );
  
  dimMap[POLYGON ]      = 2;
  dimMap[TRIANGLE]      = 2;
  dimMap[QUADRILATERAL] = 2;
  dimMap[POLYHEDRON]    = 3;
  dimMap[TETRAHEDRON]   = 3;
  dimMap[HEXAHEDRON]    = 3;
  dimMap[PRISM]         = 3;
  dimMap[PYRAMID]       = 3;
  dimMap[SEPTAHEDRON]   = 3;
  
  adjMap[TRIANGLE][0] = 3;
  adjMap[TRIANGLE][1] = 3;
  
  adjMap[QUADRILATERAL][0] = 4;
  adjMap[QUADRILATERAL][1] = 4;
  
  adjMap[TETRAHEDRON][0] = 4;
  adjMap[TETRAHEDRON][1] = 6;
  adjMap[TETRAHEDRON][2] = 4;
  
  adjMap[HEXAHEDRON][0] = 8;
  adjMap[HEXAHEDRON][1] = 12;
  adjMap[HEXAHEDRON][2] = 6;
  
  adjMap[PRISM][0] = 6;
  adjMap[PRISM][1] = 9;
  adjMap[PRISM][2] = 5;
  
  adjMap[PYRAMID][0] = 5;
  adjMap[PYRAMID][1] = 8;
  adjMap[PYRAMID][2] = 5;
  
  adjMap[SEPTAHEDRON][0] = 7;
  adjMap[SEPTAHEDRON][1] = 11;
  adjMap[SEPTAHEDRON][2] = 6;  /* See description in TSTT mesh interface doc */

  int side;
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TRIANGLE-FIRST_FACE][side][0] = side;
    edgeMap[TRIANGLE-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[QUADRILATERAL-FIRST_FACE][side][0] = side;
    edgeMap[QUADRILATERAL-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[TETRAHEDRON-FIRST_FACE][side][0] = side -3 ;
    edgeMap[TETRAHEDRON-FIRST_FACE][side][1] = 3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = side;
  }
  for (side = 8; side < 12; ++side)
  {
    edgeMap[HEXAHEDRON-FIRST_FACE][side][0] = side - 4;
    edgeMap[HEXAHEDRON-FIRST_FACE][side][1] = 4+(side+1)%4;
  }
  for (side = 0; side < 3; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side;
    edgeMap[PRISM-FIRST_FACE][side][1] = (side+1)%3;
  }
  for (side = 3; side < 6; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side - 3;
    edgeMap[PRISM-FIRST_FACE][side][1] = side;
  }
  for (side = 6; side < 9; ++side)
  {
    edgeMap[PRISM-FIRST_FACE][side][0] = side-3;
    edgeMap[PRISM-FIRST_FACE][side][1] = 3+(side+1)%3;
  }
  for (side = 0; side < 4; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = (side+1)%4;
  }
  for (side = 4; side < 8; ++side)
  {
    edgeMap[PYRAMID-FIRST_FACE][side][0] = side - 4;
    edgeMap[PYRAMID-FIRST_FACE][side][1] = 4;
  }
  
  for (side = 0; side < 3; ++side)
  {
    faceMap[TETRAHEDRON-FIRST_VOL][side][0] = 3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[TETRAHEDRON-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[TETRAHEDRON-FIRST_VOL][side][3] = 3;
  }
  faceMap[TETRAHEDRON-FIRST_VOL][3][0] = 3;
  faceMap[TETRAHEDRON-FIRST_VOL][3][1] = 2;
  faceMap[TETRAHEDRON-FIRST_VOL][3][2] = 1;
  faceMap[TETRAHEDRON-FIRST_VOL][3][3] = 0;

  for (side = 0; side < 4; ++side)
  {
    faceMap[HEXAHEDRON-FIRST_VOL][side][0] = 4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][1] = side;
    faceMap[HEXAHEDRON-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][3] = 4+(side+1)%4;
    faceMap[HEXAHEDRON-FIRST_VOL][side][4] = side + 4;
  }
  faceMap[HEXAHEDRON-FIRST_VOL][4][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][4][1] = 3;
  faceMap[HEXAHEDRON-FIRST_VOL][4][2] = 2;
  faceMap[HEXAHEDRON-FIRST_VOL][4][3] = 1;
  faceMap[HEXAHEDRON-FIRST_VOL][4][4] = 0;
  faceMap[HEXAHEDRON-FIRST_VOL][5][0] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][1] = 4;
  faceMap[HEXAHEDRON-FIRST_VOL][5][2] = 5;
  faceMap[HEXAHEDRON-FIRST_VOL][5][3] = 6;
  faceMap[HEXAHEDRON-FIRST_VOL][5][4] = 7;
  
  for (side = 0; side < 4; ++side)
  {
    faceMap[PYRAMID-FIRST_VOL][side][0] = 3;
    faceMap[PYRAMID-FIRST_VOL][side][1] = side;
    faceMap[PYRAMID-FIRST_VOL][side][2] = (side+1)%4;
    faceMap[PYRAMID-FIRST_VOL][side][3] = 4;
  }
  faceMap[PYRAMID-FIRST_VOL][4][0] = 4;
  faceMap[PYRAMID-FIRST_VOL][4][1] = 3;
  faceMap[PYRAMID-FIRST_VOL][4][2] = 2;
  faceMap[PYRAMID-FIRST_VOL][4][3] = 1;
  faceMap[PYRAMID-FIRST_VOL][4][4] = 0;

  for (side = 0; side < 3; ++side)
  {
    faceMap[PRISM-FIRST_VOL][side][0] = 4;
    faceMap[PRISM-FIRST_VOL][side][1] = side;
    faceMap[PRISM-FIRST_VOL][side][2] = (side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][3] = 3+(side+1)%3;
    faceMap[PRISM-FIRST_VOL][side][4] = side + 3;
  }
  faceMap[PRISM-FIRST_VOL][3][0] = 3;
  faceMap[PRISM-FIRST_VOL][3][1] = 2;
  faceMap[PRISM-FIRST_VOL][3][2] = 1;
  faceMap[PRISM-FIRST_VOL][3][3] = 0;
  faceMap[PRISM-FIRST_VOL][4][0] = 3;
  faceMap[PRISM-FIRST_VOL][4][1] = 3;
  faceMap[PRISM-FIRST_VOL][4][2] = 4;
  faceMap[PRISM-FIRST_VOL][4][3] = 5;
}

void TopologyInfo::higher_order( EntityTopology topo,
                                     unsigned num_nodes,
                                     bool& midedge,
                                     bool& midface,
                                     bool& midvol,
                                     MsqError& err )
{
  midedge = midface = midvol = false;
  if (topo >= MIXED || num_nodes < instance.adjMap[topo][0])
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_ARG);
    return;
  }
  
  unsigned dim = instance.dimMap[topo];
  assert( num_nodes >= instance.adjMap[topo][0] );
  unsigned nodes = num_nodes - instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  if (edges && nodes >= edges)
  {
    nodes -= edges;
    midedge = true;
  }
  if (faces && nodes >= faces)
  {
    nodes -= faces;
    midface = true;
  }
  if (1 == nodes)
  {
    if (2 == dim)
    {
      nodes -= 1;
      midface = true;
    }
    else if(3 == dim)
    {
      nodes -= 1;
      midvol = true;
    }
  }
  
  if (nodes)
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_STATE);
  }
}


const unsigned*  TopologyInfo::edge_vertices( EntityTopology topo,
                                              unsigned edge, 
                                              MsqError& err)
{
  if (topo < (EntityTopology)FIRST_FACE || 
      topo > (EntityTopology)LAST_VOL || 
      edge >= edges( topo ) )
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    topo = (EntityTopology)FIRST_FACE;
    edge = 0;
  }
  
  return instance.edgeMap[topo-FIRST_FACE][edge];
}

const unsigned* TopologyInfo::face_vertices( EntityTopology topo,
                                             unsigned face,
                                             unsigned& length,
                                             MsqError& err )
{
  if (topo < (EntityTopology)FIRST_VOL || 
      topo > (EntityTopology)LAST_VOL || 
      face >= faces( topo ) )
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    topo = (EntityTopology)FIRST_VOL;
    face = 0;
  }
  
  length = instance.faceMap[topo-FIRST_VOL][face][0];
  return instance.faceMap[topo-FIRST_VOL][face] + 1;
}




const unsigned* TopologyInfo::side_vertices( EntityTopology topo,
                                             unsigned dim,
                                             unsigned side,
                                             unsigned& count_out,
                                             MsqError& err )
{
  static const unsigned all[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  const unsigned* result;
  
  if (dim != 0 && dim == dimension(topo))
  {
    count_out = corners( topo );
    result = all;
  }
  else if (dim == 1)
  {
    count_out = 2;
    result = edge_vertices( topo, side, err );
  }
  else if( dim == 2)
  {
    result = face_vertices( topo, side, count_out, err );
  } 
  else
  {
    MSQ_SETERR(err)(MsqError::INVALID_ARG);
    count_out = 0;
    result = 0;
  }    
  return result;
}
      


void TopologyInfo::side_number( EntityTopology topo,
                                    unsigned num_nodes,
                                    unsigned node_index,
                                    unsigned& side_dim_out,
                                    unsigned& side_num_out,
                                    MsqError& err )
{
  if (topo >= (EntityTopology)MIXED || num_nodes < instance.adjMap[topo][0])
  {
    MSQ_SETERR(err)("Invalid element topology", MsqError::INVALID_ARG);
    return;
  }
  
  unsigned nodes = instance.adjMap[topo][0];
  unsigned edges = instance.adjMap[topo][1];
  unsigned faces = instance.adjMap[topo][2];
  side_num_out = node_index;

  if (side_num_out < nodes)
  {
    side_dim_out = 0;
    return;
  }
  num_nodes -= nodes;
  side_num_out -= nodes;
  
  if (edges && num_nodes >= edges)
  {
    if (side_num_out < edges)
    {
      side_dim_out = 1;
      return;
    }
    num_nodes -= edges;
    side_num_out -= edges;
  }
  if (faces && num_nodes >= faces)
  {
    if (side_num_out < faces)
    {
      side_dim_out = 2;
      return;
    }
    num_nodes -= faces;
    side_num_out -= faces;
  }
  if (side_num_out == 0)
  {
    side_dim_out = instance.dimMap[topo];
    side_num_out = 0;
    return;
  }
  
  MSQ_SETERR(err)(MsqError::INVALID_ARG);
}
  
  




} //namepsace Mesquite
