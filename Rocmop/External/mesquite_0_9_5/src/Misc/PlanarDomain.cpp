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
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov,
    kraftche@cae.wisc.edu   
   
  ***************************************************************** */
#include "PlanarDomain.hpp"

void Mesquite::PlanarDomain::set_plane( const Mesquite::Vector3D& normal, 
                                        const Mesquite::Vector3D& point)
{
  mNormal = normal;
  mNormal.normalize();
  mCoeff = -(mNormal % point);
}

void Mesquite::PlanarDomain::snap_to(Mesquite::Mesh::EntityHandle
                                       entity_handle,
                                     Vector3D &coordinate) const
{
  coordinate -= mNormal * ( mNormal % coordinate + mCoeff );
}


void Mesquite::PlanarDomain::normal_at(
  Mesquite::Mesh::EntityHandle /*entity_handle*/,
  Mesquite::Vector3D &coordinate) const
{
  coordinate = mNormal;
}


void Mesquite::PlanarDomain::normal_at( Mesquite::Mesh::EntityHandle ,
                                        Vector3D coords[],
                                        unsigned count,
                                        Mesquite::MsqError& ) const
{
  for (unsigned i = 0; i < count; ++i)
    coords[i] = mNormal;
}

void Mesquite::PlanarDomain::closest_point( Mesquite::Mesh::EntityHandle ,
                                            const Mesquite::Vector3D& position,
                                            Mesquite::Vector3D& closest,
                                            Mesquite::Vector3D& normal,
                                            Mesquite::MsqError& ) const
{
  normal = mNormal;
  closest = position - mNormal * (mNormal % position + mCoeff);
}
