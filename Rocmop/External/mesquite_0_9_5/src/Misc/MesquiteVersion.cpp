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

#include "Mesquite.hpp"

// When changing version number/type, change these #defines, as well
// as what is returned from Mesquite::release_type().
#define MSQ_MAJOR_VERSION 0
#define MSQ_MINOR_VERSION 9
#define MSQ_BUILD_NUMBER 0
#define MSQ_VERSION_STRING "Mesquite 0.9.5"
#define MSQ_BUILD_STRING "Build Number 0"

const char* Mesquite::version_string(bool include_build_number)
{
  if (include_build_number)
    return MSQ_VERSION_STRING MSQ_BUILD_STRING;
  return MSQ_VERSION_STRING;
}

unsigned int Mesquite::major_version_number()
{
  return MSQ_MAJOR_VERSION;
}

unsigned int Mesquite::minor_version_number()
{
  return MSQ_MINOR_VERSION;
}

unsigned int Mesquite::build_number()
{
  return MSQ_BUILD_NUMBER;
}

Mesquite::ReleaseType Mesquite::release_type()
{  
  return Mesquite::ALPHA;
}

