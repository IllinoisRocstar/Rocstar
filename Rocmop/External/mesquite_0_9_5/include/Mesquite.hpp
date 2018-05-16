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
#ifndef MESQUITE_HPP
#define MESQUITE_HPP
#include <cstring>
#ifdef _MSC_VER
#  pragma warning ( 4 : 4786)
//#  include "mesquite_config.win.h"
//#else
//#  include "mesquite_config.h"
#endif

#ifdef MSQ_USE_OLD_C_HEADERS
#  define msq_stdc 
#else
#  define msq_stdc std
#endif

#ifdef MSQ_USE_OLD_IO_HEADERS
#  define msq_stdio
#else
#  define msq_stdio std
#endif

#ifdef MSQ_USE_OLD_STD_HEADERS
#  define msq_std
#else
#  define msq_std std
#endif


#ifdef MSQ_USE_OLD_C_INCLUDES
#  include <math.h>
#  include <float.h>
#  include <limits.h>
#else
#  include <cmath>
#  include <cfloat>
#  include <climits>
#  ifdef HAVE_CBRT
#    include <math.h>
#  endif
#endif

/*! \file Mesquite.hpp
 */

/*!
  \namespace Mesquite
  Copyright 2003 Sandia Corporation and the University of Chicago. Under
  the terms of Contract DE-AC04-94AL85000 with Sandia Corporation and
  Contract W-31-109-ENG-38 with the University of Chicago, the U.S.
  Government retains certain rights in this software.

*/
namespace Mesquite
{
  typedef int StatusCode;

  typedef double real;

  enum StatusCodeValues
  {
    MSQ_FAILURE = 0,
    MSQ_SUCCESS
  };

  enum EntityTopology
  {
    POLYGON =7,
    TRIANGLE =8,
    QUADRILATERAL =9,
    POLYHEDRON =10,
    TETRAHEDRON =11,
    HEXAHEDRON =12,
    PRISM =13,
    PYRAMID =14,
    SEPTAHEDRON =15,
    MIXED
  };

    // Version information
  const char* version_string(bool include_build_number = false);
  unsigned int major_version_number();
  unsigned int minor_version_number();
  unsigned int build_number();
  enum ReleaseType
  {
    STABLE_RELEASE,
    BETA,
    ALPHA
  };
  Mesquite::ReleaseType release_type();
  
    //GLOBAL variables
  const int MSQ_MAX_NUM_VERT_PER_ENT=8;
  const int MSQ_HIST_SIZE=7;//number of division in histogram
  static const double MSQ_SQRT_TWO = msq_stdc::sqrt(2.0);
  static const double MSQ_SQRT_THREE = msq_stdc::sqrt(3.0);
  static const double MSQ_SQRT_THREE_DIV_TWO=MSQ_SQRT_THREE/2.0;
  static const double MSQ_SQRT_THREE_INV=1.0/MSQ_SQRT_THREE;
  static const double MSQ_SQRT_TWO_INV=1.0/MSQ_SQRT_TWO;
  static const double MSQ_SQRT_TWO_DIV_SQRT_THREE=MSQ_SQRT_TWO/MSQ_SQRT_THREE;
  static const double MSQ_ONE_THIRD = 1.0 / 3.0;
  static const double MSQ_TWO_THIRDS = 2.0 / 3.0;
  static const double MSQ_3RT_2_OVER_6RT_3 = msq_stdc::pow( 2/MSQ_SQRT_THREE, MSQ_ONE_THIRD );

#ifdef UINT_MAX
  const unsigned MSQ_UINT_MAX = UINT_MAX;
#else
  const unsigned MSQ_UINT_MAX = ~(unsigned)0;
#endif

#ifdef INT_MAX
  const int MSQ_INT_MAX = INT_MAX;
#else     
  const int MSQ_INT_MAX = MSQ_UINT_MAX >> 1;
#endif

#ifdef INT_MIN
  const int MSQ_INT_MIN = INT_MIN;
#else
  const int MSQ_INT_MIN = ~MSQ_INT_MAX;
#endif

#ifdef DBL_MIN
  const double MSQ_DBL_MIN = DBL_MIN;
#else
  /* This value is rather large - DBL_MIN is normally about 2e-308 
     Put an error here to see if any platform really doesn't
     have DBL_MIN or DBL_MAX defined, and evaluate what to do then.
  */
  #error DBL_MIN not defined
  const double MSQ_DBL_MIN = 1.0E-30;
#endif
  const double MSQ_MIN = MSQ_DBL_MIN;
  
#ifdef DBL_MAX
  const double MSQ_DBL_MAX = DBL_MAX;
#else
  /* This value is rather small - DBL_MAX is normally about 2e308 
     Put an error here to see if any platform really doesn't
     have DBL_MIN or DBL_MAX defined, and evaluate what to do then.
  */
  #error DBL_MAX not defined
  const double MSQ_DBL_MAX = 1.0E30;
#endif    
  const double MSQ_MAX = MSQ_DBL_MAX;
  const double MSQ_MAX_CAP = 1.e6;

    //macro to return the min/max of a set of arguements.  The integer
    // (e.g., '2') tells how many arguements should be passed for comparison.
template <class T> inline T MSQ_MIN_2(T a, T b) { return a < b ? a : b; }
template <class T> inline T MSQ_MAX_2(T a, T b) { return a > b ? a : b; }


  // Utility functions
inline double cbrt( double d ) 
{
#ifdef HAVE_CBRT
  return ::cbrt( d );
#else
  return msq_stdc::pow( d, MSQ_ONE_THIRD );
#endif
}

inline double cbrt_sqr( double d )
{
#ifdef HAVE_CBRT
  return ::cbrt(d*d);
#else
  return msq_stdc::pow( d, MSQ_TWO_THIRDS );
#endif
}
  
}


#endif
