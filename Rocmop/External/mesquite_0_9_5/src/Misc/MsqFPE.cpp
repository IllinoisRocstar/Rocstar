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

#include "MsqFPE.hpp"




/* First check for BSD-style fpsetmask, which is the closest
   to a standard across unix-like OSs */

#if defined (HAVE_FPSETMASK)

#include <ieeefp.h>

bool Mesquite::MsqFPE::fpe_trap_supported()
  { return true; }

int Mesquite::MsqFPE::get_current_fpe_state()
  { return (int)fpgetmask(); }

void Mesquite::MsqFPE::set_current_fpe_state(int state)
  { fpsetmask( (fp_except)state ); }

void Mesquite::MsqFPE::enable_trap_fpe()
  { fpsetmask( fpgetmask()|FP_X_INV|FP_X_OFL|FP_X_DZ ); }





/* Next try GNU-C feenableexcept mechanism */

#elif defined (HAVE_FEENABLEEXCEPT)

#ifndef _GNU_SOURCE
#  define MSQ_SET_GNU_SOURCE
#  define _GNU_SOURCE
#endif
#include <fenv.h>
#ifdef MSQ_SET_GNU_SOURCE
#  undef _GNU_SOURCE
#  undef MSQ_SET_GNU_SOURCE
#endif


bool Mesquite::MsqFPE::fpe_trap_supported() 
  { return true; }

int Mesquite::MsqFPE::get_current_fpe_state()
  { return fegetexcept(); }

void Mesquite::MsqFPE::set_current_fpe_state(int state)
  { feenableexcept( state ); }

void Mesquite::MsqFPE::enable_trap_fpe()
{ 
  const int flags = FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW;
  feclearexcept( flags );
  feenableexcept( flags );
}





/* Next try Microsoft */

#elif defined (_MSC_VER)

#include <float.h>

bool Mesquite::MsqFPE::fpe_trap_supported() 
  { return true; }

int Mesquite::MsqFPE::get_current_fpe_state()
  { return _MCW_EM & ~_controlfp(0,0); }

void Mesquite::MsqFPE::set_current_fpe_state(int state)
  { _controlfp( state, _MCW_EM ); }

void Mesquite::MsqFPE::enable_trap_fpe()
{ 
  const int flags = _EM_ZERODIVIDE|_EM_INVALID|_EM_OVERFLOW;
  _controlfp( _controlfp(0,0) & ~flags, _MCW_EM );
}



/* Unsupported platform */
#else

bool Mesquite::MsqFPE::fpe_trap_supported() 
  { return false; }

int Mesquite::MsqFPE::get_current_fpe_state()
  { return 0; }

void Mesquite::MsqFPE::set_current_fpe_state(int )
  { }

void Mesquite::MsqFPE::enable_trap_fpe()
  { }

#endif




Mesquite::MsqFPE::MsqFPE( bool enabled ) : isEnabled( enabled )
{
  if (isEnabled) 
  {
    prevState = get_current_fpe_state();
    enable_trap_fpe();
  }
}

Mesquite::MsqFPE::~MsqFPE( )
{
  if (isEnabled)
  {
    set_current_fpe_state( prevState );
  }
}



