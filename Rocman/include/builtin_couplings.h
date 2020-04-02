/* *******************************************************************
 * Rocstar Simulation Suite                                          *
 * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
 *                                                                   *
 * Illinois Rocstar LLC                                              *
 * Champaign, IL                                                     *
 * www.illinoisrocstar.com                                           *
 * sales@illinoisrocstar.com                                         *
 *                                                                   *
 * License: See LICENSE file in top level of distribution package or *
 * http://opensource.org/licenses/NCSA                               *
 *********************************************************************/
/* *******************************************************************
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************/
// $Id: builtin_couplings.h,v 1.8 2008/12/06 08:45:22 mtcampbe Exp $

#ifndef _BUILTIN_COUPLINGS_H_
#define _BUILTIN_COUPLINGS_H_

#include "RocstarCoupling.h"

DECLARE_NEW_COUPLING_SCHEME_NO_BURN(FluidAlone);

DECLARE_NEW_COUPLING_SCHEME_WITH_BURN(FluidBurnAlone);

DECLARE_NEW_COUPLING_SCHEME_NO_BURN(SolidAlone);

DECLARE_NEW_COUPLING_SCHEME_WITH_BURN(SolidBurnAlone);

DECLARE_NEW_FULLY_COUPLING_SCHEME_NO_BURN(SolidFluidSPC);

DECLARE_NEW_FULLY_COUPLING_SCHEME_WITH_BURN(SolidFluidBurnSPC);

#endif
