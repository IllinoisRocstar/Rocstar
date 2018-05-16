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
// $Id: derived_couplings.h,v 1.17 2008/12/06 08:45:22 mtcampbe Exp $

#ifndef _DERIVED_COUPLINGS_H_
#define _DERIVED_COUPLINGS_H_

#include "Coupling.h"

#define DECLARE_NEW_COUPLING_SCHEME( New_scheme)			\
  class New_scheme : public Coupling {					\
  public:								\
    New_scheme(MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);		\
    New_scheme(const char *, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);\
    New_scheme(const char *, const char *, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);		\
    New_scheme(const char *, const char *, const char*, MPI_Comm com, Control_parameters *p, const RocmanControl_parameters *mp);	\
  };

DECLARE_NEW_COUPLING_SCHEME( FluidSolidISS);

DECLARE_NEW_FULLY_COUPLING_SCHEME( SolidFluidBurnEnergySPC);

#endif






