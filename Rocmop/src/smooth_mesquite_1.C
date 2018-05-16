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
#ifdef MESQUITE
#define USE_STD_INCLUDES 1
#define USE_C_PREFIX_INCLUDES 1
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprovementWrapper.hpp"

// algorithms
#include "MeanRatioFunctions.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MesqPane_95.h"

using namespace Mesquite;
#endif

#include "Rocmop_1.h"
#include "com.h"
#include "Pane.hpp"
#include "Rocblas.h"
#include "Rocmap.h"
#include "Geometric_Metrics_3.h"

MOP_BEGIN_NAMESPACE

using MAP::Rocmap;

void Rocmop::smooth_vol_mesq_wg(){
  
  print_legible(1,"      Entering Rocmop::smooth_vol_mesq_wg");

  print_legible(2,"        Updating ghost node positions.");
  Rocmap::update_ghosts(_buf_window->dataitem(COM::COM_NC),
			_buf_window->dataitem(COM::COM_PCONN));

  std::vector<COM::Pane*> allpanes;
  _buf_window->panes(allpanes);
  
  // smooth the panes via MESQUITE with ghosts.
  smooth_mesquite(allpanes,1);
  
  print_legible(2,"        Updating shared and ghost node positions.");
  Rocmap::reduce_average_on_shared_nodes(_buf_window->dataitem(COM::COM_NC));
  Rocmap::update_ghosts(_buf_window->dataitem(COM::COM_NC),
			_buf_window->dataitem(COM::COM_PCONN));
  
  print_legible(1,"      Exiting Rocmop::smooth_vol_mesq_wg");
}

MOP_END_NAMESPACE






