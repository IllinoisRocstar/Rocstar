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
#include "MesquiteError.hpp"
#include "InstructionQueue.hpp"
#include "MeshSet.hpp"
#include "TerminationCriterion.hpp"
#include "QualityAssessor.hpp"
#include "PlanarDomain.hpp"
#include "ShapeImprovementWrapper.hpp"

// algorithms
#include "MeanRatioQualityMetric.hpp"
#include "EdgeLengthQualityMetric.hpp"
#include "LPtoPTemplate.hpp"
#include "FeasibleNewton.hpp"
#include "ConjugateGradient.hpp"
#include "MsqMessage.hpp"
#include "MesqPane.h"

using namespace Mesquite;
#endif

#include "Rocmop.h"
#include "roccom.h"
#include "Pane.h"
#include "Rocblas.h"
#include "Rocmap.h"
#include "Geometric_Metrics_3.h"

MOP_BEGIN_NAMESPACE

using MAP::Rocmap;

void Rocmop::smooth_vol_mesq_wg(double initial_quality){
   if(_verb) 
    std::cout << "      Entering Rocmop::smooth_vol_mesq_wg" << std::endl;

  std::vector<COM::Pane*> allpanes;
  _wrk_window->panes(allpanes);
  int total_npanes = (int)allpanes.size();
  
  // Get the worst dihedral angle of all of the elements
  double pre_worst = initial_quality, post_worst = 180.0;
  double pre_worst_all = 180.0;
  double post_worst_all = 0.0;
  int to_cycle = 1;
  
  //std::string msg("\n  Initial shared quality = ");
  //print_mquality(msg, _is_shared_node);
  
  for(int i=0; (i<_ncycle && to_cycle); ++i){
    if(_verb > 2)
      std::cout << "        Subcycle " << i << std::endl;

    // smooth the panes via MESQUITE with ghosts.
    smooth_mesquite(allpanes,1);

    //msg = "\n  Post-mesquite quality = ";
    //print_mquality(msg, _is_shared_node);

    //msg = "\n  Post-mesquite all quality = ";
    //print_quality(msg);

    //get the worst quality of all the shared elements
    pre_worst = check_marked_elem_quality(_is_shared_elem,allpanes);
    
    if(_verb > 2)
      std::cout << "        Updating shared and ghost node positions.\n";
    Rocmap::reduce_average_on_shared_nodes(_wrk_window->attribute(COM::COM_NC));
    Rocmap::update_ghosts(_wrk_window->attribute(COM::COM_NC),
    			  _wrk_window->attribute(COM::COM_PCONN));

    //msg = "\n  Post-communicated quality = ";
    //print_mquality(msg, _is_shared_node);

    //msg = "\n  Post-communicated all quality = ";
    //print_quality(msg);

    if(_ctol != 0.0){
      post_worst = check_marked_elem_quality(_is_shared_elem,allpanes);
      if(pre_worst/post_worst > _ctol)
        to_cycle = 0;
    }
    agree_int(to_cycle, MPI_MAX);
  }
  
  if(_verb > 1) 
    std::cout << "      Exiting Rocmop::smooth_vol_mesq_wg" << std::endl;
}

MOP_END_NAMESPACE






