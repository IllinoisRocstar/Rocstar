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
// $Id: Interpolate.C,v 1.39 2008/12/06 08:45:22 mtcampbe Exp $

#include "rocman.h"
#include "FluidAgent.h"
#include "BurnAgent.h"
#include "Interpolate.h"

InterpolateBase::InterpolateBase(Agent *ag, Agent *bkag, int cond): 
  Action(0, (const char**)NULL, NULL, NULL, (char *)"InterpolateBase"), 
             agent(ag), bkagent(bkag), conditional(cond)
{
  // register to Agent for create_dataitem() and backup() callback
  if (bkagent) bkagent->register_interpolate(this);

  int i;
  for (i=0; i<4; i++) attr_hdls[i] = -1;
  for (i=0; i<3; i++) bkup_hdls[i] = -1;
  on = 1;
}

void InterpolateBase::init(double t)
{
  MAN_DEBUG(3, ("%s::init called.\n", name()));
  // INIT_INTERP_HANDLES() in "man_basic.f90"
  attr_hdls[0] = get_dataitem_handle(0);
  attr_hdls[1] = COM_get_dataitem_handle_const(attr[1]);    // old
  attr_hdls[2] = COM_get_dataitem_handle(attr[2]);          // alp
  attr_hdls[3] = COM_get_dataitem_handle(attr[3]);          // grad
  for (int i=0; i<4; i++) {
      // init to 0
    if (get_io(i) == OUT) {
      double zero = 0.0;
      COM_call_function( RocBlas::copy_scalar, &zero, &attr_hdls[i]);
    }
  }
  if ( attr_hdls[1] <= 0 && attr_hdls[2] <= 0) 
       std::cout << "Rocstar Warning: Could not find dataitem " << attr[2] << std::endl;
  else if (attr_hdls[2] <= 0) { 
       std::cout << "Rocstar Error: Could not find dataitem " << attr[2] << std::endl;
       COM_assertion_msg(0, "ERROR: Abort!");
  }

  // backup handles
  bkup_hdls[0] = get_dataitem_handle_const(0);
  bkup_hdls[1] = get_dataitem_handle(1);   // old
  bkup_hdls[2] = COM_get_dataitem_handle(attr[3]);   // grad
}

void InterpolateBase::backup()
{
  MAN_DEBUG(3, ("Rocstar: InterpolateBase::backup (%s) called with dt_old:%e.\n", attr[0], agent->get_old_dt()));

  // BACKUP() in "man_basic.f90"
  if (bkup_hdls[0] > 0 && bkup_hdls[1] > 0) {
    // Compute gradient
    if ( bkup_hdls[2] > 0) {
      double dt_old = agent->get_old_dt();
      if (dt_old > 0.) {
        COM_call_function( RocBlas::sub, &bkup_hdls[0], &bkup_hdls[1], &bkup_hdls[2]);
        COM_call_function( RocBlas::div_scalar, &bkup_hdls[2], &dt_old, &bkup_hdls[2]);
      }
      else {
        double v = 0.0;
        COM_call_function( RocBlas::copy_scalar, &v, &bkup_hdls[2]);
      }
    } 
    COM_call_function( RocBlas::copy, &bkup_hdls[0], &bkup_hdls[1]);
  }
}

// INTERPOLATE_LINEAR
void InterpolateBase::extrapolate_Linear(double dt, double dt_old, double time_old, int a_old, double time_new, int a_new, double time_out, int a_out, int a_grad)
{

  //printf("extrapolate_Linear: %f %d %f %d %f %d %d\n", time_old, a_old, time_new, a_new, time_out, a_out, a_grad);  

  if (time_out == time_new) {
    COM_call_function( RocBlas::copy, &a_new, &a_out);
  }
  else if (a_old <= 0) {
    printf("Rocstar Error: Could not find the old dataitem correspond to the attribute with handle %d\n", a_out);
    COM_assertion_msg(0, "ERROR: Abort!");
  }
  else if (time_out == time_old) {
    COM_call_function( RocBlas::copy, &a_old, &a_out);
  }
  else {
      // See the interpolation section in developers' guide for the algorithm
    COM_call_function( RocBlas::sub, &a_new, &a_old, &a_out);
    double a;
    if (time_old == 0.) 
      a = time_out-1.;
    else if (time_old == -0.5) {
      if (a_grad != -100 && a_grad>0) {
        double time = (dt_old+dt)/2.0;
        COM_call_function( RocBlas::div_scalar, &a_out, &time, &a_out);
        COM_call_function( RocBlas::limit1, &a_grad, &a_out, &a_out);
        a = (time_out - 0.5)*dt;
      }
      else
        a = 2.*(time_out - 0.5)*dt/(dt_old+dt);
    }
    else if (time_old == -1) {
      COM_assertion_msg(0, "ERROR: Abort!");
    }
    else  {
      printf("Rocstar Error: Unsupported interpolation mode with old time stamp %f\n", time_old);
      COM_assertion_msg(0, "ERROR: Abort!");
    }
    COM_call_function( RocBlas::axpy_scalar, &a, &a_out, &a_new, &a_out);
  }
}


Extrapolate_Linear::Extrapolate_Linear(Agent *ag, Agent *bkag, const std::string attr, const std::string alpsuf, int cond): InterpolateBase(ag, bkag, cond)
{
  action_name = (char *)"Extrapolate_Linear";
    // INIT_INTERP_HANDLES() in "man_basic.f90"
  std::string atts[4];
  atts[0] = attr;
  atts[1] = attr + "_old";
  atts[2] = attr +  alpsuf;
  atts[3] = attr + "_grad";

  set_attr(4, atts);

  // Define I/O of the arguments.
  int io[] = {IN, IN, OUT, IN};
  set_io( 4, io);
}

void Extrapolate_Linear::init(double t)
{
  // check condition
  attr_hdls[2] = COM_get_dataitem_handle(attr[2]);          // alp
  if (conditional) {
    if (attr_hdls[2] <= 0) return;
  }

  InterpolateBase::init(t);
}

void Extrapolate_Linear::run(double t, double dt, double alpha)
{
  if (attr_hdls[2] <= 0) return;

  MAN_DEBUG(3, ("Extrapolate_Linear::run (%s) called with t:%e dt:%e alpha:%e.\n", attr[0], t, dt, alpha));

  COM_assertion_msg(alpha>=-1.e-6 && alpha <= 1+1.e-6, "ERROR: Abort!");

  int order = agent->get_coupling()->get_rocmancontrol_param()->order;
  if (order == 0 || agent->get_coupling()->new_start(t)) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }

  if (agent->get_coupling()->get_ipc() <= 1 && attr_hdls[3] <= 0) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }
  double base = 0.0;
  // linear interpolation.
  extrapolate_Linear(dt, agent->get_old_dt(), base, attr_hdls[1], base+1.0, attr_hdls[0], alpha, attr_hdls[2], attr_hdls[3]);
}

Extrapolate_Central::Extrapolate_Central(Agent *ag, Agent *bkag, const std::string attr, const std::string alpsuf, int cond): InterpolateBase(ag, bkag, cond)
{
  action_name = (char *)"Extrapolate_Central";
    // INIT_INTERP_HANDLES() in "man_basic.f90"
  std::string atts[4];
  atts[0] = attr;
  atts[1] = attr + "_old";
  atts[2] = attr +  alpsuf;
  atts[3] = attr + "_grad";

  set_attr(4, atts);

  // Define I/O of the arguments.
  int io[] = {IN, IN, OUT, IN};
  set_io( 4, io);

/*
    // register dataitems
  std::string surf_win = agent->get_surface_window();
  std::string::size_type pos = attr.find(".");
  std::string wname, aname;
  if ( pos == std::string::npos) {
    COM_assertion_msg(0, "ERROR: Here!");
  }
  else {
    wname = attr.substr( 0, pos);
    aname = attr.substr( pos, attr.size()-pos);
  }

  agent->register_use_dataitem( wname, aname, ((BurnAgent*)agent)->parentWin, aname);   // attr
  agent->register_clone_dataitem( conditional, wname, aname+"_old", surf_win, aname+"_alp");   // attr_old
  agent->register_clone_dataitem( conditional, wname, aname+"_grad", surf_win, aname+"_alp");  // attr_grad
*/
}

void Extrapolate_Central::init(double t)
{
  // check condition
  // first check if window exists
  std::string dataitem = attr[0];
  std::string::size_type pos = dataitem.find(".");
  COM_assertion_msg(pos != std::string::npos, "ERROR: Here!");
  std::string wname = dataitem.substr( 0, pos);
  if (COM_get_window_handle(wname) <= 0) return;

    // check if dataitem exists
  attr_hdls[2] = COM_get_dataitem_handle(attr[2]);          // alp
  if (conditional) {
    if (attr_hdls[2] <= 0) return;
  }

  InterpolateBase::init(t);
}

void Extrapolate_Central::run(double t, double dt, double alpha)
{
  if (attr_hdls[2] <= 0) return;

  MAN_DEBUG(3, ("Extrapolate_Central::run (%s) called with t:%e dt:%e alpha:%e.\n", attr[0], t, dt, alpha));

  COM_assertion_msg(alpha>=0.0, "ERROR: Extrapolate_Central called with invalid alpha!");

  COM_assertion_msg(alpha>=-1.e-6 && alpha <= 1+1.e-6, "ERROR: Abort!");

  int order = agent->get_coupling()->get_rocmancontrol_param()->order;
  if (order == 0 || agent->get_coupling()->new_start(t)) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }

  if (agent->get_coupling()->get_ipc() <= 1 && attr_hdls[3] <= 0) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }
  double base = -0.5;
  // linear interpolation.
  extrapolate_Linear(dt, agent->get_old_dt(), base, attr_hdls[1], base+1.0, attr_hdls[0], alpha, attr_hdls[2], attr_hdls[3]);
}


Interpolate_Linear::Interpolate_Linear(Agent *ag, Agent *bkag, const std::string attr, const std::string alpsuf): InterpolateBase(ag, bkag)
{
  action_name = (char *)"Interpolate_Linear";
    // INIT_INTERP_HANDLES() in "man_basic.f90"
  std::string atts[4];
  atts[0] = attr;
  atts[1] = attr + "_old";
  atts[2] = attr +  alpsuf;
  atts[3] = attr + "_grad";

  set_attr(4, atts);

  // Define I/O of the arguments.
  int io[] = {IN, IN, OUT, IN};
  set_io( 4, io);
}

void Interpolate_Linear::run(double t, double dt, double alpha)
{
  if (attr_hdls[2] <= 0) return;

  MAN_DEBUG(3, ("Interpolate_Linear::run (%s) called with t:%e dt:%e alpha:%e.\n", attr[0], t, dt, alpha));

  COM_assertion_msg(alpha>=0.0, "ERROR: Interpolate_Linear called with invalid alpha!");

  COM_assertion_msg(alpha>=-1.e-6 && alpha <= 1+1.e-6, "ERROR: Abort!");
  
  int order = agent->get_coupling()->get_rocmancontrol_param()->order;
  if (order == 0 || agent->get_coupling()->new_start(t)) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }

  double base = 0.0;
  // linear interpolation.
  extrapolate_Linear(dt, agent->get_old_dt(), base, attr_hdls[1], base+1.0, attr_hdls[0], alpha, attr_hdls[2], attr_hdls[3]);
}


Interpolate_Constant::Interpolate_Constant(Agent *ag, Agent *bkag, const std::string attr, const std::string alpsuf): InterpolateBase(ag, bkag)
{
  action_name = (char *)"Interpolate_Constant";

  std::string atts[4];
  atts[0] = attr;
  atts[1] = attr + "_old";
  atts[2] = attr +  alpsuf;
  atts[3] = attr + "_grad";
  set_attr(4, atts);

  // Define I/O of the arguments.
  int io[] = {IN, IN, OUT, IN};
  set_io( 4, io);

    // register dataitems
  std::string::size_type pos = attr.find(".");
  std::string wname, aname;
  if ( pos == std::string::npos) {
    COM_assertion_msg(0, "ERROR: Here!");
  }
  else {
    wname = attr.substr( 0, pos);
    aname = attr.substr( pos, attr.size()-pos);
  }
  std::string surf_win = agent->get_surface_window();

}

void Interpolate_Constant::run(double t_dummy, double dt_dummy, double alpha)
{
  if (attr_hdls[2] <= 0) return;

  MAN_DEBUG(3, ("Interpolate_Constant::run (bc) called with alpha:%e.\n", alpha));

  COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
}

Interpolate_Central::Interpolate_Central(Agent *ag, Agent *bkag, const std::string attr, const std::string alpsuf): InterpolateBase(ag, bkag)
{
  action_name = (char *)"Interpolate_Central";
    // INIT_INTERP_HANDLES() in "man_basic.f90"
  std::string atts[4];
  atts[0] = attr;
  atts[1] = attr + "_old";
  atts[2] = attr +  alpsuf;
  atts[3] = attr + "_grad";

  set_attr(4, atts);

  // Define I/O of the arguments.
  int io[] = {IN, IN, OUT, IN};
  set_io( 4, io);

/*
    // register dataitems
  std::string surf_win = agent->get_surface_window();
  std::string::size_type pos = attr.find(".");
  std::string wname, aname;
  if ( pos == std::string::npos) {
    COM_assertion_msg(0, "ERROR: Here!");
  }
  else {
    wname = attr.substr( 0, pos);
    aname = attr.substr( pos, attr.size()-pos);
  }

  if (attr_action[0] == 'n')
    agent->register_new_dataitem( wname, aname, 'e', COM_DOUBLE, 1, "kg/(m^2 s)");

  if (attr_action[1] == 'n')
    agent->register_new_dataitem( wname, aname+"_old", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  else if (attr_action[1] == 'c')
    agent->register_clone_dataitem( 0, wname, aname+"_old", surf_win, aname);   // attr_old

  if (attr_action[2] == 'c')
    agent->register_clone_dataitem( 0, wname, aname+alpsuf, surf_win, aname);  // attr_grad

  if (attr_action[3] == 'n')
    agent->register_new_dataitem( wname, aname+"_grad", 'e', COM_DOUBLE, 1, "kg/(m^2 s)");
  else if (attr_action[3] == 'c')
    agent->register_clone_dataitem( 0, wname, aname+"_grad", surf_win, aname);  // attr_grad
*/
}

void Interpolate_Central::run(double t, double dt, double alpha)
{
  if (attr_hdls[2] <= 0) return;

  MAN_DEBUG(3, ("Interpolate_Central::run (%s) called with t:%e dt:%e alpha:%e.\n", attr[0], t, dt, alpha));

  COM_assertion_msg(alpha>=0.0, "ERROR: Extrapolate_Central called with invalid alpha!");

  COM_assertion_msg(alpha>=-1.e-6 && alpha <= 1+1.e-6, "ERROR: Abort!");

  int order = agent->get_coupling()->get_rocmancontrol_param()->order;
  if (order == 0 || agent->get_coupling()->new_start(t)) {
    COM_call_function( RocBlas::copy, &attr_hdls[0], &attr_hdls[2]);
    return;
  }

  double base = -0.5;

  // linear interpolation.
  extrapolate_Linear(dt, agent->get_old_dt(), base, attr_hdls[1], base+1.0, attr_hdls[0], alpha, attr_hdls[2], attr_hdls[3]);
}








