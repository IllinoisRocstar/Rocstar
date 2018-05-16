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

/** \file Agent.C
 *  * Contains the base implementation of Agent. Each agent represent
 *  * one instance of a physics module (e.g. rocflo, rocflu, rocsolid or rocfrac)
 *       */

/* Author: Gengbin Zheng */

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "rocman.h"
#include "Action.h"
#include "Scheduler.h"
#include "Agent.h"
#include "Interpolate.h"
#include "basic_actions.h"

int Agent::read_by_control_handle=0;
int Agent::read_files_handle=0;
int Agent::obtain_attr_handle=0;
int Agent::write_attr_handle=0;
int Agent::write_ctrl_handle=0;

// agent's main function
void PhysicsAction::run(double t, double dt, double alpha)
{
  MAN_DEBUG(3, ("[%d] Rocstar: Agent %s::PhysicsAction run with t:%e dt:%e alpha:%e.\n", agent->comm_rank, agent->get_agent_name().c_str(), t, dt, alpha));

    // store it, gm callback needs it
  agent->current_deltatime = dt;
  agent->timestamp = t;

  if (agent->get_coupling()->get_ipc() == 1)
      agent->old_dt = dt;

  // include Interpolate::backup
  agent->run_bcinitaction(t, dt);

  // 
  if(!agent->withgm)
    COM_call_function(agent->update_handle, &t, &dt, 
		      &agent->bc_handle);
  else
    COM_call_function(agent->update_handle, &t, &dt, 
		      &agent->bc_handle, &agent->gm_handle);

  agent->current_time = t;
}

void PhysicsAction::print(FILE *f)
{
  agent->print(f);
}

char *PhysicsAction::name() 
{ 
  return (char *)agent->get_agent_name().c_str(); 
}

void AttributeBase::create(std::string bufname)
{
}

void AttributeBase::assign(std::string bufname)
{
  // check if target_window is created
//  if (COM_get_window_handle( target_window) <= 0) {
//    COM_new_window( target_window);
//    COM_use_dataitem( target_window, bufname+".all");
//  }
  COM_assertion_msg(COM_get_window_handle( target_window)>0, "ERROR: assign error.\n");
  
  COM_use_dataitem( target_window, bufname+attr);
  // when to call init_done ??????????????
}

NewAttribute::NewAttribute(Agent *ag, std::string target_window_, std::string attr_, char loc_, int type_, int ncomp_, const char * unit_) :
		AttributeBase(ag, target_window_, attr_),
		loc(loc_), type(type_), ncomp(ncomp_), unit(unit_)
{
}

void NewAttribute::create(std::string bufname)
{
  // test if already created
  if (COM_get_dataitem_handle( bufname+attr) > 0) {
    char loc_;
    int type_;
    int ncomp_;
    std::string unit_;
    COM_get_dataitem( bufname+attr, &loc_, &type_, &ncomp_, &unit_);
    if (loc_ != loc || type_ != type || ncomp_ != ncomp || unit_ != unit) {
      std::cerr << "Rocstar Error: NewAttribute::create(): Could not create " << bufname+attr << " in two different ways " << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    else
      return;
  }
  MAN_DEBUG(3, ("NewAttribute::create: %s%s.\n", bufname.c_str(), attr.c_str()));
  COM_new_dataitem( bufname+attr, loc, type, ncomp, unit);
  COM_resize_array( bufname+attr);
}

CloneAttribute::CloneAttribute(Agent *ag, int cond, 
                            std::string target_window_, std::string attr_, 
                            std::string parent_window_, std::string parent_attr_,
                            int wg_, const char *ptnname_, int val_):  
 		AttributeBase(ag, target_window_, attr_),
		parent_window(parent_window_), parent_attr(parent_attr_),
		wg(wg_), val(val_), condition(cond)
{
  ptnname = ptnname_?strdup(ptnname_):NULL;
}

void CloneAttribute::create(std::string bufname)
{
  if (parent_window=="") {   // filling
    parent_window = agent->get_surface_window();
  }
  if (COM_get_dataitem_handle( bufname+attr) > 0) {
    char loc1, loc2;
    int type1, type2;
    int ncomp1, ncomp2;
    std::string unit1, unit2;
    COM_get_dataitem( bufname+attr, &loc1, &type1, &ncomp1, &unit1);
    COM_get_dataitem( parent_window+parent_attr, &loc2, &type2, &ncomp2, &unit2);
    if (loc1 != loc2 || type1 != type2 || ncomp1 != ncomp2 || unit1 != unit2) {
      std::cerr << "Rocstar Error: CloneAttribute::create(): Could not create " << target_window+attr << " in two different ways " << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    else {
      MAN_DEBUG(3, ("CloneAttribute::create: %s%s %s%s condition:%d handle:%d SKIPPED.\n", bufname.c_str(), attr.c_str(), parent_window.c_str(), parent_attr.c_str(), condition, COM_get_dataitem_handle( parent_window+parent_attr)));
      return;
    }
  }
  MAN_DEBUG(3, ("CloneAttribute::create: %s%s %s%s condition:%d handle:%d.\n", bufname.c_str(), attr.c_str(), parent_window.c_str(), parent_attr.c_str(), condition, COM_get_dataitem_handle( parent_window+parent_attr)));
  if (condition) 
    if (COM_get_dataitem_handle( parent_window+parent_attr) <= 0)  return;
  if (ptnname) {
    if (ptnname[0] == '.')
      ptnname = strdup((agent->get_surface_window() + ptnname).c_str());
  }
  COM_clone_dataitem( bufname+attr, parent_window+parent_attr, wg, ptnname, val);
}

UseAttribute::UseAttribute(Agent *ag, std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_, const char *ptnname_, int val_):  
 		AttributeBase(ag, target_window_, attr_),
		parent_window(parent_window_), parent_attr(parent_attr_),
		wg(wg_), ptnname(ptnname_), val(val_)
{
}

void UseAttribute::create(std::string bufname)
{
  if (parent_window=="") {
    parent_window = agent->get_surface_window();
  }
  if (COM_get_dataitem_handle( bufname+attr) > 0) {
    char loc1, loc2;
    int type1, type2;
    int ncomp1, ncomp2;
    std::string unit1, unit2;
    COM_get_dataitem( bufname+attr, &loc1, &type1, &ncomp1, &unit1);
    COM_get_dataitem( parent_window+parent_attr, &loc2, &type2, &ncomp2, &unit2);
    if (loc1 != loc2 || type1 != type2 || ncomp1 != ncomp2 || unit1 != unit2) {
      std::cerr << "Rocstar Error: NewAttribute::create(): Could not create " << target_window+attr << " in two different ways " << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    else
      return;
  }
  MAN_DEBUG(3, ("UseAttribute::create: %s%s %s%s.\n", bufname.c_str(), attr.c_str(), parent_window.c_str(), parent_attr.c_str()));
  COM_use_dataitem( bufname+attr, parent_window+parent_attr, wg, ptnname, val);
}

Agent::Agent(Coupling *cp, std::string mod, std::string obj, 
	     const char *agentname, MPI_Comm com, bool wgm, bool skipio)
  : agent_name(agentname), communicator(com), coupling(cp), action(this), 
    withgm(wgm), rocmod_name( mod), mod_instance(obj), skipInputIO(skipio)
{
  MPI_Comm_rank( com, &comm_rank);
 
  if(comm_rank == 0) 
    MAN_DEBUG(3, ("Rocstar: Agent::Agent create window %s.\n", get_agent_name().c_str()));
  create_window(agent_name.c_str());

  inDir = obj;  inDir.append("/Rocin/");
  outDir = obj; outDir.append("/Rocout/");
  suffix = ".hdf";

    // set proper names for visu
  std::string physicsname = agentname;
  physicsname = physicsname+"-PhysicsAction";
  action.set_name(physicsname.c_str());

  bcInitScheduler.set_name("bcInitSched");
  gmScheduler.set_name("gmSched");

  timestamp = cp->get_control_param()->current_time;;
  old_dt = 0.0;
  dobackup = 1;
}

Agent::~Agent()
{
}

// mod_instance can be RocburnAPN, etc
void Agent::init_function_handles()
{
  // init function handle
  char fname[128];
  sprintf(fname, "%s%s", mod_instance.c_str(), ".initialize");
  init_handle = COM_get_function_handle(fname);
  COM_assertion_msg(init_handle>0, "ERROR: Agent::init_function_handles.\n");

  sprintf(fname, "%s%s", mod_instance.c_str(), ".update_solution");
  update_handle = COM_get_function_handle(fname);
  COM_assertion_msg(update_handle>0, "ERROR: Agent::init_function_handles.\n");

  sprintf(fname, "%s%s", mod_instance.c_str(), ".finalize");
  finalize_handle = COM_get_function_handle(fname);
  COM_assertion_msg(finalize_handle>0, "ERROR: Agent::init_function_handles.\n");

  // optional
  sprintf(fname, "%s%s", mod_instance.c_str(), ".pre_hdf_output");
  pre_hdf_handle = COM_get_function_handle(fname);

  // optional
  sprintf(fname, "%s%s", mod_instance.c_str(), ".post_hdf_output");
  post_hdf_handle = COM_get_function_handle(fname);

  COM_set_profiling_barrier( init_handle, communicator);
  COM_set_profiling_barrier( update_handle, communicator);
  COM_set_profiling_barrier( finalize_handle, communicator);
  if (pre_hdf_handle != -1)
    COM_set_profiling_barrier( pre_hdf_handle, communicator);
  if (post_hdf_handle != -1)
    COM_set_profiling_barrier( post_hdf_handle, communicator);

  sprintf(fname, "%s%s", mod_instance.c_str(), ".compute_integrals");
  compute_integrals_handle = COM_get_function_handle(fname);

  read_files_handle = COM_get_function_handle( "IN.read_window");
  read_by_control_handle = COM_get_function_handle( "IN.read_by_control_file");
  obtain_attr_handle = COM_get_function_handle( "IN.obtain_dataitem");
  write_attr_handle = COM_get_function_handle( "OUT.write_dataitem");
  write_ctrl_handle = COM_get_function_handle( "OUT.write_rocin_control_file");

  COM_set_profiling_barrier( read_files_handle, communicator);
  COM_set_profiling_barrier( read_by_control_handle, communicator);
  COM_set_profiling_barrier( obtain_attr_handle, communicator);
  COM_set_profiling_barrier( write_attr_handle, communicator);
}

void Agent::create_window(const char *window_name)
{
  // create window and global variables
  char global[128];
  COM_new_window(window_name);
  sprintf(global, "%s%s", window_name, ".global");
  COM_new_dataitem(global, 'w', COM_VOID, 1, "");
  COM_set_object(global, 0, this);

  // 
  COM_Type types[4];
  types[0] = COM_RAWDATA;

  char window_func[128];

  sprintf(window_func, "%s.init_callback", window_name);
  types[1] = COM_STRING;
  types[2] = COM_STRING;
  types[3] = COM_VOID;
  COM_set_member_function(window_func, (Member_func_ptr)&Agent::init_callback,
			  global, "biiI", types);
  ic_handle = COM_get_function_handle(window_func);

  types[1] = COM_DOUBLE;
  types[2] = COM_INT;
  sprintf(window_func, "%s.obtain_bc", window_name);
  COM_set_member_function(window_func, (Member_func_ptr)&Agent::obtain_bc, global, "biI", types);
  bc_handle = COM_get_function_handle(window_func);

  sprintf(window_func, "%s.obtain_gm", window_name);
  COM_set_member_function(window_func, (Member_func_ptr)&Agent::obtain_gm, global, "bi", types);
  gm_handle = COM_get_function_handle(window_func);

  COM_window_init_done(window_name);
}

double Agent::max_timestep(double t, double dt)
{
  return dt;
}

void Agent::register_new_dataitem(std::string target_window_, std::string attr_, char loc_, int type_, int ncomp_, const char * unit_)
{
  NewAttribute *newAttr = new NewAttribute(this, target_window_, attr_, loc_, type_, ncomp_, unit_);
  dataitemList.push_back(newAttr);
}

void Agent::register_clone_dataitem(int cond, std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_, const char *ptnname_, int val_)
{
  CloneAttribute *newAttr = new CloneAttribute(this, cond, target_window_, attr_, parent_window_, parent_attr_, wg_, ptnname_, val_);
  dataitemList.push_back(newAttr);
}

void Agent::register_use_dataitem(std::string target_window_, std::string attr_, std::string parent_window_, std::string parent_attr_, int wg_, const char *ptnname_, int val_)
{
  UseAttribute *newAttr = new UseAttribute(this, target_window_, attr_, parent_window_, parent_attr_, wg_, ptnname_, val_);
  dataitemList.push_back(newAttr);
}

void Agent::create_registered_dataitems(std::string tmpBuf)
{
  unsigned int n = dataitemList.size();
  for (unsigned int i=0; i<n; i++) {
    MAN_DEBUG(3, ("[%d] creating %s %s\n", comm_rank, dataitemList[i]->target_window.c_str(), dataitemList[i]->attr.c_str()));
    dataitemList[i]->create(dataitemList[i]->target_window);
  }
  COM_window_init_done( tmpBuf);
}

void Agent::create_registered_window_dataitems(std::string target_window)
{
  unsigned int n = dataitemList.size();
  for (unsigned int i=0; i<n; i++) {
    if (dataitemList[i]->target_window == target_window) {
      MAN_DEBUG(3, ("[%d] creating %s %s\n", comm_rank, dataitemList[i]->target_window.c_str(), dataitemList[i]->attr.c_str()));
      dataitemList[i]->create(dataitemList[i]->target_window);
    }
  }
}

void Agent::create_buffer_all()
{
  MAN_DEBUG(3, ("Rocstar: Agent %s::create_buffer_all called %s.\n", get_agent_name().c_str(), tmp_window.c_str()));
  COM_new_window( tmp_window);
//  COM_use_dataitem( tmp_window, surf_window+".all", 1);
  COM_use_dataitem( tmp_window, surf_window+".mesh", 1);
  COM_use_dataitem( tmp_window, surf_window+".data", 1);
  //create_registered_dataitems(tmp_window);
  create_registered_window_dataitems( tmp_window);
  COM_window_init_done( tmp_window);
}

void Agent::assign_dataitems()
{
  unsigned int n = dataitemList.size();
  for (unsigned int i=0; i<n; i++) {
    dataitemList[i]->assign(tmp_window);
  }
}

void Agent::add_icaction(Action *act)
{
  icScheduler.add_action(act);
}

// l starts from 1
void Agent::add_bcinitaction(Action *act)
{
  bcInitScheduler.add_action(act);
}

// l starts from 1
void Agent::add_bcaction(Action *act, int l)
{
  int cur =  bcScheduler.size();
  COM_assertion_msg(l>=1 && l <= cur+1, "ERROR: add_bcaction.\n");
  UserScheduler *sched;
  if (cur == l-1) {
    sched = new UserScheduler;
    char str[128];
    sprintf(str, "%s%d", "bcSched", l);
    sched->set_name(str);
    bcScheduler.push_back(sched);
  }
  else
    sched = bcScheduler[l-1];
  COM_assertion_msg(sched, "ERROR: add_bcaction with no scheduler.\n");
  sched->add_action(act);
}

void Agent::add_gmaction(Action *act)
{
  gmScheduler.add_action(act);
}

void Agent::init_module(double t, double dt)
{
 MAN_DEBUG(3, ("[%d] Rocstar: %s::init_module with current_time:%e current_deltatime:%e.\n", comm_rank, get_agent_name().c_str(), t, dt));

  current_time = initial_time = t;
  current_deltatime = dt;
  old_dt = 0.0;
}

void Agent::init_subscheduler(double t)
{
  MAN_DEBUG(3, ("Rocstar: %sAgent::init_subscheduler called.\n", get_agent_name().c_str()));

  callMethod(&Scheduler::init_actions, t);

  MAN_DEBUG(3, ("Rocstar: %sAgent::init_subscheduler done.\n", get_agent_name().c_str()));
}

void Agent::callMethod(Scheduler_voidfn1_t fn, double t)
{
  unsigned int i;

  (bcInitScheduler.*fn)( t);
  for (i=0; i<bcScheduler.size(); i++)  {
	(bcScheduler[i]->*fn)( t);
  }
  (gmScheduler.*fn)( t);
}

void Agent::init_bcactions(double t)
{
  for (unsigned int i=0; i<bcScheduler.size(); i++) 
	bcScheduler[i]->init_actions(t);
}

void Agent::init_gmactions(double t)
{
  gmScheduler.init_actions(t);
}

void Agent::schedule()
{
  icScheduler.schedule();
  bcInitScheduler.schedule();
  for (unsigned int i=0; i<bcScheduler.size(); i++)  {
	bcScheduler[i]->schedule();
  }
  gmScheduler.schedule();
}

void Agent::finalize()
{
  MAN_DEBUG(3, ("[%d] Rocstar: %s::finalize called.\n", comm_rank, get_agent_name().c_str()));

  icScheduler.finalize_actions();
  bcInitScheduler.finalize_actions();
  for (unsigned int i=0; i<bcScheduler.size(); i++) 
	bcScheduler[i]->finalize_actions();
  gmScheduler.finalize_actions();

    // at restart, we don't want windows in physics modules deleted
  if (! get_coupling()->in_restart()) 
    COM_call_function(finalize_handle);
}

// part of INITIALIZE_XXX, customized portion
void Agent::init_callback( const char *surf_win, const char *vol_win,
			   void *option)
{
  MAN_DEBUG(3, ("Rocstar: %s::init_callback called: surfwin: %s volwin: %s.\n", get_agent_name().c_str(), surf_win, vol_win));

  surf_window = surf_win;
  vol_window = vol_win;
  option_data = option;

   // creat a window buffer to inherit from surface window
   // they are like ifluid_all, iburn_all, etc
  create_buffer_all();

  if (icScheduler.isEmpty()) {
    if (comm_rank==0)
      std::cerr << "Rocstar: Warning: No actions defined for "
	      << "Initialization callback for module "
	      << rocmod_name << std::endl;
  }
  else {
    MAN_DEBUG(3, ("[%d] Rocstar: %s::init_callback called IC with current_time:%e current_deltatime:%e.\n", comm_rank, get_agent_name().c_str(), current_time, current_deltatime));
    icScheduler.set_alpha(0.0);
    icScheduler.init_actions(current_time);
    icScheduler.run_actions(current_time, current_deltatime);
    MAN_DEBUG(3, ("[%d] Rocstar: %s::init_callback called IC DONE.\n", comm_rank, get_agent_name().c_str()));
  }

  MAN_DEBUG(3, ("Rocstar: Agent::init_callback called: surfwin: %s volwin: %s DONE.\n", surf_win, vol_win));
}

/*
void Agent::run_initactions(double t, double dt)
{
  MAN_DEBUG(3, ("Rocstar: %s::run_initactions called with t:%e dt:%e.\n", get_agent_name().c_str(), t, dt));
  initScheduler.set_alpha(0.0);
  //initScheduler.init_actions(t);
  initScheduler.run_actions(t, dt);
}
*/

// called implicitly in update_solution in PhysicsAction::run()
// a is alpha time
void Agent::obtain_bc(double *a, int *l)
{
  if ( l!=NULL && (unsigned int)(*l) > bcScheduler.size() && comm_rank==0) {
    std::cerr << "Rocstar: Warning: No actions defined for "
	      << "level-" << *l << " boundary condition "
	      << "for module " << rocmod_name << std::endl;
    return;
  }
  
  int level = (l==NULL)?1:*l;
  MAN_DEBUG(3, ("[%d] Rocstar: Agent %s::obtain_bc called at level: %d t:%e.\n", comm_rank, get_agent_name().c_str(), level, timestamp))
  Scheduler *sched = get_bcScheduler(level);
  if (sched) {
    sched->set_alpha(*a);
    sched->run_actions(timestamp, current_deltatime);
  }
}

// callback. da is alpha time
void Agent::obtain_gm(double *da)
{
  if (gmScheduler.isEmpty() && comm_rank==0) {
    std::cerr << "Rocstar: Warning: No actions defined for "
	      << "Grid Motion "
	      << "for module " << rocmod_name << std::endl;
    return;
  }
  MAN_DEBUG(3, ("[%d] Rocstar: Agent::obtain_gm called with alpha: %e.\n", comm_rank, *da))
    // ???????????????
  gmScheduler.set_alpha(*da);
  gmScheduler.run_actions(0.0, current_deltatime);
}

void Agent::init_bcinitaction(double t)
{
  bcInitScheduler.init_actions(t);
}

void Agent::run_bcinitaction(double t, double dt)
{
  MAN_DEBUG(3, ("Rocstar: %s::run_bcinitaction called t=%e dt=%e.\n", get_agent_name().c_str(), t, dt));

  // backup first
  if (get_coupling()->get_ipc() <= 1 && dobackup == 1)
  for (unsigned int i=0; i<interpolateList.size(); i++)
    interpolateList[i]->backup();

  bcInitScheduler.run_actions(t, dt);
}

//  Get a string that encodes the given time. If the string is
//    xx.yyyyyy, it means 0.yyyyyy*10^xx nanoseconds. This format
//    ensures that sorting the strings alphabetically gives the right
//    ordering of time. The string is null-terminated.
void Agent::get_time_string( double t, std::string &s) {
  char chrstring[15];

  std::sprintf( chrstring, "%.5e", t*1.e10);

  s = ""; s.append( &chrstring[9]); s.append( ".");
  s.append( &chrstring[0], 1); s.append( &chrstring[2], 5);
}
 
int Agent::
read_by_control_file( double t, const std::string base, const std::string window) {

  if (skipInputIO) return 0;    // physics module does input itself

  // Read in surface window
  std::string fname = inDir; fname.append( base); fname.append("_in_");

  std::string timeLevel; get_time_string( t, timeLevel);
  fname.append( timeLevel); fname.append( ".txt");

  MAN_DEBUG(3, ("[%d] Agent %s::read_by_control_file run with file: %s t:%e.\n", comm_rank, get_agent_name().c_str(), fname.c_str(), t));

  bool ctrl_exist;
  FILE *fp = fopen( fname.c_str(), "r");
  ctrl_exist = (fp!=NULL);
  if ( !ctrl_exist && t!=0.0 ) {
    fname = outDir + base + "_in_" + timeLevel + ".txt";
    fp = fopen( fname.c_str(), "r");
    ctrl_exist = (fp!=NULL);
  }

  int comm_rank = COMMPI_Comm_rank(communicator);
  if ( ctrl_exist) {
    if (comm_rank==0 && man_verbose > 2)
      std::cout << "Rocstar: Found control file " << fname << std::endl;
    
    COM_call_function( read_by_control_handle, fname.c_str(), 
		       window.c_str(), &communicator);
    return 0;
  }
  else {
    if (comm_rank==0 && man_verbose > 2)
      std::cout << "Rocstar: Did not find control file " << fname << std::endl;
    return -1;
  }
}

void Agent::register_interpolate(InterpolateBase *ip)
{
  interpolateList.push_back(ip);
}

void Agent::
write_data_files( double t, const std::string base, const std::string attr, 
		  const char *ref) {
  int separate_out = get_coupling()->get_rocmancontrol_param()->separate_out;
  // Read in surface window
  std::string fname = outDir; 
/*
  if (separate_out)  { 
    std::ostringstream rank;
    rank << comm_rank;
    fname += rank.str(); 
    fname.append( "/"); 
  }
*/
  fname.append( base); fname.append("_");

  std::string timeLevel; get_time_string( t, timeLevel);
  fname.append( timeLevel); fname.append("_");

  MAN_DEBUG(3, ("[%d] Rocstar: Agent %s::write_data_file with file prefix %s at t:%e.\n", comm_rank, get_agent_name().c_str(), fname.c_str(), t));

  int attr_hdl = COM_get_dataitem_handle( attr);
  COM_call_function( write_attr_handle, fname.c_str(), &attr_hdl, 
		     base.c_str(), timeLevel.c_str(), ref, &communicator);
}

void Agent::
write_control_file( double t, const std::string base, const std::string window) {
  int separate_out = get_coupling()->get_rocmancontrol_param()->separate_out;
  // Read in surface window
  std::string fname;  // = base; fname.append("_");
/*
  if (separate_out)  { 
    std::ostringstream rank;
    rank << comm_rank;
    fname += rank.str(); fname.append( "/"); 
  }
*/
  fname.append( base);
  fname.append("_");

  std::string timeLevel; get_time_string( t, timeLevel);
  fname.append( timeLevel); fname.append("_");

  std::string ctrl_fname = outDir; ctrl_fname.append( base);
  // If base contains '*', then remove it
  std::string::size_type pos = ctrl_fname.find('*');
  if ( pos != std::string::npos) ctrl_fname.erase( pos);
  ctrl_fname.append("_in_"); ctrl_fname.append( timeLevel); 
  ctrl_fname.append(".txt");
  
  MAN_DEBUG(3, ("[%d] Rocstar: Agent %s::write_control_file %s with file prefix %s at t:%e.\n", comm_rank, get_agent_name().c_str(), ctrl_fname.c_str(), fname.c_str(), t));

  COM_call_function( write_ctrl_handle, window.c_str(), fname.c_str(), 
		     ctrl_fname.c_str());
}

// phase 1:   new, use Surf_win
// phase 2:        use 
void Agent::
split_surface_window( const std::string surfAll,
		      const std::string surf_i, 
		      const std::string surf_nb, 
		      const std::string surf_b, 
		      const std::string surf_ni) {
  std::string bcflag = surfAll+".bcflag";
  std::string atts = surfAll+".all";
  
  if (COM_get_window_handle(surf_i) <= 0)
  COM_new_window( surf_i);
  COM_use_dataitem( surf_i, atts, 0, bcflag.c_str(), 0);
  COM_use_dataitem( surf_i, atts, 0, bcflag.c_str(), 1);
  COM_window_init_done( surf_i);

  if (COM_get_window_handle(surf_nb) <= 0)
  COM_new_window( surf_nb); // bcflag ==0
  COM_use_dataitem( surf_nb, atts, 1, bcflag.c_str(), 0);
  COM_window_init_done( surf_nb);
  
  if (COM_get_window_handle(surf_b) <= 0)
  COM_new_window( surf_b); // bcflag ==1
  COM_use_dataitem( surf_b, atts, 1, bcflag.c_str(), 1);
  COM_window_init_done( surf_b);
  
  if (COM_get_window_handle(surf_ni) <= 0)
  COM_new_window( surf_ni); // bcflag ==2
  COM_use_dataitem( surf_ni, atts, 1,  bcflag.c_str(), 2);
  COM_window_init_done( surf_ni);
}

void Agent::init_convergence( int iPredCorr)
{
  store_solutions (iPredCorr == 1);
}

void Agent::store_solutions( int converged)
{
  int i;
  if(comm_rank == 0)
    MAN_DEBUG(3, ("[%d] Agent %s::store_solutions converged=%d.\n", comm_rank, get_agent_name().c_str(), converged));

  if ( converged) {   //  Store internal data
     for (i=0; i<pc_count; i++)
       if ( pc_hdls[0][i]>0) 
               COM_copy_dataitem( pc_hdls[1][i], pc_hdls[0][i]);
  }
  else {
     for (i=0; i<pc_count; i++)
       if ( pc_hdls[0][i]>0)
               COM_copy_dataitem( pc_hdls[0][i], pc_hdls[1][i]);
  }
}

int Agent::check_convergence_help( int vcur, int vpre, double tol, std::string str)
{
  double nrm_val=0.0, nrm_diff=0.0, ratio;

  COM_call_function( RocBlas::sub,  &vcur, &vpre, &vpre);
  COM_call_function( RocBlas::nrm2_scalar_MPI, &vcur, &nrm_val, &communicator);
  COM_call_function( RocBlas::nrm2_scalar_MPI, &vpre, &nrm_diff, &communicator);
                                                                                
  if (nrm_val != 0.) 
       ratio = nrm_diff/nrm_val;
  else
       ratio = nrm_diff;
  int result = (ratio <= tol);
                                                                                
  if (comm_rank==0 && man_verbose > 1)
         std::cout << " Convergence ratio of " << str << " is " << ratio << std::endl;
  return result;
}

void Agent::print(FILE *f)
{
  unsigned int i;

    // agent
  fprintf(f, "graph: { title: \"%s\" label: \"%s\" \n\
        status: folded \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        color :  lightyellow           \n\
        node.color     : lightblue   \n\
        node.textcolor : lightblue   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : lightblue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n", get_agent_name().c_str(), get_agent_name().c_str());

  char *bcinit = bcInitScheduler.print(f, get_agent_name().c_str());

    // main physics
  std::string main_action = get_agent_name() + "-main";
  fprintf(f, "graph: { title: \"%s\" label: \"%s\" \n\
        status: folded \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        color :  lightred           \n\
        node.color     : lightblue   \n\
        node.textcolor : black   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : lightblue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n", main_action.c_str(), main_action.c_str());

  for (i=0; i<bcScheduler.size(); i++) 
    bcScheduler[i]->print(f, get_agent_name().c_str());

  gmScheduler.print(f, get_agent_name().c_str());

  std::string main_physics = get_agent_name() + "-physics";
  fprintf(f, "node: { title:\"%s\" label:\"%s\"}\n", main_physics.c_str(), main_physics.c_str());

  fprintf(f, "}\n");
  
  if (bcinit) {
    fprintf(f, "edge: { sourcename: \"%s\" targetname: \"%s\" label: \"\"}\n", bcinit, main_action.c_str());
    free(bcinit);
  }

  fprintf(f, "}\n");
}









