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
// $Id: Action.C,v 1.15 2008/12/06 08:45:22 mtcampbe Exp $

#include <typeinfo>

#include "rocman.h"
#include "Action.h"
#include "Scheduler.h"
#include "com.h"

Action::Action(void *p, char *name): action_name(name), attr(NULL),
				     idx(NULL), count(0), usr_ptr(p) 
{
}

Action::Action(int n, const char *at[], int *id, void *p, char *name): action_name(name), count(n), usr_ptr(p) 
{
  set_attr(n, at, id);
  if (action_name == NULL) action_name = (char*)typeid(*this).name();
}

Action::Action(int n, const std::string at[], int *id, void *p, char *name): action_name(name), count(n), usr_ptr(p) 
{
  set_attr(n, at, id);
  if (action_name == NULL) action_name = (char*)typeid(*this).name();
}

Action::~Action()
{
  for (int i=0; i<count; i++) free(attr[i]);
  delete [] attr;  
  delete [] idx;
}

void Action::set_attr(int n, const char *at[], int *id)
{
  int i;
  count = n;
  attr = new char *[count];
  for (i=0; i<count; i++) attr[i] = strdup(at[i]);
  idx = new int[count];
  for (i=0; i<count; i++) idx[i] = id?id[i]:0;
}

void Action::set_attr(int n, const std::string at[], int *id)
{
  int i;
  count = n;
  attr = new char *[count];
  for (i=0; i<count; i++) attr[i] = strdup(at[i].c_str());
  idx = new int[count];
  for (i=0; i<count; i++) idx[i] = id?id[i]:0;
}

// Declare the input and ouput variables
void Action::declare(Scheduler& sched) {
  //CmiAssert(count>0);
  for ( int i=0; i<count; ++i) {
    if ( inout[i] & IN) 
      sched.reads(this, attr[i], idx[i]); 
    if ( inout[i] & OUT)
      sched.writes(this, attr[i], idx[i]);
  }  
}

void Action::print(FILE *f)
{
  fprintf(f, "node: { title:\"%s\" label:\"%s\"}\n", name(), name());
}

// Obtain dataitem handle of ith attribute
int Action::get_dataitem_handle( int i) {
  COM_assertion_msg ( attr[i]!=NULL, (std::string("Attribute \"")
			      +"\" does not exist").c_str());
  int hdl = COM_get_dataitem_handle( attr[i]);
  COM_assertion_msg ( hdl>0, (std::string("Attribute \"")+attr[i]
			      +"\" does not exist").c_str());
  return hdl;
}

int Action::get_dataitem_handle( const std::string str) {
  int hdl = COM_get_dataitem_handle( str);
  COM_assertion_msg ( hdl>0, (std::string("Attribute \"")+str
			      +"\" does not exist").c_str());
  return hdl;
}

// Obtain dataitem handle of ith attribute
int Action::get_dataitem_handle_const( int i) {
  int hdl = COM_get_dataitem_handle_const( attr[i]);
  COM_assertion_msg ( hdl>0, (std::string("Attribute \"")+attr[i]
			      +"\" does not exist").c_str());
  return hdl;
}

// SchedulerAction

SchedulerAction::~SchedulerAction() { 
  if ( sched) delete sched; 
}

void SchedulerAction::init(double t)
{
  sched->init_actions(t);
}

void SchedulerAction::schedule()  { 
  sched->schedule(); 
}

void SchedulerAction::run(double t, double dt, double alpha)
{
  //sched->set_alpha(alpha);
  sched->run_actions(t, dt);
}

void SchedulerAction::finalize()
{
  sched->finalize_actions();
}

void SchedulerAction::print(FILE *f, char *container_name)
{
  fprintf(f, "graph: { title: \"%s\" label: \"%s\" \n\
        status: folded \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        color :  red           \n\
        node.color     : black   \n\
        node.textcolor : red   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : lightblue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n", name(), name());
  sched->print(f, container_name);
  fprintf(f, "}\n");
}

void SchedulerAction::print_toposort(FILE *f) 
{ 
  fprintf(f, "( ");
  sched->print_toposort(f); 
  fprintf(f, ") ");
}


// UserSchedulerAction

UserSchedulerAction::~UserSchedulerAction() { 
  if ( sched) delete sched; 
}

void UserSchedulerAction::init(double t)
{
  sched->init_actions(t);
}

void UserSchedulerAction::schedule()  { 
  sched->schedule(); 
}

void UserSchedulerAction::run(double t, double dt, double alpha)
{
  sched->set_alpha(alpha);
  sched->run_actions(t, dt);
}

void UserSchedulerAction::finalize()
{
  sched->finalize_actions();
}

void UserSchedulerAction::print(FILE *f, char *container_name)
{
  fprintf(f, "graph: { title: \"%s\" label: \"%s\" \n\
        status: folded \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        color :  red           \n\
        node.color     : black   \n\
        node.textcolor : red   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : lightblue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n", name(), name());
  sched->print(f, container_name);
  fprintf(f, "}\n");
}

void UserSchedulerAction::print_toposort(FILE *f) 
{ 
  fprintf(f, "( ");
  sched->print_toposort(f); 
  fprintf(f, ") ");
}








