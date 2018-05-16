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
// $Id: Scheduler.C,v 1.14 2008/12/06 08:45:22 mtcampbe Exp $

#include <assert.h>

#include "rocman.h"
#include "com.h"
#include "Action.h"
#include "Scheduler.h"

/**
   Scheduler
*/

Scheduler::Scheduler(): alphaT(-1), scheduled(0), inited(0)
{
  scheduler_name = "Scheduler";
}

void Scheduler::add_action(Action *action)
{
  ActionItem *aitem = new ActionItem(action);
  actions.push_back(aitem);
  // call back action
  action->declare(*this);
}

void Scheduler::reads(Action *a, const char *attr, int idx)
{
  // FIXME: locate this action, doing linear search now
  ActionList::iterator aitem;
  for (aitem = actions.begin(); aitem != actions.end(); aitem++)
    if ((*aitem)->action() == a) break;
  if (aitem == actions.end()) {
    printf("ERROR: action '%s' not registered to scheduler. \n", a->name());
    exit(1);
  }
  (*aitem)->read_attr.push_back(strdup(attr));     // copy the string
  (*aitem)->read_idx.push_back(idx);
  (*aitem)->input.push_back(NULL);
}

void Scheduler::writes(Action *a, const char *attr, int idx)
{
  // locate this action
  ActionList::iterator aitem = actions.begin();
  for (; aitem != actions.end(); aitem++)
    if ((*aitem)->action() == a) break;
  if (aitem == actions.end()) {
    printf("ERROR: action '%s' not registered to scheduler. \n", a->name());
    COM_assertion(false);
  }
  (*aitem)->write_attr.push_back(strdup(attr));   // copy the string
  (*aitem)->write_idx.push_back(idx);
  (*aitem)->output.push_back(NULL);
}

void Scheduler::buildDDG()
{
#if 0
   ActionList pool;
   pool = actions;

   while (pool.size()!=0) {
     // find one leaf with no output or output already fullfiled
     ActionList::iterator act = pool.end()-1;
     for (; act>=pool.begin(); act--) 
     {
       if ((*act)->fullfilled()== 1) {
//printf("FOUND LEAF: %s\n", (*act)->name());
         break;
       }
     }
     if (act < pool.begin()) {   // error with loop
       printf("Loop detected!\n");
       COM_assertion(false);
     }

     ActionItem *item = *act;
     pool.erase(act);		// don't use iterator act afterwards

     if (item->n_read() == 0) {    // this is a tree root
       roots.push_back(item);
     }
     else
       for (unsigned int i=0; i<item->n_read(); i++) {
         if (item->input[i] == NULL) {
           // find matching input actions for item from pool
           ActionList::iterator a = pool.end()-1;
           for (; a>=pool.begin(); a--)  {
             int inIdx = (*a)->hasOutput(item->read_attr[i], item->read_idx[i]);
	     if (inIdx != -1) {     // set up double link
//printf("updated output of %s for %s %d\n", (*a)->name(), item->read_attr[i], item->read_idx[i]);
               item->input[i] = *a;
               (*a)->output[inIdx] = item;
	       break;			// only do once
             }
           }
           if (a < pool.begin()) { // a read without any matching write
             printf("Error: can not find matching input for action %s of attr:%s index:%d. \n", (*a)->name(), item->read_attr[i], item->read_idx[i]);
           }
         }       // end of for i
     }
   }    // end of while
#else
   //  make links
   ActionList pool;
   pool = actions;

   unsigned int idx = 1;
   while (idx<pool.size()) {
     
     ActionItem *item = pool[idx]; 
     // search all input
     for (unsigned int i=0; i<item->n_read(); i++) {
       for (int j=idx-1; j>=0; j--)  {
         ActionItem *aitem = pool[j]; 
         int inIdx = aitem->hasOutput(item->read_attr[i], item->read_idx[i]);
         if (inIdx != -1) {
           item->input[i] = aitem;
           aitem->output[inIdx] = item;
         }
       }
     }
     idx ++;
   }
   //  identify root nodes
   for (idx=0; idx<pool.size(); idx++) {
     ActionItem *item = pool[idx]; 
     int isroot = 1;
     for (unsigned int i=0; i<item->n_read(); i++) 
       if (item->input[i] != NULL) {
         isroot = 0;
         break;
       }
     if (isroot) roots.push_back(item);
   }
#endif
}

void Scheduler::schedule()
{
   MAN_DEBUG(1, ("Scheduler::schedule called.\n"));
  COM_assertion_msg(!scheduled, "ERROR: Scheduler has already been scheduled.\n");

  // schedule all sub-schedulers
  for (unsigned int i=0; i<actions.size(); i++)
    actions[i]->action()->schedule();

  // debugging
  printActions();

  buildDDG();

  sanityCheck();

  topological_sort();

  printf("Topological sort:\n");
  print_toposort(stdout);
  printf("\n");

  scheduled = 1;
}

void Scheduler::printActions()
{
  for (unsigned int i=0; i<actions.size(); i++)
    actions[i]->print(stdout);
  printf("\n");
}

void Scheduler::print_helper(FILE *f, ActionItem *aitem)
{
  unsigned int i;
  if (aitem->print_flag == 1) return;
  aitem->print_flag = 1;
  aitem->action()->print(f);
  for (i=0; i<aitem->n_read(); i++) {
    if (aitem->input[i])
      fprintf(f, "edge: { sourcename: \"%s\" targetname: \"%s\" label: \"%s,%d\"}\n", aitem->input[i]->name(), aitem->name(), aitem->read_attr[i], aitem->read_idx[i]);
  }
  for (i=0; i<aitem->n_write(); i++)
    if ( aitem->output[i]) print_helper(f, aitem->output[i]);
}

char * Scheduler::print(FILE *f, const char *container_name)
{
  if (actions.size() == 0) return NULL;

  std::string sched_name = container_name;
  sched_name = sched_name + "-" + name();
  fprintf(f, "graph: { title: \"%s\" label: \"%s\" \n\
        status: folded \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        color :  white           \n\
        node.color     : lightblue   \n\
        node.textcolor : black   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : lightblue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n", sched_name.c_str(), sched_name.c_str());

  unsigned int i;
  for (i=0; i<actions.size(); i++)
    actions[i]->print_flag = 0;

  for (i=0; i<roots.size(); i++) {
    print_helper(f, roots[i]);
  }

  fprintf(f, "}\n");

  return strdup(sched_name.c_str());
}

// print in GDL
void Scheduler::print(const char *fname)
{
  FILE *f = fopen(fname, "w");

  fprintf(f, "graph: { title: \"DependenceTree\" \n\
        display_edge_labels: yes \n\
        layoutalgorithm: tree   \n\
        scaling: maxspect   \n\
        node.color     : green   \n\
        node.textcolor : black   \n\
        node.bordercolor: black \n\
        node.borderwidth: 1    \n\
        edge.color     : blue   \n\
        edge.arrowsize : 7   \n\
        edge.thickness : 2   \n\
        edge.fontname:\"helvO08\"  \n\
        node.label: \"no type\" \n");

  print(f, "scheduler");

  fprintf(f, "} \n");
  fclose(f);
}

void Scheduler::topological_sort()
{
  ActionList pool;
  pool = actions;
  while (1) {
    ActionList::iterator aitem = pool.begin();
    for (; aitem!=pool.end(); aitem++) {
      if ((*aitem)->n_read() == 0) { break; }
      // check every input
      int flag = 1;
      for (unsigned int i=0; i<(*aitem)->n_read(); i++) {
        ActionItem *in = (*aitem)->input[i];
        // search in sort
        ActionList::iterator s;
        for (s = sort.begin(); s!=sort.end(); ++s)
          if ((*s)->action() == in->action()) break;
        if (s == sort.end()) { flag = 0; break; }
      }
      if (flag) break;
    }
    if (aitem != pool.end()) {
      sort.push_back(*aitem);
      pool.erase(aitem);
    }
    else
      break;
  }
  
  if (pool.size() != 0) {
    printf("ERROR in sorting!\n");
    exit (1);
  }
}

void Scheduler::print_toposort(FILE *f)
{
  for (unsigned int i=0; i<sort.size(); i++)  sort[i]->action()->print_toposort(f);
}

void Scheduler::sanityCheck()
{
  // make sure all input and output are not null
  ActionList::iterator aitem;
  for (aitem = actions.begin(); aitem!=actions.end(); aitem++) 
  {
    ActionItem *item = *aitem;
    unsigned int i;
    COM_assertion_msg(item->n_read() == item->input.size(), "ERROR: Scheduler::sanityCheck failed.\n");
    for (i=0; i<item->n_read(); i++) 
      COM_assertion_msg(item->input[i] != NULL, "ERROR: Scheduler::sanityCheck failed at 2.\n");
    COM_assertion_msg(item->n_write() == item->output.size(), "ERROR: Scheduler::sanityCheck failed at 3.\n");
    for (i=0; i<item->n_write(); i++) 
      COM_assertion_msg(item->output[i] != NULL, "ERROR: Scheduler::sanityCheck failed at 4.\n");
  }
}

void Scheduler::init_actions(double t)
{
  COM_assertion_msg(scheduled, "RROR: Scheduler has not been scheduled.\n");
  if (inited) return;
  inited = 1;

  // do at sorted order
  ActionList::iterator aitem;
  for (aitem = sort.begin(); aitem!=sort.end(); aitem++) 
  {
    ActionItem *item = *aitem;
    item->action()->init(t);
  }
}

void Scheduler::run_actions(double t, double dt)
{
  COM_assertion_msg(scheduled, "RROR: Scheduler has not been scheduled when calling run_actions.\n");
  // do at sorted order
  ActionList::iterator aitem;
  for (aitem = sort.begin(); aitem!=sort.end(); aitem++) 
  {
    ActionItem *item = *aitem;
    item->action()->run(t, dt, alphaT);
  }
}

void Scheduler::finalize_actions()
{
  COM_assertion_msg(scheduled, "ERROR: Scheduler has not been scheduled when calling finalize_actions.\n");
  if (sort.size()==0) return;
  // do at reversed order
  ActionList::iterator aitem;
  for (aitem = sort.end()-1; aitem>=sort.begin(); aitem--) 
  {
    ActionItem *item = *aitem;
    item->action()->finalize();
  }
}



// return true if all output actions satisfied
int Scheduler::ActionItem::fullfilled()
{
  for (unsigned int i=0; i<n_write(); i++)
    if (output[i] == NULL) return 0;
  return 1;
}

int Scheduler::ActionItem::hasInput(const char *attr, int idx)
{
  for (unsigned int i=0; i<n_read(); i++) 
    if (strcasecmp(attr, read_attr[i]) == 0 && idx == read_idx[i]) return i;
  return -1;
}

int Scheduler::ActionItem::hasOutput(const char *attr, int idx)
{
  for (unsigned int i=0; i<n_write(); i++) 
    if (strcasecmp(attr, write_attr[i]) == 0 && idx == write_idx[i]) return i;
  return -1;
}


void Scheduler::ActionItem::print(FILE *f)
{
  unsigned int i;
  fprintf(f, "=========== Action %s =============\n", name());
  for (i=0; i<n_read(); i++) 
    printf("reads: %s %d \n", read_attr[i], read_idx[i]);
  for (i=0; i<n_write(); i++) 
    printf("writes: %s %d \n", write_attr[i], write_idx[i]);
}


// scheduler is done by user specified add action sequence
void UserScheduler::add_action(Action *action)
{
  Scheduler::add_action(action);
  if (roots.size() == 0) roots.push_back(actions[actions.size()-1]);
  sort.push_back(actions[actions.size()-1]);
}

void UserScheduler::schedule() 
{ 
/*
  for (unsigned int i=0; i<actions.size(); i++)  {
    if (i+1<actions.size()) {
      actions[i]->output.push_back(actions[i+1]);
      actions[i+1]->input.push_back(actions[i]);
    }
  }  
*/
  buildDDG();
  scheduled = 1; 
}






