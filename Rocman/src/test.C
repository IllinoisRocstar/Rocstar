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

#include "rocman.h"

class ActionA: public Action {
public:
  ActionA(int n, char *at[], int i[], void *p=0, char *name=NULL): 
				Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.writes(this, "b", 1);
    sched.writes(this, "c", 1);
  }
};

class ActionB: public Action {
public:
  ActionB(int n, char *at[], int i[], void *p=0, char *name=NULL): Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.reads(this, "b", 1);
    sched.writes(this, "d", 1);
  }
};

class ActionC: public Action {
public:
  ActionC(int n, char *at[], int i[], void *p=0, char *name=NULL): Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.reads(this, "c", 1);
    sched.writes(this, "c", 1);
    sched.writes(this, "e", 1);
  }
};

class ActionD: public Action {
public:
  ActionD(int n, char *at[], int i[], void *p=0, char *name=NULL): Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.reads(this, "d", 1);
    sched.reads(this, "c", 1);
  }
};

class ActionE: public SchedulerAction
{
public:
   ActionE(Scheduler *s, int n, char *at[], int i[], void *p=0, char *name=NULL): SchedulerAction(s, n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.reads(this, "e", 1);
  }
};

class ActionF: public Action {
public:
  ActionF(int n, char *at[], int i[], void *p=0, char *name=NULL): Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.writes(this, "f", 1);
  }
};

class ActionG: public Action {
public:
  ActionG(int n, char *at[], int i[], void *p=0, char *name=NULL): Action(n, at, i, p, name) {}
  virtual void declare(Scheduler& sched) {
    sched.reads(this, "f", 1);
  }
};

int main()
{
  Scheduler * sched = new Scheduler;
 
  // root action
  ActionA *a = new ActionA(0, NULL, NULL, NULL, "A");
  sched->add_action(*a);
  a->declare(*sched);

  char *attr_c[] = { "c" };
  int idx_c[] = { 1 };
  ActionC *c = new ActionC(1, attr_c, idx_c, NULL, "C");
  sched->add_action(*c);
  c->declare(*sched);

  char *attr_d[] = { "d", "c" };
  int idx_d[] = { 1, 1 };
  ActionD *d = new ActionD(2, attr_d, idx_d, NULL, "D");
  sched->add_action(*d);
  d->declare(*sched);

  char *attr_b[] = { "b" };
  int idx_b[] = { 1 };
  ActionB *b = new ActionB(1, attr_b, idx_b, NULL, "B");
  sched->add_action(*b);
  b->declare(*sched);

  // action E is a container action that has its own sub-schduler
  // the sub scheduler contains action F and G
  Scheduler * subsched = new Scheduler;

  ActionF *f = new ActionF(0, NULL, NULL, NULL, "F");
  subsched->add_action(*f);
  f->declare(*subsched);

  char *attr_g[] = { "f" };
  int idx_g[] = { 1 };
  ActionG *g = new ActionG(1, attr_g, idx_g, NULL, "G");
  subsched->add_action(*g);
  g->declare(*subsched);

  char *attr_e[] = { "f" };
  int idx_e[] = { 1 };
  ActionE *e = new ActionE(subsched, 1, attr_e, idx_e, NULL, "E");
  sched->add_action(*e);
  e->declare(*sched);

  // all creation done now
  sched->schedule();

  // generate GDL file for visualization with iaSee tool
  sched->print("out.gdl");

  sched->init_actions(1);

  sched->run_actions(1, 0.1);

  sched->finalize_actions();
}








