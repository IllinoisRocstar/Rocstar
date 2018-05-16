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

#include <cstdlib>

// Utility Functions
vector<string> Vectize(const char **);
vector<string> Vectize(const char **,int);
const string stripdir(const string &);

class
clop
{
protected:
  string       _longop;
  char         _shortop;
public:
  clop()
    : _shortop((char)0)
  {};
  clop(const clop &op)
    : _longop(op._longop),_shortop(op._shortop)
  {};
  clop(const string &lo,const char &so)
    : _longop(lo),_shortop(so)
  {};
  bool operator==(const clop &op) const
  {
    return(_longop  == op._longop  && 
           _shortop == op._shortop);
  };
  bool operator==(const string &s) const
  {
    return(s == _longop);
  };
  clop &operator=(const clop &op)
  {
    _longop = op._longop;
    _shortop = op._shortop;
    return(*this);
  };
  string &longop()
  {
    return(_longop);
  };
  const string &longop() const
  {
    return(_longop);
  };
  char &shortop()
  {
    return(_shortop);
  };
  const char &shortop() const
  {
    return(_shortop);
  };
};

bool
GetOp(const string &ops,const vector<string> &args);

bool
GetOp(const string &ops,string &rv,const vector<string> &args);

bool
CheckOps(const vector<string> &args,string &badop);

void
AddOp(const string &lo,const char &so);

// ADD UNDERSCORE IF NEEDED
#ifdef FORTRAN_UNDERSCORE
#  define F_FUNC(x)  x ## _
#else
#  define F_FUNC(x)  x
#endif

#ifdef FORTRAN_LOWERCASE 

// FUNCTIONS FOR LOWERCASE FORTRAN SYMBOLS
#define RFLUCONV            F_FUNC(rfluconv)
#define RFLUINIT            F_FUNC(rfluinit)
#define RFLUMAP             F_FUNC(rflumap)
#define RFLUMP              F_FUNC(rflump)
#define RFLUPART            F_FUNC(rflupart)
#define RFLUPICK            F_FUNC(rflupick)
#define RFLUPOST            F_FUNC(rflupost)

#else 

// FUNCTIONS FOR UPPERCASE FORTRAN SYMBOLS
#define RFLUCONV            F_FUNC(RFLUCONV)
#define RFLUINIT            F_FUNC(RFLUINIT)
#define RFLUMAP             F_FUNC(RFLUMAP)
#define RFLUMP              F_FUNC(RFLUMP)
#define RFLUPART            F_FUNC(RFLUPART)
#define RFLUPICK            F_FUNC(RFLUPICK)
#define RFLUPOST            F_FUNC(RFLUPOST)

#endif







