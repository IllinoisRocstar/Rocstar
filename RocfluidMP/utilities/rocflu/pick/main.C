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
#include <iostream>
#include <vector>
#include <list>
#include <sstream>

using namespace std;

#include "clop.H"

extern "C" {
  void RFLUPICK(const char *,const char *,int *,long int,long int);
}

void
Usage(const string &pn)
{
  cout << endl << "Usage: " << pn << " -c <casename> -s <stamp> [-v 0-2]" 
       << endl << endl
       << "   -c | --casename  :  Specifies the casename" << endl
       << "   -s | --stamp     :  Iteration or time stamp" << endl
       << "   -v | --verbosity :  Verbosity level:" << endl
       << "                       0 - Nothing but warnings" << endl
       << "                       1 - Normal verbosity" << endl
       << "                       2 - More verbose" << endl
       << "                       3 - Ridiculously verbose" << endl
       << endl;
}

int
main(int argc,char *argv[])
{
  // Get the commandline into a string vector - it's easier
  // to deal with that way.
  vector<string> args = Vectize((const char **)argv,argc);

  // Get the name of the executable by stripping off any leading
  // directory names
  string program_name(stripdir(args[0]));

  // Specify the allowable options to the program
  AddOp("casename",'c');
  AddOp("verbosity",'v');
  AddOp("help",'h');
  AddOp("stamp",'s');

  // Declare some variables for command line argument handling
  string casename;
  string sverb;
  string stamp;
  int verbosity;
  bool help;
  bool isset;
  
  // See if the help option is specified, if so give'm the usage text
  if(help = GetOp("help",args)){
    Usage(program_name);
    exit(0);
  }

  // Process casename option, if it's not set then fail
  if(GetOp("casename",casename,args)){
    if(casename.empty()){ // casename was empty
      cerr << program_name 
	   << ": Expected casename after casename option."
	   << " Use -h for usage instructions."
	   << endl;
      exit(1);
    }
  }
  else{ // option not specified (but it's required!)
    cerr << program_name
	 << ": Missing required casename option."
	 << "  Use -h for usage instructions."
	 << endl;
    exit(1);
  }      

  // Process stamp option, if it's not set then fail
  if(GetOp("stamp",stamp,args)){
    if(casename.empty()){ // stamp was empty
      cerr << program_name 
	   << ": Expected stamp after stamp option."
	   << " Use -h for usage instructions."
	   << endl;
      exit(1);
    }
  }
  else{ // option not specified (but it's required!)
    cerr << program_name
	 << ": Missing required stamp option."
	 << "  Use -h for usage instructions."
	 << endl;
    exit(1);
  }      

  // Process verbosity option
  if(GetOp("verbosity",sverb,args)){
    if(sverb.empty()){
      cerr << program_name
	   << ": Expected verbosity level. "
	   << "Use -h for usage instructions." << endl;
      exit(1);
    }
    istringstream Istr(sverb);
    Istr >> verbosity;
    if(verbosity < 0 || verbosity > 3){ // Some jerk specified a non numeric or negative
      cerr << program_name
	   << ": Invalid verbosity value.  Use -h for usage "
	   << "instructions." << endl;
      exit(1);
    }
  }
  else{ // Default verbosity
    verbosity = 1;
  }

  RFLUPICK(casename.c_str(),stamp.c_str(),&verbosity,casename.length(),stamp.length());
}






