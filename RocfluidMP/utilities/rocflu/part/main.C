/* *******************************************************************
 * Illinois Open Source License                                      *
 *                                                                   *
 * University of Illinois/NCSA                                       * 
 * Open Source License                                               *
 *                                                                   *
 * Copyright@2008, University of Illinois.  All rights reserved.     *
 *                                                                   *
 *  Developed by:                                                    *
 *                                                                   *
 *     Center for Simulation of Advanced Rockets                     *
 *                                                                   *
 *     University of Illinois                                        *
 *                                                                   *
 *     www.csar.uiuc.edu                                             *
 *                                                                   *
 * Permission is hereby granted, free of charge, to any person       *
 * obtaining a copy of this software and associated documentation    *
 * files (the "Software"), to deal with the Software without         *
 * restriction, including without limitation the rights to use,      *
 * copy, modify, merge, publish, distribute, sublicense, and/or      *
 * sell copies of the Software, and to permit persons to whom the    *
 * Software is furnished to do so, subject to the following          *
 * conditions:                                                       *
 *                                                                   *
 *                                                                   *
 * @ Redistributions of source code must retain the above copyright  * 
 *   notice, this list of conditions and the following disclaimers.  *
 *                                                                   * 
 * @ Redistributions in binary form must reproduce the above         *
 *   copyright notice, this list of conditions and the following     *
 *   disclaimers in the documentation and/or other materials         *
 *   provided with the distribution.                                 *
 *                                                                   *
 * @ Neither the names of the Center for Simulation of Advanced      *
 *   Rockets, the University of Illinois, nor the names of its       *
 *   contributors may be used to endorse or promote products derived * 
 *   from this Software without specific prior written permission.   *
 *                                                                   *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
 * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
 * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
 *********************************************************************
 * Please acknowledge The University of Illinois Center for          *
 * Simulation of Advanced Rockets in works and publications          *
 * resulting from this software or its derivatives.                  *
 *********************************************************************/
#include <iostream>
#include <vector>
#include <list>
#include <sstream>

using namespace std;

#include "clop.H"
#include "FC.h"

extern "C" {
  void FC_GLOBAL(rflupart,RFLUPART)(const char *,int *,long int);
}

void
Usage(const string &pn)
{
  cout << endl << "Usage: " << pn << " -c <casename> [-v 0-2]" << endl << endl
       << "   -c | --casename  :  Specifies the casename" << endl
       << "   -v | --verbosity :  Verbosity level:" << endl
       << "                       0 - Quiet" << endl
       << "                       1 - Moderately verbose" << endl
       << "                       2 - Ridiculously verbose" << endl
       << "                       3 - Exteremly verbose" << endl
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
  
  // Declare some variables for command line argument handling
  string casename;
  string sverb;
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
    if(verbosity < 0 || verbosity > 4){ // Some jerk specified a non numeric or negative
      cerr << program_name
	   << ": Invalid verbosity value.  Use -h for usage "
	   << "instructions." << endl;
      exit(1);
    }
  }
  else{ // Default verbosity
    verbosity = 1;
  }

  FC_GLOBAL(rflupart,RFLUPART)(casename.c_str(),&verbosity,casename.length());
}

