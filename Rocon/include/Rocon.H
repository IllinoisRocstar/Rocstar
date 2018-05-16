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

#ifndef _ROCON_H_
#define _ROCON_H_

#include "com.h"
#include "com_devel.hpp"
#include "BndSurface.h"

class Rocon : public COM_Object {
public:
  /// Default constructor.
  Rocon(): _TOL(1E-14){};
  virtual ~Rocon(){};
  /// Initialize Rocon from the given file
  void init_from_file( const char *inp,const int *ndiv);

  /// Initialize Rocon with given mesh.
  void initialize( const COM::DataItem *pmesh,const int *ndiv);

  /// Displace the points in pmesh by disp and calcuate the intersections
  /// (if any) with the constraint mesh set in "initialize".  The int IO 
  /// attribute "constr" indicates (by a 1 or 0) which nodes are already 
  /// on the constraint surface on input.  On output, "constr" will be
  /// updated to include the newly constrained nodes as well. 
  void find_intersections( const COM::DataItem *pmesh, const COM::DataItem *disp,
			   COM::DataItem *pos,COM::DataItem *constr);

  /// Constrain the displacements, disp, to the appropriate positions, pos.
  void constrain_displacements( const COM::DataItem *pmesh, COM::DataItem *disp,
				const COM::DataItem *pos, const COM::DataItem *constr);

  /// Sets bflag = 0 for every element for which every node's cflag = 1
  void burnout( const COM::DataItem *pmesh, const COM::DataItem *cflag,
		COM::DataItem *bflag);
  
  /// Sets target = 0 for every burning pane element on which bflag = 0
  void burnout_filter( const COM::DataItem *bflag,
		       COM::DataItem *target);
  
  std::string &WindowName() { return (_wname); };

  static void load( const std::string &mname);

  static void unload( const std::string &mname);

  void SetVerbosity(int v){ _verbose = v; };
protected:
  int                     _verbose;       // Verbose level
  int                     _rank;       // process rank.
  MPI_Comm                _comm;
  double                  _TOL;
  std::string             _wname;
  MeshBndSurf            *_mbs;
  bool                    _checked;
  double                  _transform_x;
  double                  _transform_y;
  double                  _transform_z;
  bool                    _transform;
};

#endif // _ROCON_H_

