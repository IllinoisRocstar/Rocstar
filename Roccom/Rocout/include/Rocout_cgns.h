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
/** \file Rocout_cgns.h
 *  Declaration of Rocout CGNS routines.
 */
/*  Author John Norris
 *  Initial date:   Oct. 19, 2004
 */
#ifndef _ROCOUT_CGNS_H
#define _ROCOUT_CGNS_H

#include "roccom.h"
#include <string>

/**
 ** Write the data for the given attribute to file.
 **
 ** Write the given attribute to file using the CGNS format.  The attribute
 ** may be a "mesh", "all" or some other predefined attribute.
 **
 ** \param fname The name of the main datafile. (Input)
 ** \param mname The name of the optional mesh datafile. (Input)
 ** \param attr The attribute to write out. (Input)
 ** \param material The name of the material. (Input)
 ** \param timelevel The simulation time for this data. (Input)
 ** \param pane_id The id for the local pane. (Input)
 ** \param ghosthandle "ignore" or "write" on ghost data.
 ** \param errorhandle "ignore", "warn", or "abort" on errors.
 ** \param mode Write == 0, append == 1. (Input)
 **/
void write_attr_CGNS(const std::string& fname, const std::string& mfile, 
                     const COM::Attribute* attr, const char* material, 
                     const char* timelevel, int pane_id,
                     const std::string& ghosthandle,
                     const std::string& errorhandle, int mode);

#endif // !defined(_ROCOUT_CGNS_H)






