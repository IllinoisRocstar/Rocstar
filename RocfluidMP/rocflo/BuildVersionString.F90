! *********************************************************************
! * Rocstar Simulation Suite                                          *
! * Copyright@2015, Illinois Rocstar LLC. All rights reserved.        *
! *                                                                   *
! * Illinois Rocstar LLC                                              *
! * Champaign, IL                                                     *
! * www.illinoisrocstar.com                                           *
! * sales@illinoisrocstar.com                                         *
! *                                                                   *
! * License: See LICENSE file in top level of distribution package or *
! * http://opensource.org/licenses/NCSA                               *
! *********************************************************************
! *********************************************************************
! * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
! * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
! * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
! * NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
! * COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
! * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
! * Arising FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
! * USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
! *********************************************************************
!******************************************************************************
!
! Purpose: Build version string for printing in header.
!
! Description: none.
!
! Input: none.
!
! Output: 
!   versionString = string containing version number and date.
!
! Notes: 
!   1. The strings are NOT to be edited by anyone except the main code 
!      developer. 
!   2. Marks Rocbuild program will edit the build string to insert a 
!      build number or identifier
!
!******************************************************************************
!
! $Id: BuildVersionString.F90,v 1.50 2008/12/06 08:44:25 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE BuildVersionString( versionString )

  USE ModDataTypes
  IMPLICIT NONE

! ... parameters
  CHARACTER(*) :: versionString

! ... local variables
  CHARACTER(LEN=2)  :: major, minor, patch
  CHARACTER(LEN=3)  :: build
  CHARACTER(CHRLEN) :: date

!******************************************************************************
! set strings: DO NOT EDIT UNLESS YOU ARE MAIN DEVELOPER

  major = '4'
  minor = '2'
  patch = '0'
  build = '0' ! to be edited by Rocbuild

  date  = '05/11/05'

! write into string

  WRITE(versionString,'(A)') TRIM(major)//'.'//TRIM(minor)//'.'//TRIM(patch)
  WRITE(versionString,'(A)') 'Version: '//TRIM(versionString)//'-'//TRIM(build)
  WRITE(versionString,'(A)') TRIM(versionString)//', Date: '//TRIM(date)

END SUBROUTINE BuildVersionString

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BuildVersionString.F90,v $
! Revision 1.50  2008/12/06 08:44:25  mtcampbe
! Updated license.
!
! Revision 1.49  2008/11/19 22:17:36  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.48  2005/05/11 20:21:18  wasistho
! changed version
!
! Revision 1.47  2004/11/29 22:44:35  wasistho
! update v-s
!
! Revision 1.46  2003/11/20 16:40:36  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.42  2003/10/01 23:52:10  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.41  2003/08/28 20:05:39  jblazek
! Added acceleration terms.
!
! Revision 1.40  2003/08/25 21:51:24  jblazek
! Full version of global grid motion scheme.
!
! Revision 1.39  2003/08/11 21:51:18  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.38  2003/07/03 21:48:45  jblazek
! Implemented dual-time stepping.
!
! Revision 1.37  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.36  2003/05/29 17:28:43  jblazek
! Implemented Roe scheme.
!
! Revision 1.35  2003/05/20 20:46:57  jblazek
! Values in edge & corner cells now corrected at noslip and symmetry walls.
!
! Revision 1.34  2003/05/19 21:18:21  jblazek
! Automated switch to 0th-order extrapolation at slip walls and injection boundaries.
!
! Revision 1.33  2003/05/15 02:57:03  jblazek
! Inlined index function.
!
! Revision 1.32  2003/05/06 20:05:39  jblazek
! Corrected bug in grid motion (corner "averaging").
!
! Revision 1.31  2003/04/10 00:30:06  jblazek
! Corrected bug in grid movement algorithm.
!
! Revision 1.30  2003/04/04 21:05:00  jblazek
! Corrected bug in dumping out the solution.
!
! Revision 1.29  2003/03/14 22:05:11  jblazek
! Improved mesh motion algorithm - node movement exchaged between blocks.
!
! Revision 1.28  2003/03/04 21:56:47  jblazek
! Corrected bug for predictor-corrector iterations.
!
! Revision 1.27  2003/02/28 21:04:27  jblazek
! Corrected bug in send/recv of CE cells on single processor.
!
! Revision 1.26  2003/02/14 22:32:36  jblazek
! Finished implementation of corener and edge cells.
!
! Revision 1.25  2003/02/03 19:20:47  jblazek
! Added treatment of edge and corner cells for one processor.
!
! Revision 1.24  2003/01/23 17:48:53  jblazek
! Changed algorithm to dump convergence, solution and probe data.
!
! Revision 1.23  2003/01/10 17:58:43  jblazek
! Added missing explicit interfaces.
!
! Revision 1.22  2002/12/27 22:07:14  jblazek
! Splitted up RFLO_ModInterfaces and ModInterfaces.
!
! Revision 1.21  2002/12/06 22:29:26  jblazek
! Corrected bug for geometry exchange between minimal patches.
!
! Revision 1.20  2002/10/25 18:36:47  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.19  2002/10/23 18:43:24  jblazek
! Changed temporary pointer arrays into allocatable arrays
! in grid and solution I/O routines.
!
! Revision 1.18  2002/10/16 18:30:38  jblazek
! Within GenX, BC data at t=0 are updated in FlowSolver before calling
! the time-stepping routine.
!
! Revision 1.17  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.16  2002/10/04 19:33:18  jblazek
! Corrected one more bug in GenX restart.
!
! Revision 1.15  2002/10/03 21:25:39  jblazek
! Changed init. of burnig boundaries for GenX.
!
! Revision 1.14  2002/09/27 00:57:09  jblazek
! Changed makefiles - no makelinks needed.
!
! Revision 1.13  2002/09/20 22:22:36  jblazek
! Finalized integration into GenX.
!
! Revision 1.12  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.11  2002/08/30 19:08:58  jblazek
! Dimensions of work arrays now set in derivedInputValues.
!
! Revision 1.10  2002/08/30 01:47:58  jblazek
! Added support for moving grids.
!
! Revision 1.9  2002/08/16 21:35:59  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.8  2002/08/15 19:48:05  jblazek
! Implemented grid deformation capability.
!
! Revision 1.7  2002/08/07 21:07:55  jblazek
! Executables are left in the top-level directory instead of $(HOME)/bin.
!
! Revision 1.6  2002/07/25 01:00:08  jblazek
! Option for TVD type pressure switch.
!
! Revision 1.5  2002/07/22 22:59:11  jblazek
! Some more clean up.
!
! Revision 1.4  2002/07/16 21:34:37  jblazek
! Prefixed screen output with SOLVER_NAME.
!
! Revision 1.3  2002/07/05 23:20:46  jblazek
! Corrected bug in perfgasDependentVars.F90; did some cosmetics.
!
! Revision 1.2  2002/06/22 01:13:38  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/13 23:05:42  jblazek
! Added version string.
!
!******************************************************************************






