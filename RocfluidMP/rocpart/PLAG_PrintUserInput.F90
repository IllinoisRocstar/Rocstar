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
! Purpose: print user input for PLAG.
!
! Description: none.
!
! Input: region = current region.
!
! Output: none.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_PrintUserInput.F90,v 1.5 2008/12/06 08:44:34 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_PrintUserInput( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModGlobal, ONLY          : t_global
  USE ModMixture, ONLY         : t_mixt_input
  USE ModPartLag, ONLY         : t_plag_input
  USE ModError
  USE ModParameters
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region) :: region

! ... loop variables
  INTEGER :: iCont

! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(CHRLEN) :: msg

  INTEGER :: iMat

  TYPE(t_global),     POINTER :: global

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_PrintUserInput.F90,v $ $Revision: 1.5 $'

  global => region%global

  CALL RegisterFunction( global, 'PLAG_PrintUserInput',&
  'PLAG_PrintUserInput.F90' )

! Print all PLAG input values -------------------------------------------------

  WRITE(STDOUT,1015) SOLVER_NAME//'        nPclsMax         ', &
                     region%plagInput%nPclsMax
  WRITE(STDOUT,1010) SOLVER_NAME//'        ejecModel        ', &
                     region%plagInput%ejecModel
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcVelRatio     ', &
                     region%plagInput%injcVelRatio
  WRITE(STDOUT,1020) SOLVER_NAME//'        spLoad           ', &
                     region%plagInput%spLoad
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcBeta         ', &
                     region%plagInput%injcBeta
  WRITE(STDOUT,1010) SOLVER_NAME//'        injcDiamDist     ', &
                     region%plagInput%injcDiamDist
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcDiamMean     ', &
                     region%plagInput%injcDiamMean
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcDiamMin      ', &
                     region%plagInput%injcDiamMin
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcDiamMax      ', &
                     region%plagInput%injcDiamMax
  WRITE(STDOUT,1020) SOLVER_NAME//'        injcStdDev       ', &
                     region%plagInput%injcStdDev
  WRITE(STDOUT,1010) SOLVER_NAME//'        intrplMixtModel  ', &
                     region%plagInput%intrplMixtModel
  WRITE(STDOUT,1010) SOLVER_NAME//'        breakupModel     ', &
                     region%plagInput%breakupModel
  WRITE(STDOUT,1020) SOLVER_NAME//'        breakupFac       ', &
                     region%plagInput%breakupFac
  WRITE(STDOUT,1010) SOLVER_NAME//'        breakupWebSwi     ', &
                     region%plagInput%breakupWebSwi
  WRITE(STDOUT,1015) SOLVER_NAME//'        nPclsBuffTot      ', &
                     region%plagInput%nPclsBuffTot
  WRITE(STDOUT,1015) SOLVER_NAME//'        nPclsBuffCECellsMax ', &
                     region%plagInput%nPclsBuffCECellsMax
  WRITE(STDOUT,1010) SOLVER_NAME//'        findPclModel      ', &
                     region%plagInput%findPclMethod                

! Print all PLAG input for Constituents values --------------------------------

  WRITE(STDOUT,1010) SOLVER_NAME//'        nCont            ', &
                     region%plagInput%nCont

  DO iCont = 1, region%plagInput%nCont
    iMat = region%plagInput%materialIndex(iCont)
    WRITE(STDOUT,2020) SOLVER_NAME//'          name(',iCont,')             =  '&
                       //TRIM(global%materials(iMat)%name)
    WRITE(STDOUT,2025) SOLVER_NAME//'          materialIndex(',iCont,')    =', &
                       region%plagInput%materialIndex(iCont)
    WRITE(STDOUT,2010) SOLVER_NAME//'          molw(',iCont,')             =', &
                       region%plagInput%molw(iCont)
    WRITE(STDOUT,2015) SOLVER_NAME//'          dens(',iCont,')             =', &
                       region%plagInput%dens(iCont)
    WRITE(STDOUT,2015) SOLVER_NAME//'          spht(',iCont,')             =', &
                       region%plagInput%spht(iCont)
    WRITE(STDOUT,2015) SOLVER_NAME//'          surfTens(',iCont,')         =', &
                       region%plagInput%surftens(iCont)
    WRITE(STDOUT,2015) SOLVER_NAME//'          injcMassFluxRatio(',iCont,')=', &
                       region%plagInput%injcMassFluxRatio(iCont)
    WRITE(STDOUT,2015) SOLVER_NAME//'          injcTemp(',iCont,')         =', &
                       region%plagInput%injcTemp(iCont)
  ENDDO ! iCont

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

1000 FORMAT(/,A,1X,80('-'))
1005 FORMAT(/,A)
1010 FORMAT(A,' = ',I2)
1015 FORMAT(A,' = ',I8)
1020 FORMAT(A,' = ',E12.5)

2010 FORMAT(A,I2,A,ES15.5)
2015 FORMAT(A,I2,A,EN15.5)
2020 FORMAT(A,I2,A)
2025 FORMAT(A,I2,A,I4)

END SUBROUTINE PLAG_PrintUserInput

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_PrintUserInput.F90,v $
! Revision 1.5  2008/12/06 08:44:34  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2007/03/06 23:15:32  fnajjar
! Renamed nPclsTot to nPclsMax
!
! Revision 1.2  2006/09/18 20:31:01  fnajjar
! Added injcTemp in printout
!
! Revision 1.1  2004/12/01 20:58:02  fnajjar
! Initial revision after changing case
!
! Revision 1.12  2004/10/08 22:13:02  haselbac
! Added printing of findPclMethod
!
! Revision 1.11  2004/06/17 15:19:03  fnajjar
! Added infrastructure for ejection model
!
! Revision 1.10  2004/06/16 23:07:17  fnajjar
! Renamed variabled for CRE kernel
!
! Revision 1.9  2004/03/10 23:11:36  fnajjar
! Included maximum buffer sizes to output file
!
! Revision 1.8  2004/03/05 22:09:03  jferry
! created global variables for peul, plag, and inrt use
!
! Revision 1.7  2003/09/17 21:05:13  fnajjar
! Added infrastructure for skewed Log distribution in injection model
!
! Revision 1.6  2003/09/15 20:24:57  fnajjar
! Fixed screen value for surface tension
!
! Revision 1.5  2003/09/15 15:50:38  fnajjar
! Added printing surface tension in outputfile
!
! Revision 1.4  2003/09/13 20:14:22  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.3  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.2  2003/03/12 21:21:46  fnajjar
! Use Material datastructure for Particle properties
!
! Revision 1.1  2002/10/25 14:18:35  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







