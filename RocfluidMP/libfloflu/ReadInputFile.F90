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
! ******************************************************************************
!
! Purpose: Read in user input (done on all processors).
!
! Description: None.
!
! Input: 
!   regions     Region data
!
! Output: 
!   None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadInputFile.F90,v 1.7 2008/12/06 08:44:09 mtcampbe Exp $
!
! Copyright: (c) 2001-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadInputFile(regions)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModBuildFileNames, ONLY: BuildFileNamePlain  
  
  USE ModInterfaces, ONLY: ReadFormatsSection, &
                           ReadReferenceSection, &
                           ReadFlowSection, & 
                           ReadProbeSection, & 
                           ReadForcesSection, &
                           ReadTimestepSection, &
                           ReadMixtureSection, & 
                           ReadMultigridSection, &
                           ReadNumericsSection, & 
                           ReadTransformSection, &
                           ReadViscositySection, &
                           ReadThrustSection, &
                           ReadInitFlowSection, &
                           ReadGridMotionSection, &
                           ReadAccelerationSection, &
                           ReadRandomSection, &
                           ReadPostSection
#ifdef INRT
  USE ModInterfacesInteract, ONLY: INRT_ReadMaterialInput
#endif
#ifdef RFLU
  USE ModInterfaces, ONLY: ReadPrepSection, &
                           ReadTimeZoomingSection, &
                           ReadRocketSection
#endif
#ifdef STATS
  USE ModInterfacesStatistics, ONLY: ReadStatisticSection
#endif
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: regions(:)

! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: fname
  CHARACTER(256) :: line
  INTEGER :: errorFlag
  TYPE(t_global), POINTER :: global

! ******************************************************************************
! Start
! ******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction(global,'ReadInputFile',&
  'ReadInputFile.F90')

! ******************************************************************************
! Open file
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.inp',fname)
  
  OPEN(IF_INPUT,FILE=TRIM(fname),FORM='FORMATTED',STATUS='OLD',IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! ******************************************************************************
! Read file looking for keywords
! ******************************************************************************

  DO
    READ(IF_INPUT,'(A256)',ERR=10,END=86) line

    SELECT CASE( TRIM(line) )
      CASE ('# FORMATS')
        CALL ReadFormatsSection(global)
      CASE ('# REFERENCE')
        CALL ReadReferenceSection(global)
      CASE ('# RANDOM')
        CALL ReadRandomSection(global)
      CASE ('# FLOWMODEL')
        CALL ReadFlowSection(regions)
#ifdef RFLU
      CASE ('# MIXTURE')
        CALL ReadMixtureSection(regions)
#endif
      CASE ('# VISCMODEL')
        CALL ReadViscositySection(regions)
      CASE ('# ACCELERATION')
        CALL ReadAccelerationSection(global)
      CASE ('# PROBE')
        CALL ReadProbeSection(global)
      CASE ('# FORCES')
        CALL ReadForcesSection(global)
      CASE ('# TIMESTEP')
        CALL ReadTimestepSection(global)
      CASE ('# MULTIGRID')
        CALL ReadMultigridSection(global)
      CASE ('# NUMERICS')
        CALL ReadNumericsSection(regions)
      CASE ('# THRUST')
        CALL ReadThrustSection(global)
      CASE ('# INITFLOW')
        CALL ReadInitFlowSection(regions) 
      CASE ('# POST') 
        CALL ReadPostSection(global)
#ifdef RFLO
      CASE ('# GRIDMOTION')
        CALL ReadGridMotionSection(global)
#endif
#ifdef RFLU
      CASE ('# TRANSFORM')
        CALL ReadTransformSection(global)      
      CASE ('# GRIDMOTION')
        CALL ReadGridMotionSection(regions)
      CASE ('# PREP')
        CALL ReadPrepSection(global)
      CASE ('# ROCKET')
        CALL ReadRocketSection(global)
      CASE ('# TIMEZOOMING')
        CALL ReadTimeZoomingSection(global)
#endif
#ifdef STATS
      CASE ('# STATISTICS')
        CALL ReadStatisticSection(global)
#endif
    END SELECT ! TRIM(line)
  END DO ! <empty>

86   CONTINUE

! ******************************************************************************
! Close file
! ******************************************************************************

  CLOSE(IF_INPUT,IOSTAT=errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN
    CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
  END IF ! global%error

! ******************************************************************************
! Read material sections
! ******************************************************************************

#ifdef INRT
  CALL INRT_ReadMaterialInput(global)
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE

END SUBROUTINE ReadInputFile

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadInputFile.F90,v $
! Revision 1.7  2008/12/06 08:44:09  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:23  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/04/14 14:24:12  mtcampbe
! Updated for TZ and Rocket case constraints
!
! Revision 1.4  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.3  2005/11/10 01:57:00  haselbac
! Added reading of MIXTURE section for Rocflu, clean-up
!
! Revision 1.2  2005/10/31 19:24:37  haselbac
! Added call to ReadMixtureSection, clean-up
!
! Revision 1.1  2004/12/01 16:50:28  haselbac
! Initial revision after changing case
!
! Revision 1.36  2004/11/29 17:13:57  wasistho
! use ModInterfacesStatistics
!
! Revision 1.35  2004/07/24 03:45:28  wasistho
! release readPostSection to flo and flu
!
! Revision 1.34  2004/04/08 03:15:13  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.33  2003/11/21 22:33:10  fnajjar
! Updated Random Number Generator
!
! Revision 1.32  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.29  2003/09/26 21:45:34  fnajjar
! Modified ModInterfaces calls to proper physical modules
!
! Revision 1.28  2003/08/28 20:05:38  jblazek
! Added acceleration terms.
!
! Revision 1.27  2003/08/11 21:51:17  jblazek
! Added basic global grid smoothing scheme.
!
! Revision 1.26  2003/06/02 17:12:00  jblazek
! Added computation of thrust.
!
! Revision 1.25  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.24  2003/04/28 22:39:12  haselbac
! Added readPostSection and readPrepSection for RFLU
!
! Revision 1.23  2003/04/10 23:25:53  fnajjar
! Added readViscositySection
!
! Revision 1.22  2003/03/24 23:25:48  jferry
! moved some routines from libfloflu to rocinteract
!
! Revision 1.21  2003/03/17 19:37:02  jblazek
! Added inDir to path to the input file.
!
! Revision 1.20  2003/03/11 16:04:19  jferry
! Enclosed USE statements for multi-physics routines within ifdefs
!
! Revision 1.19  2003/02/26 23:38:30  jferry
! eliminated end=999 convention to ensure that files get closed
!
! Revision 1.18  2003/01/28 16:16:32  haselbac
! Read transform section only for RFLU
!
! Revision 1.17  2003/01/28 16:08:31  haselbac
! Added new calls transform, initial solution (RFLU), and grid motion (RFLU)
!
! Revision 1.16  2002/10/12 03:20:50  jblazek
! Replaced [io]stat=global%error with local errorFlag for Rocflo.
!
! Revision 1.15  2002/10/08 15:48:35  haselbac
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem
!
! Revision 1.14  2002/09/20 22:22:35  jblazek
! Finalized integration into GenX.
!
! Revision 1.13  2002/09/09 14:02:21  haselbac
! mixtInput now under regions
!
! Revision 1.12  2002/09/05 17:40:20  jblazek
! Variable global moved into regions().
!
! Revision 1.10  2002/08/18 02:27:47  wasistho
! Added ReadTurbulenceSection
!
! Revision 1.9  2002/06/14 21:17:01  wasistho
! Added time avg statistics
!
! Revision 1.8  2002/06/07 16:40:37  jblazek
! Grid & solution for all regions in one file.
!
! Revision 1.7  2002/05/04 16:20:55  haselbac
! Cosmetic changes
!
! Revision 1.6  2002/03/26 19:03:37  haselbac
! Added ROCFLU functionality
!
! Revision 1.5  2002/02/21 23:25:05  jblazek
! Blocks renamed as regions.
!
! Revision 1.4  2002/01/28 23:55:22  jblazek
! Added flux computation (central scheme).
!
! Revision 1.3  2002/01/11 17:18:31  jblazek
! Updated description of I/O variables.
!
! Revision 1.2  2002/01/02 16:00:03  jblazek
! Added input for multigrid parameters.
!
! Revision 1.1  2001/12/07 16:54:31  jblazek
! Added files to read user input.
!
! ******************************************************************************







