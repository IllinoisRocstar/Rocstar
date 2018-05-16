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
! *****************************************************************************
!
! Purpose: Read in user input related to physical materials.
!
! Description: None.
!
! Input: 
!   global 	Pointer to global data
!
! Output: None.
!
! Notes: None.
!
! *****************************************************************************
!
! $Id: INRT_ReadMaterialInput.F90,v 1.5 2008/12/06 08:44:32 mtcampbe Exp $
!
! Copyright: (c) 2003-2005 by the University of Illinois
!
! *****************************************************************************

SUBROUTINE INRT_ReadMaterialInput(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModMaterials
  USE ModError
  USE ModParameters

  USE ModBuildFileNames, ONLY: BuildFileNamePlain
    
  USE ModInterfaces, ONLY: ReadBothSection

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER, PARAMETER :: NSTRKEYS_MAX = 5
  INTEGER, PARAMETER :: NKEYS_MAX = 20
  LOGICAL :: strDefined(NSTRKEYS_MAX),defined(NKEYS_MAX)
  CHARACTER(20) :: strKeys(NSTRKEYS_MAX),keys(NKEYS_MAX)
  CHARACTER(256) :: line  
  CHARACTER(CHRLEN) :: RCSIdentString,strVals(NSTRKEYS_MAX)
  CHARACTER(CHRLEN) :: fname
  INTEGER :: errorFlag,nMat,iKeyMolw,iKeyDens,iKeySpht,iKeySurfTens, &
             iKeyTboil,iKeyTmelt,iMat,iPass,iStrKeyName,iStrKeyPhase, &
             nStrKeys,nKeys
  REAL(RFREAL) :: vals(NKEYS_MAX)
  TYPE(t_material), POINTER :: material

! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: INRT_ReadMaterialInput.F90,v $ $Revision: 1.5 $'

  CALL RegisterFunction( global,'INRT_ReadMaterialInput',&
  'INRT_ReadMaterialInput.F90' )

! ******************************************************************************
! Define keys for string-valued quantities
! ******************************************************************************

  iStrKeyName  = 1
  iStrKeyPhase = 2
  nStrKeys     = 2

  IF ( nStrKeys > NSTRKEYS_MAX ) THEN
    CALL ErrorStop(global,ERR_EXCEEDS_DECL_MEM,__LINE__)
  END IF ! nStrKeys

  strKeys(iStrKeyName)  = 'NAME'
  strKeys(iStrKeyPhase) = 'PHASE'

! ******************************************************************************
! Define keys for real-valued quantities
! ******************************************************************************

  iKeyMolw     = 1
  iKeyDens     = 2
  iKeySpht     = 3
  iKeySurfTens = 4
  iKeyTboil    = 5
  iKeyTmelt    = 6
  nKeys        = 6

  IF ( nKeys > NKEYS_MAX ) THEN 
    CALL ErrorStop(global,ERR_EXCEEDS_DECL_MEM,__LINE__)
  END IF ! nKeys

  keys(iKeyMolw)     = 'MOLW'
  keys(iKeyDens)     = 'DENS'
  keys(iKeySpht)     = 'SPHT'
  keys(iKeySurfTens) = 'SURFTENS'
  keys(iKeyTboil)    = 'TBOIL'
  keys(iKeyTmelt)    = 'TMELT'

! ******************************************************************************
! Search for MATERIAL sections
! ******************************************************************************

  CALL BuildFileNamePlain(global,FILEDEST_INDIR,'.inp',fname)

  DO iPass = 1,2

! ==============================================================================
!   Open file
! ==============================================================================

    OPEN(IF_INPUT,FILE=TRIM(fname),FORM='FORMATTED',STATUS='OLD', &
         IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error

! ==============================================================================
!   Read file looking for keywords
! ==============================================================================

    SELECT CASE ( iPass )

! -----------------------------------------------------------------------------
!     First pass: Count number of material sections and allocate materials
! -----------------------------------------------------------------------------

      CASE ( 1 ) 
        nMat = 0 ! initialize count of materials

        DO
          READ(IF_INPUT,'(A256)',ERR=10,END=86) line
          IF ( TRIM(line) == '# MATERIAL' ) THEN 
            nMat = nMat + 1
          END IF ! TRIM
        END DO ! <empty>

86      CONTINUE

        IF ( nMat > 0 ) THEN
          ALLOCATE(global%materials(nMat),STAT=errorFlag)
          global%error = errorFlag
          IF ( global%error /= ERR_NONE ) THEN
            CALL ErrorStop(global,ERR_ALLOCATE,__LINE__)
          END IF ! global%error
        ELSE 
          NULLIFY(global%materials)
        END IF ! nMat

        global%nMaterials = nMat

! -----------------------------------------------------------------------------
!     Second pass: Search for real- and string-valued keys
! -----------------------------------------------------------------------------

      CASE ( 2 ) 
        iMat = 0 ! initialize index of materials

        DO
          READ(IF_INPUT,'(A256)',ERR=10,END=87) line
          
          IF ( TRIM(line) == '# MATERIAL' ) THEN
            iMat = iMat + 1
            material => global%materials(iMat)

! --------- Set values to indicate the status of undefined --------------------

            material%molw     = -1.0_RFREAL
            material%dens     = -1.0_RFREAL
            material%spht     = -1.0_RFREAL
            material%surftens = -1.0_RFREAL
            material%Tboil    = -1.0_RFREAL
            material%Tmelt    = -1.0_RFREAL
            material%Phase    = -1

! --------- Read real- and string-valued keys ---------------------------------

            CALL ReadBothSection(global,IF_INPUT,nKeys,nStrKeys,keys, &
                                 strKeys,vals,strVals,defined,strDefined)

! --------- ensure that material is named -------------------------------------

            IF ( strDefined(iStrKeyName) .EQV. .TRUE. ) THEN
              material%name = TRIM(strVals(iStrKeyName))
            ELSE
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__)
            END IF ! strDefined(iStrKeyName)

! --------- Set material phase ------------------------------------------------

            IF ( strDefined(iStrKeyPhase) .EQV. .TRUE. ) THEN                                              
              SELECT CASE ( strVals(iStrKeyPhase)(1:1) )
                CASE ( 'G','g' )
                  material%phase = 1
                CASE ( 'L','l' )
                  material%phase = 2
                CASE ( 'S','s' )
                  material%phase = 3
                CASE DEFAULT
                  CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
                END SELECT ! strVals
            END IF ! strDefined

! --------- Set material index to iMat (its index in global%materials(:)) -----

            material%index = iMat

! --------- Set real-valued quantities ----------------------------------------

            IF ( defined(iKeyMolw) .EQV. .TRUE. ) THEN 
              material%molw = vals(iKeyMolw)
            ELSE 
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'MOLW')
            END IF ! defined(iKeyMolw)

            IF ( defined(iKeyDens) .EQV. .TRUE. ) THEN 
              material%dens = vals(iKeyDens)
            ELSE 
              IF ( material%phase /= 1 ) THEN 
                CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'DENS')
              END IF ! material%phase
            END IF ! defined(iKeyDens)

            IF ( defined(iKeySpht) .EQV. .TRUE. ) THEN
              material%spht = vals(iKeySpht)
            ELSE
              CALL ErrorStop(global,ERR_VAL_UNDEFINED,__LINE__,'SPHT')
            END IF ! defined(iKeySpht)

            IF ( defined(iKeySurfTens) .EQV. .TRUE. ) THEN 
              material%surftens = vals(iKeySurfTens)
            END IF ! defined(iKeySurfTens)

            IF ( defined(iKeyTboil) .EQV. .TRUE. ) THEN
              material%Tboil = vals(iKeyTboil)
            END IF ! defined(iKeyTboil)

            IF ( defined(iKeyTmelt) .EQV. .TRUE. ) THEN
              material%Tmelt = vals(iKeyTmelt)
            END IF ! defined(iKeyTmelt)
          END IF ! line
        END DO ! <empty>

87      CONTINUE

! -----------------------------------------------------------------------------
!     Default
! -----------------------------------------------------------------------------

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! iPass

! ==============================================================================
!   Close file
! ==============================================================================

    CLOSE(IF_INPUT,IOSTAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /=  ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,'File: '//TRIM(fname))
    END IF ! global%error
  END DO ! iPass

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)
  
  GOTO 999

10   CONTINUE
  CALL ErrorStop(global,ERR_FILE_READ,__LINE__,'File: '//TRIM(fname))

999  CONTINUE

END SUBROUTINE INRT_ReadMaterialInput

! *****************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_ReadMaterialInput.F90,v $
! Revision 1.5  2008/12/06 08:44:32  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:44  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2005/11/14 17:01:44  haselbac
! Added check for DENS, clean-up
!
! Revision 1.2  2005/11/10 02:32:35  haselbac
! Clean-up
!
! Revision 1.1  2004/12/01 21:56:39  fnajjar
! Initial revision after changing case
!
! Revision 1.6  2004/07/27 21:28:32  jferry
! minor bug fix
!
! Revision 1.5  2004/07/23 22:43:17  jferry
! Integrated rocspecies into rocinteract
!
! Revision 1.4  2004/04/15 16:04:21  jferry
! minor formatting (removed trailing spaces)
!
! Revision 1.3  2004/03/02 21:49:46  jferry
! Added melting and boiling point to material definitions
!
! Revision 1.2  2003/09/13 20:17:31  fnajjar
! Added surface tension to Materials datastructure
!
! Revision 1.1  2003/03/24 23:23:25  jferry
! converted from libfloflu routine to rocinteract routine
!
! *****************************************************************************







