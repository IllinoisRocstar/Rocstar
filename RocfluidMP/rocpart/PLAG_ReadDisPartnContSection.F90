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
! Purpose: read in user input related to discrete particle module.
!
! Description: none.
!
! Input: user input file.
!
! Output: regions = number of constituents and their properties.
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_ReadDisPartnContSection.F90,v 1.5 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_ReadDisPartnContSection( regions )

  USE ModDataTypes 
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModPartLag, ONLY    : t_plag_input
  USE ModInterfaces, ONLY : ReadListSection, ReadPrefixedListSection
  USE ModInterfacesInteract, ONLY : INRT_SetMaterial
  USE ModError
  USE ModParameters
  USE ModMaterials, ONLY  : t_material
  USE PLAG_ModParameters
  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iVal, iReg
  
! ... local variables
  CHARACTER(CHRLEN) :: RCSIdentString
  CHARACTER(15)     :: keysnCont
  CHARACTER(CHRLEN), POINTER :: strValsnCont(:)

  INTEGER :: brbeg, brend, errorFlag, nCols, nCont, nRegions, nRows

  LOGICAL :: definednCont
    
  REAL(RFREAL), POINTER :: valsnCont(:,:)
  
  TYPE(t_global),   POINTER :: global 
  TYPE(t_material), POINTER :: pMaterial  

!******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_ReadDisPartnContSection.F90,v $ $Revision: 1.5 $'

  global => regions(1)%global
  
  CALL RegisterFunction( global, 'PLAG_ReadDisPartnContSection',&
  'PLAG_ReadDisPartnContSection.F90' )

! Read section pertinent to constituents --------------------------------------

  definednCont = .FALSE.
  keysnCont    = 'NCONT'
!  nCols        = 1
  nCols        = 2
  
  nRegions = global%nRegions
   
#ifdef RFLO    
  brbeg = 1
  brend = nRegions
#endif
#ifdef RFLU
  brbeg = LBOUND(regions,1)
  brend = UBOUND(regions,1)
#endif

!  CALL ReadListSection( global, IF_INPUT,keysnCont,nCols,nRows, &
!                        valsnCont,definednCont )

  CALL ReadPrefixedListSection( global, IF_INPUT,keysnCont,nCols,nRows, &
                                valsnCont,strValsnCont,definednCont )  

  IF (definednCont) THEN
    regions(brbeg:brend)%plagInput%nCont        = nRows

    DO iReg = brbeg,brend
      ALLOCATE( regions(iReg)%plagInput%molw(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'regions%plagInput%molw' ) 
      END IF ! global%error
    
      ALLOCATE( regions(iReg)%plagInput%dens(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'regions%plagInput%dens' ) 
      END IF ! global%error
    
      ALLOCATE( regions(iReg)%plagInput%spht(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'regions%plagInput%spht' ) 
      END IF ! global%error
    
      ALLOCATE( regions(iReg)%plagInput%surftens(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ ,'regions%plagInput%surftens' ) 
      END IF ! global%error
    
      ALLOCATE( regions(iReg)%plagInput%injcMassFluxRatio(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                        'regions%plagInput%injcMassFluxRatio' )
      END IF ! global%error
    
      ALLOCATE( regions(iReg)%plagInput%materialIndex(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                        'regions%plagInput%imaterialIndex' )
      END IF ! global%error

      ALLOCATE( regions(iReg)%plagInput%materialName(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                       'regions%plagInput%materialName' ) 
      END IF ! global%error
      
      ALLOCATE( regions(iReg)%plagInput%injcTemp(nRows),stat=errorFlag )
      global%error = errorFlag
      IF (global%error /= ERR_NONE) THEN
        CALL ErrorStop( global, ERR_ALLOCATE,__LINE__ , &
                        'regions%plagInput%injcTemp' )
      END IF ! global%error
    
    END DO !iReg
    
    DO iReg= brbeg,brend
      nCont = regions(iReg)%plagInput%nCont   
      DO iVal=1, nCont  
        CALL INRT_SetMaterial(global,pMaterial,strValsnCont(iVal))
        
        regions(iReg)%plagInput%materialIndex(iVal) = pMaterial%index
        regions(iReg)%plagInput%materialName(iVal)  = pMaterial%name

        regions(iReg)%plagInput%molw(iVal) = pMaterial%molw
        regions(iReg)%plagInput%dens(iVal) = pMaterial%dens
        regions(iReg)%plagInput%spht(iVal) = pMaterial%spht
        regions(iReg)%plagInput%surftens(iVal) = pMaterial%surftens
        regions(iReg)%plagInput%injcMassFluxRatio(iVal) = ABS(valsnCont(iVal,1))
        regions(iReg)%plagInput%injcTemp(iVal)          = ABS(valsnCont(iVal,2))
      ENDDO ! ival
    ENDDO ! iReg
    
  ELSE
    CALL ErrorStop( global, ERR_MISSING_VALUE,__LINE__) 
  ENDIF ! definednCont  

! deallocate pointer arrays

  DEALLOCATE( strValsnCont, stat=errorFlag ) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'strValsnCont' ) 
  END IF ! global%error
  
  DEALLOCATE( valsnCont, stat=errorFlag ) 
  global%error = errorFlag
  IF (global%error /= ERR_NONE) THEN
    CALL ErrorStop( global, ERR_DEALLOCATE,__LINE__ ,'valsnCont' ) 
  END IF ! global%error
   
! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_ReadDisPartnContSection

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ReadDisPartnContSection.F90,v $
! Revision 1.5  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/09/18 20:37:35  fnajjar
! Added injcTemp in particle input
!
! Revision 1.2  2005/12/01 18:39:15  fnajjar
! Added error trap for missing input value
!
! Revision 1.1  2004/12/01 20:58:04  fnajjar
! Initial revision after changing case
!
! Revision 1.7  2004/02/26 21:02:20  haselbac
! Added RFLU support
!
! Revision 1.6  2003/09/26 21:48:20  fnajjar
! Changed interface call for INRT_SetMaterial to ModInterfacesInteract
!
! Revision 1.5  2003/09/13 20:14:22  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.4  2003/03/24 23:28:41  jferry
! converted SetMaterial to INRT_SetMaterial
!
! Revision 1.3  2003/03/12 21:21:46  fnajjar
! Use Material datastructure for Particle properties
!
! Revision 1.2  2003/01/10 17:03:41  f-najjar
! Bug fix to read properly input data for all regions
!
! Revision 1.1  2002/10/25 14:19:16  f-najjar
! Initial Import of Rocpart
!
!
!******************************************************************************







