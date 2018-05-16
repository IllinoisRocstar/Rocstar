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
! Purpose: Correct the mixture properties using higher-order
!          interpolation schemed onto particle locations
!
! Description: none.
!
! Input: pRegion = current region.
!
! Output: region%levels%plag%dv 
!
! Notes: none.
!
!******************************************************************************
!
! $Id: PLAG_RFLU_CorrectMixtProperties.F90,v 1.4 2008/12/06 08:44:35 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
!******************************************************************************

SUBROUTINE PLAG_RFLU_CorrectMixtProperties(pRegion)

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModGrid, ONLY: t_grid
  USE ModPartLag, ONLY    : t_plag, t_plag_input
  USE ModMixture, ONLY    : t_mixt
  USE ModError  
  USE ModParameters
  USE PLAG_ModParameters

  USE ModInterfaces, ONLY: MixtPerf_R_M, & 
                           MixtPerf_T_DPR

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================  
! Arguments 
! ==============================================================================  

  TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
! Locals
! ==============================================================================  

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: icg,indMol,intrplMixtModel,iPcl,nPcls
  INTEGER, POINTER, DIMENSION(:,:)    :: pAiv
  
  REAL(RFREAL) :: dx,dy,dz,gc,mm,p,r
  REAL(RFREAL), POINTER, DIMENSION(:,:)    :: pCv,pDv,pDvMixt,pGvMixt
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: pGrad
  
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag 
  TYPE(t_mixt),   POINTER :: pMixt 
 
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = '$RCSfile: PLAG_RFLU_CorrectMixtProperties.F90,v $ $Revision: 1.4 $'

  global => pRegion%global
  
  CALL RegisterFunction( global, 'PLAG_RFLU_CorrectMixtProperties',&
  'PLAG_RFLU_CorrectMixtProperties.F90' )

! ******************************************************************************
! Set variables and pointers
! ******************************************************************************

  pMixt => pRegion%mixt 
  pGrid => pRegion%grid
  pPlag => pRegion%plag

  pDvMixt => pMixt%dv
  pGvMixt => pMixt%gv

  nPcls  = pRegion%plag%nPcls
  indMol = pRegion%mixtInput%indMol

  intrplMixtModel = pRegion%plagInput%intrplMixtModel

! ******************************************************************************
! Set pointers for mixture 
! ******************************************************************************

  pDvMixt => pMixt%dv
  pGvMixt => pMixt%gv

  pGrad => pRegion%mixt%gradCell

! ******************************************************************************
! Set pointers for discrete particles 
! ******************************************************************************

  pAiv    => pPlag%aiv
  pCv     => pPlag%cv  
  pDv     => pPlag%dv

! ******************************************************************************
! Check that have primitive state vector
! ******************************************************************************

  IF ( pMixt%cvState /= CV_MIXT_STATE_DUVWP ) THEN 
    CALL ErrorStop(global,ERR_CV_STATE_INVALID,__LINE__)
  END IF ! pMixt%cvState

! ******************************************************************************
! Trap error for inconsistent input deck
! ******************************************************************************

  IF ( intrplMixtModel /= ZEROTH_ORDER .AND. &
       pRegion%mixtInput%spaceOrder == 1     ) THEN 
    CALL ErrorStop(global,ERR_PLAG_INTRPLMODEL,__LINE__)
  END IF ! intrplMixtModel

! ******************************************************************************        
! Correct discrete particle dv
! ******************************************************************************
    
  SELECT CASE ( intrplMixtModel )    
    CASE ( ZEROTH_ORDER )
    CASE ( FIRST_ORDER  )        
      DO iPcl = 1, nPcls
        icg = pAiv(AIV_PLAG_ICELLS,iPcl)

        dx = pCv(CV_PLAG_XPOS,iPcl) - pGrid%cofg(XCOORD,icg)
        dy = pCv(CV_PLAG_YPOS,iPcl) - pGrid%cofg(YCOORD,icg)
        dz = pCv(CV_PLAG_ZPOS,iPcl) - pGrid%cofg(ZCOORD,icg)                

        pDv(DV_PLAG_UVELMIXT,iPcl) = pDv(DV_PLAG_UVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_XVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_XVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_XVEL,icg)*dz
        pDv(DV_PLAG_VVELMIXT,iPcl) = pDv(DV_PLAG_VVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_YVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_YVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_YVEL,icg)*dz
        pDv(DV_PLAG_WVELMIXT,iPcl) = pDv(DV_PLAG_WVELMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_ZVEL,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_ZVEL,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_ZVEL,icg)*dz
        pDv(DV_PLAG_PRESMIXT,iPcl) = pDv(DV_PLAG_PRESMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_PRES,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_PRES,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_PRES,icg)*dz
        pDv(DV_PLAG_DENSMIXT,iPcl) = pDv(DV_PLAG_DENSMIXT,iPcl) & 
                                   + pGrad(XCOORD,GRC_MIXT_DENS,icg)*dx &
                                   + pGrad(YCOORD,GRC_MIXT_DENS,icg)*dy &
                                   + pGrad(ZCOORD,GRC_MIXT_DENS,icg)*dz 
      END DO ! iPcl
      
      SELECT CASE ( pRegion%mixtInput%gasModel ) 
        CASE ( GAS_MODEL_TCPERF,GAS_MODEL_MIXT_TCPERF )
          DO iPcl = 1,nPcls      
            icg = pAiv(AIV_PLAG_ICELLS,iPcl)
          
            mm = pGvMixt(GV_MIXT_MOL,indMol*icg)
            gc = MixtPerf_R_M(mm)
                  
            r = pDv(DV_PLAG_DENSMIXT,iPcl)
            p = pDv(DV_PLAG_PRESMIXT,iPcl)        
                  
            pDv(DV_PLAG_TEMPMIXT,iPcl) = MixtPerf_T_DPR(r,p,gc)
          END DO ! iPcl
        CASE DEFAULT
          CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
      END SELECT ! pRegion%mixtInput%gasModel
      
    CASE ( SECOND_ORDER )
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )
    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )     
  END SELECT  ! intrplMixtModel

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction( global )

END SUBROUTINE PLAG_RFLU_CorrectMixtProperties

!******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_CorrectMixtProperties.F90,v $
! Revision 1.4  2008/12/06 08:44:35  mtcampbe
! Updated license.
!
! Revision 1.3  2008/11/19 22:17:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.2  2006/04/07 15:19:24  haselbac
! Removed tabs
!
! Revision 1.1  2005/11/30 22:26:54  fnajjar
! Initial version
!
!******************************************************************************







