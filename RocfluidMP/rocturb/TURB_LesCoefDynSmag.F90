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
! Purpose: Obtain model coefficient for dynamic Smagorinsky model.
!
! Description: none. 
!
! Input: region  = data of current region
!        ibn,ien = begin and end node index
!        ijk     = ijk-face is being treated
!
! Output: cModel = turb%coef = model coefficient
!
! Notes: none.
!
!******************************************************************************
!
! $Id: TURB_LesCoefDynSmag.F90,v 1.7 2008/12/06 08:44:41 mtcampbe Exp $
!
! Copyright: (c) 2001 by the University of Illinois
!
!******************************************************************************

SUBROUTINE TURB_LesCoefDynSmag( region,ibn,ien,ijk )

  USE ModDataTypes
  USE ModDataStruct, ONLY      : t_region
  USE ModTurbulence, ONLY      : t_turb
  USE ModGlobal, ONLY          : t_global
#ifdef RFLU
  USE ModBndPatch, ONLY        : t_patch
  USE TURB_ModInterfaces, ONLY : TURB_FluLesBLij, TURB_FluLesBMij
  USE TURB_ModInterfaces, ONLY : TURB_FluLesC2F
#endif                                 
#ifdef RFLO
  USE TURB_ModInterfaces, ONLY : TURB_FloLesGenC2F
#endif                                 
  USE TURB_ModInterfaces, ONLY : TURB_LesLij, TURB_LesMij, TURB_LesContract
  USE ModError
  USE ModParameters
  USE TURB_ModParameters
  IMPLICIT NONE

! ... parameters
#ifdef RFLO
  TYPE(t_region), TARGET  :: region
#endif
#ifdef RFLU
  TYPE(t_region), POINTER :: region
#endif
  INTEGER                 :: ibn,ien,ijk

! ... loops variables
  INTEGER :: iPatch

! ... local variables
  CHARACTER(CHRLEN)       :: RCSIdentString
  TYPE(t_global), POINTER :: global
  TYPE(t_turb), POINTER   :: turb
#ifdef RFLU
  TYPE(t_patch), POINTER  :: patch
#endif

  INTEGER               :: tNDel(DIRI:DIRK), errorFlag
  REAL(RFREAL), POINTER :: lij(:,:)
#ifdef RFLO
  INTEGER               :: iLev, filterType
#endif
#ifdef RFLU
  INTEGER               :: nPatches, nBFaces
  REAL(RFREAL), POINTER :: bLij(:,:)
#endif

!******************************************************************************

  RCSIdentString = '$RCSfile: TURB_LesCoefDynSmag.F90,v $'

  global => region%global
  CALL RegisterFunction( global,'Turb_LesCoefDynSmag',&
  'TURB_LesCoefDynSmag.F90' )

! get some parameters and pointers ------------------------------------------

#ifdef RFLO
  iLev       =  region%currLevel
  filterType =  region%turbInput%filterType
  turb       => region%levels(iLev)%turb
#endif
#ifdef RFLU
  nPatches = region%grid%nPatches
  nBFaces  = 0

  DO iPatch = 1,nPatches
    patch   => region%patches(iPatch)
    nBFaces =  nBFaces + patch%nBTris + patch%nBQuads
  END DO ! iPatch
  turb => region%turb
#endif

! settle the filter-width for the test-filter 

  tNDel(DIRI) = 2*region%turbInput%filterWidth(DIRI)
#ifdef RFLO
  tNDel(DIRJ) = 2*region%turbInput%filterWidth(DIRJ)
  tNDel(DIRK) = 2*region%turbInput%filterWidth(DIRK)
#endif
#ifdef RFLU
  tNDel(DIRJ) = 0._RFREAL
  tNDel(DIRK) = 0._RFREAL
#endif

! first transfer rho and rho*ui from the centers to the faces since we 
! require the turbulent stress tensor in these points; store the results 
! in rhof, rhou1f, rhou2f and rhou3f; face-cv, fVar, is already allocated

#ifdef RFLO
  CALL TURB_FloLesGenC2F( region,ijk )
#endif
#ifdef RFLU
  CALL TURB_FluLesC2F( region )
#endif

! allocate and point arrays required by LesLij (ffVar, lij), LesMij (ffVar),
! and by lesContract (lij)
  
  ALLOCATE( turb%ffVar(CV_TURB_NELM,ibn:ien),stat=errorFlag )    
  ALLOCATE( turb%lij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )  
  ALLOCATE( turb%mij(TENSOR_SYMM_NELM,ibn:ien),stat=errorFlag )  
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  lij => turb%lij

#ifdef RFLU  
  ALLOCATE( turb%bffVar(CV_TURB_NELM,nBFaces),stat=errorFlag )    
  ALLOCATE( turb%bLij(TENSOR_SYMM_NELM,nBFaces),stat=errorFlag )  
  ALLOCATE( turb%bMij(TENSOR_SYMM_NELM,nBFaces),stat=errorFlag )  
  global%error = errorFlag
  IF (global%error /= 0) CALL ErrorStop( global,ERR_ALLOCATE,__LINE__ )

  bLij => turb%bLij
#endif

! determine components of lij for which we use the routine TURB_LesLij
! with test-filter instead of bar-filter

  CALL TURB_LesLij( region,ijk,tNDel,lij )
#ifdef RFLU
  CALL TURB_FluLesBLij( region,tNDel,bLij )
#endif

! determine components of mij which we store in mij using rhof and rhoftb 
! already available from TURB_LesLij

  CALL TURB_LesMij( region,ijk )
  DEALLOCATE( turb%ffVar )
#ifdef RFLU
  CALL TURB_FluLesBMij( region )
  DEALLOCATE( turb%bffVar )
#endif

! contract lij and mij to obtain model-coefficient at cell faces

  CALL TURB_LesContract( region,ijk )
  DEALLOCATE( turb%lij,turb%mij )
#ifdef RFLU
  DEALLOCATE( turb%bLij,turb%bMij )
#endif

! finalize --------------------------------------------------------------------

  CALL DeregisterFunction( global )

END SUBROUTINE TURB_LesCoefDynSmag

!******************************************************************************
!
! RCS Revision history:
! 
! $Log: TURB_LesCoefDynSmag.F90,v $
! Revision 1.7  2008/12/06 08:44:41  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:54  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2004/07/30 22:34:33  wasistho
! replaced floLesUniC2F by floLesGenC2F
!
! Revision 1.4  2004/05/28 02:00:45  wasistho
! update unstructured grid LES
!
! Revision 1.3  2004/03/19 02:49:17  wasistho
! prepared for RFLU
!
! Revision 1.2  2004/03/12 02:55:35  wasistho
! changed rocturb routine names
!
! Revision 1.1  2004/03/05 04:37:00  wasistho
! changed nomenclature
!
! Revision 1.3  2003/08/29 01:40:56  wasistho
! Added TARGET attribute to region variable, since pointers are cached into it
!
! Revision 1.2  2002/10/16 01:59:27  wasistho
! Changed global%error flag
!
! Revision 1.1  2002/10/14 23:55:29  wasistho
! Install Rocturb
!
!
!******************************************************************************







