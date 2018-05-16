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
! Purpose: Suite of routines to carry out interpolation operations.
!
! Description: None.
!
! Notes: None.
!  1. Removed RFLU_InterpCells2FacesPatches because it uses bf2bg
!     and bf2bg is being removed as bGradFace is moved from pRegion to patch
!
! ******************************************************************************
!
! $Id: RFLU_ModInterpolation.F90,v 1.21 2008/12/06 08:44:22 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModInterpolation

  USE ModDataTypes
  USE ModError
  USE ModParameters  
  USE ModGlobal, ONLY: t_global
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModTools, ONLY: CompFact
  USE ModMPI

  USE RFLU_ModStencilsUtils, ONLY: RFLU_ComputeStencilWeights

  USE ModInterfaces, ONLY: RFLU_InvertMatrixSVD

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: RFLU_InterpCells2Face, &
            RFLU_InterpCells2FacePatch, &  
            RFLU_InterpCells2Faces, &
            RFLU_InterpCells2Verts, & 
            RFLU_InterpSimpleCells2Verts
  
  SAVE
  
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
       
  CHARACTER(CHRLEN) :: RCSIdentString = & 
    '$RCSfile: RFLU_ModInterpolation.F90,v $ $Revision: 1.21 $'
       
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  
  
  
! ******************************************************************************
!
! Purpose: Interpolation from cells to single face.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   ifg                 Global face index
!   src                 Source array of cell data
!
! Output:
!   dst                 Destination array of face data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InterpCells2Face(pRegion,ifg,src,dst)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: ifg
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: src
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: dst
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: errorFlag,icg,iDst,iEndDst,iEndSrc,iSrc,isl
    REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz, &
                    r11,r12,r13,r14,r22,r23,r24,r33,r34,r44,term,wti
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_InterpCells2Face',&
  'RFLU_ModInterpolation.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid        

    iEndSrc = UBOUND(src,1)
    iEndDst = UBOUND(dst,1)

! ******************************************************************************
!   Get face weights 
! ******************************************************************************

! ==============================================================================  
!   Select appropriate dimensionality
! ==============================================================================  

    SELECT CASE ( pRegion%mixtInput%dimens )
    
! ------------------------------------------------------------------------------
!     Two dimensions 
! ------------------------------------------------------------------------------    
    
      CASE ( 2 )
        r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
        r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
        r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
        r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
        r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
        r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)

        c11 = 1.0_RFREAL/r11
        c22 = 1.0_RFREAL/r22
        c33 = 1.0_RFREAL/r33

        c12 = - c11*r12
        c13 = -(c11*r13 + c12*c22*r23) 

        c23 = - c22*r23
        
! ------------------------------------------------------------------------------
!     Three dimensions 
! ------------------------------------------------------------------------------    
            
      CASE ( 3 )            
        r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
        r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
        r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
        r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
        r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
        r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)
        r14 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_14)
        r24 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_24)
        r34 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_34)
        r44 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_44)                  

        c11 = 1.0_RFREAL/r11
        c22 = 1.0_RFREAL/r22
        c33 = 1.0_RFREAL/r33
        c44 = 1.0_RFREAL/r44

        c12 = - c11*r12
        c13 = -(c11*r13 + c12*c22*r23) 
        c14 = -(c11*r14 + c12*c22*r24 + c13*c33*r34) 

        c23 = - c22*r23
        c24 = -(c22*r24 + c23*c33*r34)

        c34 = - c33*r34
        
! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------    
            
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************
!   Initialize destination array
! ******************************************************************************

    DO iDst = 1,iEndDst ! Explicit loop because of Frost problems
      dst(iDst) = 0.0_RFREAL
    END DO ! iDst

! ******************************************************************************
!   Loop over stencil members and interpolate
! ******************************************************************************

! ==============================================================================  
!   Select appropriate dimensionality
! ==============================================================================  

    SELECT CASE ( pRegion%mixtInput%dimens ) 

! ------------------------------------------------------------------------------
!     Two dimensions 
! ------------------------------------------------------------------------------    

      CASE ( 2 )
        DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
          icg = pGrid%f2cs(ifg)%cellMembs(isl)

          dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
          dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)

          term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

          dx = term*dx
          dy = term*dy        

          wti = term*c33*c33*(term + c23*dy + c13*dx)

          iDst = 1

          DO iSrc = 1,iEndSrc
            dst(iDst) = dst(iDst) + wti*src(iSrc,icg)

            iDst = iDst + 1          
          END DO ! iSrc                    
        END DO ! isl       

! ------------------------------------------------------------------------------
!     Three dimensions 
! ------------------------------------------------------------------------------    

      CASE ( 3 ) 
        DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
          icg = pGrid%f2cs(ifg)%cellMembs(isl)

          dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
          dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)
          dz = pGrid%cofg(ZCOORD,icg) - pGrid%fc(ZCOORD,ifg)   

          term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

          dx = term*dx
          dy = term*dy
          dz = term*dz        

          wti = term*c44*c44*(term + c34*dz + c24*dy + c14*dx)

          iDst = 1

          DO iSrc = 1,iEndSrc
            dst(iDst) = dst(iDst) + wti*src(iSrc,icg)

            iDst = iDst + 1          
          END DO ! iSrc                    
        END DO ! isl

! ------------------------------------------------------------------------------
!     Default
! ------------------------------------------------------------------------------    

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens 
 
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InterpCells2Face  
  
  
  
  

! ******************************************************************************
!
! Purpose: Interpolation from cells to single face on single patch.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   pPatch              Pointer to patch data
!   ifg                 Global face index
!   src                 Source array of cell data
!
! Output:
!   dst                 Destination array of face data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InterpCells2FacePatch(pRegion,pPatch,ifl,src,dst)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: ifl
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: src
    REAL(RFREAL), DIMENSION(:), INTENT(OUT) :: dst
    TYPE(t_patch), POINTER :: pPatch
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: errorFlag,icg,iDst,iEndDst,iEndSrc,iSrc,isl
    REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz, &
                    r11,r12,r13,r14,r22,r23,r24,r33,r34,r44,term,wti
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_InterpCells2FacePatch',&
  'RFLU_ModInterpolation.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid        

    iEndSrc = UBOUND(src,1)
    iEndDst = UBOUND(dst,1)

! ******************************************************************************
!   Get face weights 
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%dimens )
      CASE ( 2 )
        r11 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11)           
        r12 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_12) 
        r22 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_22)           
        r13 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_13) 
        r23 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_23)                   
        r33 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_33)

        c11 = 1.0_RFREAL/r11
        c22 = 1.0_RFREAL/r22
        c33 = 1.0_RFREAL/r33

        c12 = - c11*r12
        c13 = -(c11*r13 + c12*c22*r23)  

        c23 = - c22*r23
      CASE ( 3 ) 
        r11 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_11)           
        r12 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_12) 
        r22 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_22)           
        r13 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_13) 
        r23 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_23)                   
        r33 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_33)
        r14 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_14)
        r24 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_24)
        r34 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_34)
        r44 = pPatch%bf2cs(ifl)%xyzMoms(XYZ_MOM_44)                  

        c11 = 1.0_RFREAL/r11
        c22 = 1.0_RFREAL/r22
        c33 = 1.0_RFREAL/r33
        c44 = 1.0_RFREAL/r44

        c12 = - c11*r12
        c13 = -(c11*r13 + c12*c22*r23) 
        c14 = -(c11*r14 + c12*c22*r24 + c13*c33*r34) 

        c23 = - c22*r23
        c24 = -(c22*r24 + c23*c33*r34)

        c34 = - c33*r34
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens

! ******************************************************************************
!   Initialize destination array
! ******************************************************************************

    DO iDst = 1,iEndDst ! Explicit loop because of Frost problems
      dst(iDst) = 0.0_RFREAL
    END DO ! iDst

! ******************************************************************************
!   Loop over stencil members and interpolate
! ******************************************************************************

    SELECT CASE ( pRegion%mixtInput%dimens ) 
      CASE ( 2 ) 
        DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
          icg = pPatch%bf2cs(ifl)%cellMembs(isl)

          dx = pGrid%cofg(XCOORD,icg) - pPatch%fc(XCOORD,ifl)
          dy = pGrid%cofg(YCOORD,icg) - pPatch%fc(YCOORD,ifl)

          term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

          dx = term*dx
          dy = term*dy       

          wti = term*c33*c33*(term + c23*dy + c13*dx)

          iDst = 1

          DO iSrc = 1,iEndSrc
            dst(iDst) = dst(iDst) + wti*src(iSrc,icg)

            iDst = iDst + 1          
          END DO ! iSrc                    
        END DO ! isl      
      CASE ( 3 ) 
        DO isl = 1,pPatch%bf2cs(ifl)%nCellMembs 
          icg = pPatch%bf2cs(ifl)%cellMembs(isl)

          dx = pGrid%cofg(XCOORD,icg) - pPatch%fc(XCOORD,ifl)
          dy = pGrid%cofg(YCOORD,icg) - pPatch%fc(YCOORD,ifl)
          dz = pGrid%cofg(ZCOORD,icg) - pPatch%fc(ZCOORD,ifl)   

          term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

          dx = term*dx
          dy = term*dy
          dz = term*dz        

          wti = term*c44*c44*(term + c34*dz + c24*dy + c14*dx)

          iDst = 1

          DO iSrc = 1,iEndSrc
            dst(iDst) = dst(iDst) + wti*src(iSrc,icg)

            iDst = iDst + 1          
          END DO ! iSrc                    
        END DO ! isl
      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens
 
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InterpCells2FacePatch
  
  
  
  
  
  
  
! ******************************************************************************
!
! Purpose: Interpolation from cells to faces.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   src                 Source array of cell data
!
! Output:
!   dst                 Destination array of face data
!
! Notes: 
!   1. Could call RFLU_InterpCells2Face for each face, but do not do so for
!      performance reasons (repeated call with src array) 
!
! ******************************************************************************

  SUBROUTINE RFLU_InterpCells2Faces(pRegion,src,dst)
 
    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: src
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: dst
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: errorFlag,icg,iDst,iEndDst,iEndSrc,ifg,isl,iSrc
    REAL(RFREAL) :: c11,c12,c13,c14,c22,c23,c24,c33,c34,c44,dx,dy,dz, &
                    r11,r12,r13,r14,r22,r23,r24,r33,r34,r44,term,wti
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_InterpCells2Faces',&
  'RFLU_ModInterpolation.F90')

! ******************************************************************************
!   Set pointers and variables
! ******************************************************************************

    pGrid => pRegion%grid

    iEndSrc = UBOUND(src,1)
    iEndDst = UBOUND(dst,1)

! ******************************************************************************
!   Loop over faces and interpolate from cells in stencil
! ******************************************************************************
 
! ==============================================================================
!   Select appropriate dimensionality
! ============================================================================== 
 
    SELECT CASE ( pRegion%mixtInput%dimens ) 
    
! --- Two dimensions -----------------------------------------------------------    
    
      CASE ( 2 ) 
        DO ifg = 1,pGrid%nFaces                                                       
          r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
          r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
          r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
          r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
          r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
          r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)

          c11 = 1.0_RFREAL/r11
          c22 = 1.0_RFREAL/r22
          c33 = 1.0_RFREAL/r33

          c12 = - c11*r12
          c13 = -(c11*r13 + c12*c22*r23) 

          c23 = - c22*r23
          
          DO iDst = 1,iEndDst ! Explicit loop because of Frost problems
            dst(iDst,ifg) = 0.0_RFREAL
          END DO ! iDst

          DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
            icg = pGrid%f2cs(ifg)%cellMembs(isl)

            dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
            dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy)            

            dx = term*dx
            dy = term*dy       

            wti = term*c33*c33*(term + c23*dy + c13*dx)

            iDst = 1

            DO iSrc = 1,iEndSrc
              dst(iDst,ifg) = dst(iDst,ifg) + wti*src(iSrc,icg)

              iDst = iDst + 1          
            END DO ! iSrc                    
          END DO ! isl                  
        END DO ! ifg      

! --- Three dimensions ---------------------------------------------------------    

      CASE ( 3 )  
        DO ifg = 1,pGrid%nFaces                                                       
          r11 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_11)           
          r12 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_12) 
          r22 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_22)           
          r13 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_13) 
          r23 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_23)                   
          r33 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_33)
          r14 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_14)
          r24 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_24)
          r34 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_34)
          r44 = pGrid%f2cs(ifg)%xyzMoms(XYZ_MOM_44)                  

          c11 = 1.0_RFREAL/r11
          c22 = 1.0_RFREAL/r22
          c33 = 1.0_RFREAL/r33
          c44 = 1.0_RFREAL/r44

          c12 = - c11*r12
          c13 = -(c11*r13 + c12*c22*r23) 
          c14 = -(c11*r14 + c12*c22*r24 + c13*c33*r34) 

          c23 = - c22*r23
          c24 = -(c22*r24 + c23*c33*r34)

          c34 = - c33*r34

          DO iDst = 1,iEndDst ! Explicit loop because of Frost problems
            dst(iDst,ifg) = 0.0_RFREAL
          END DO ! iDst

          DO isl = 1,pGrid%f2cs(ifg)%nCellMembs 
            icg = pGrid%f2cs(ifg)%cellMembs(isl)

            dx = pGrid%cofg(XCOORD,icg) - pGrid%fc(XCOORD,ifg)
            dy = pGrid%cofg(YCOORD,icg) - pGrid%fc(YCOORD,ifg)
            dz = pGrid%cofg(ZCOORD,icg) - pGrid%fc(ZCOORD,ifg)   

            term = 1.0_RFREAL/SQRT(dx*dx + dy*dy + dz*dz)            

            dx = term*dx
            dy = term*dy
            dz = term*dz        

            wti = term*c44*c44*(term + c34*dz + c24*dy + c14*dx)

            iDst = 1

            DO iSrc = 1,iEndSrc
              dst(iDst,ifg) = dst(iDst,ifg) + wti*src(iSrc,icg)

              iDst = iDst + 1          
            END DO ! iSrc                    
          END DO ! isl                  
        END DO ! ifg

! --- Default -------------------------------------------------------------------    

      CASE DEFAULT
        CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END SELECT ! pRegion%mixtInput%dimens
 
! ******************************************************************************
!   End
! ******************************************************************************

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InterpCells2Faces







! ******************************************************************************
!
! Purpose: Interpolation from cells to vertices.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   orderNominal        Nominal order of accuracy
!   nVar                Number of variables
!   src                 Source data array of cell data
!
! Output:
!   dst                 Destination data array of vertex data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InterpCells2Verts(pRegion,orderNominal,nVar,src,dst)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: nVar,orderNominal
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: src
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: dst
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: errorFlag,icg,isl,ivg,iVar,nRows,orderActual,sCount        
    REAL(RFREAL), DIMENSION(:,:), ALLOCATABLE :: dr,wts
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_InterpCells2Verts',&
  'RFLU_ModInterpolation.F90')

    IF ( global%myProcid == MASTERPROC ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                                 'Interpolating from cells to vertices...'
      END IF ! global%verbLevel
      
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Method: Proper' 
      END IF ! global%verbLevel > VERBOSE_LOW       
    END IF ! global%verbLevel 

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over vertices and interpolate from cells in stencil
! ******************************************************************************
 
    DO ivg = 1,pGrid%nVertTot                                                       
      nRows = pGrid%v2cs(ivg)%nCellMembs      

      orderActual = orderNominal      

      ALLOCATE(dr(XCOORD:ZCOORD,nRows),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'dr')
      END IF ! global%error 

      ALLOCATE(wts(1,nRows),STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'wts')
      END IF ! global%error

      DO isl = 1,pGrid%v2cs(ivg)%nCellMembs   
        icg = pGrid%v2cs(ivg)%cellMembs(isl)

        dr(XCOORD,isl) = pGrid%cofg(XCOORD,icg) - pGrid%xyz(XCOORD,ivg)
        dr(YCOORD,isl) = pGrid%cofg(YCOORD,icg) - pGrid%xyz(YCOORD,ivg)
        dr(ZCOORD,isl) = pGrid%cofg(ZCOORD,icg) - pGrid%xyz(ZCOORD,ivg)                 
      END DO ! isl

! ==============================================================================  
!     Compute interpolation weights
! ==============================================================================  

      CALL RFLU_ComputeStencilWeights(global,pRegion%mixtInput%dimens, & 
                                      COMPWTS_MODE_ADAPT,COMPWTS_SCAL_INVDIST, &
                                      DERIV_DEGREE_0,orderActual,nRows,dr, &
                                      wts,sCount)

! ==============================================================================  
!     Initialize destination array
! ==============================================================================  

      DO iVar = 1,nVar ! Explicit loop because of ASCI White problems      
        dst(iVar,ivg) = 0.0_RFREAL
      END DO ! iVar

! ==============================================================================  
!     Interpolate   
! ==============================================================================  
                                                    
      DO isl = 1,pGrid%v2cs(ivg)%nCellMembs
        icg = pGrid%v2cs(ivg)%cellMembs(isl)

        DO iVar = 1,nVar                                        
          dst(iVar,ivg) = dst(iVar,ivg) + wts(1,isl)*src(iVar,icg)                                    
        END DO ! iVar                             
      END DO ! isl         
                
! ==============================================================================  
!     Deallocate memory         
! ==============================================================================  
        
      DEALLOCATE(dr,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'dr')
      END IF ! global%error 

      DEALLOCATE(wts,STAT=errorFlag)
      global%error = errorFlag   
      IF ( global%error /= ERR_NONE ) THEN 
        CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'wts')
      END IF ! global%error                                          
    END DO ! ivg
 
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Interpolating from cells to vertices done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InterpCells2Verts
    
    
    
    



! ******************************************************************************
!
! Purpose: Interpolation from cells to vertices by simple averaging.
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   nVar                Number of variables
!   src                 Source data array of cell data
!
! Output:
!   dst                 Destination data array of vertex data
!
! Notes: None.
!
! ******************************************************************************

  SUBROUTINE RFLU_InterpSimpleCells2Verts(pRegion,nVar,src,dst)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================  
!   Arguments 
! ==============================================================================  

    INTEGER, INTENT(IN) :: nVar
    REAL(RFREAL), DIMENSION(:,:), INTENT(IN) :: src
    REAL(RFREAL), DIMENSION(:,:), INTENT(OUT) :: dst
    TYPE(t_region), POINTER :: pRegion

! ==============================================================================  
!   Locals 
! ==============================================================================  

    INTEGER :: icg,icl,ivg,iVar,nCells,v2cBeg,v2cEnd 
    REAL(RFREAL) :: wt       
    TYPE(t_grid), POINTER :: pGrid
    TYPE(t_global), POINTER :: global
    
! ******************************************************************************
!   Start
! ******************************************************************************

    global => pRegion%global

    CALL RegisterFunction(global,'RFLU_InterpSimpleCells2Verts',&
  'RFLU_ModInterpolation.F90')

    IF ( global%myProcid == MASTERPROC ) THEN 
      IF ( global%verbLevel > VERBOSE_NONE ) THEN
        WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                                 'Interpolating from cells to vertices...'
      END IF ! global%verbLevel
      
      IF ( global%verbLevel > VERBOSE_LOW ) THEN
        WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Method: Simple' 
      END IF ! global%verbLevel > VERBOSE_LOW       
    END IF ! global%verbLevel  

! ******************************************************************************
!   Set pointers
! ******************************************************************************

    pGrid => pRegion%grid
    
! ******************************************************************************
!   Loop over vertices and interpolate from cells in stencil
! ******************************************************************************
 
    DO ivg = 1,pGrid%nVertTot 
      v2cBeg = pGrid%v2cInfo(V2C_BEG,ivg)
      v2cEnd = pGrid%v2cInfo(V2C_END,ivg)
                                                              
      nCells = v2cEnd - v2cBeg + 1
          
! ==============================================================================  
!     Initialize destination array
! ==============================================================================  

      DO iVar = 1,nVar ! Explicit loop because of ASCI White problems      
        dst(iVar,ivg) = 0.0_RFREAL
      END DO ! iVar

! ==============================================================================  
!     Interpolate   
! ==============================================================================  
     
      wt = 1.0_RFREAL/REAL(nCells,RFREAL)     
                                                    
      DO icl = 1,nCells
        icg = pGrid%v2c(v2cBeg+icl-1)

        DO iVar = 1,nVar                                        
          dst(iVar,ivg) = dst(iVar,ivg) + wt*src(iVar,icg)                                    
        END DO ! iVar                             
      END DO ! isl                                                            
    END DO ! ivg
 
! ******************************************************************************
!   End
! ******************************************************************************

    IF ( global%myProcid == MASTERPROC .AND. &
         global%verbLevel > VERBOSE_NONE ) THEN
      WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, & 
                               'Interpolating from cells to vertices done.'            
    END IF ! global%verbLevel  

    CALL DeregisterFunction(global)

  END SUBROUTINE RFLU_InterpSimpleCells2Verts

    
    
    


! ******************************************************************************
! End
! ******************************************************************************

END MODULE RFLU_ModInterpolation


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModInterpolation.F90,v $
! Revision 1.21  2008/12/06 08:44:22  mtcampbe
! Updated license.
!
! Revision 1.20  2008/11/19 22:17:33  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.19  2006/08/19 15:39:23  mparmar
! Removed RFLU_InterpCells2FacesPatches routine
!
! Revision 1.18  2006/04/07 15:19:19  haselbac
! Removed tabs
!
! Revision 1.17  2005/10/05 13:57:47  haselbac
! Adapted to changes in module contents
!
! Revision 1.16  2005/03/09 14:56:49  haselbac
! Added code for 2d interpolation
!
! Revision 1.15  2004/07/21 15:00:06  haselbac
! Added RFLU_InterpSimpleCells2Verts, cosmetics
!                                                    
! Revision 1.14  2004/03/18 03:31:56  haselbac                                           
! Added routines for interp from cells to faces, clean-up                                
!
! Revision 1.13  2004/01/22 16:03:59  haselbac                                           
! Made contents of modules PRIVATE, only procs PUBLIC, to avoid errors on ALC 
! and titan  
!
! Revision 1.12  2003/12/04 03:28:52  haselbac                                           
! Complete rewrite                                                                       
!
! Revision 1.11  2003/07/22 15:39:14  haselbac                                           
! Added Nullify routines, distinction betw PUBLIC and PRIVATE members                    
!
! Revision 1.10  2003/07/22 02:05:14  haselbac                                           
! Added comp of proper wghts for cell-to-vertex interp                                   
!
! Revision 1.9  2003/03/15 18:12:59  haselbac                                            
! Now also interpolate to dummy vertices                                                 
!
! Revision 1.8  2003/01/28 16:34:29  haselbac                                            
! Cosmetics only                                                                         
!
! Revision 1.7  2002/10/09 20:47:53  haselbac                                            
! Fixed bug in RFLU_DestroyInterpolant: Missing errorFlag declaration                    
!
! Revision 1.6  2002/10/08 15:49:21  haselbac                                            
! {IO}STAT=global%error replaced by {IO}STAT=errorFlag - SGI problem                     
!
! Revision 1.5  2002/09/09 15:09:41  haselbac                                            
! global now under regions                                                               
!
! Revision 1.4  2002/07/25 15:02:59  haselbac                                            
! Only write out for MASTERPROC                                                          
!
! Revision 1.3  2002/06/17 13:39:45  haselbac                                            
! Prefixed SOLVER_NAME to all screen output                                              
!
! Revision 1.2  2002/05/04 16:39:51  haselbac                                            
! Added PRIVATE dataitem to RCSIdentString                                              
!
! Revision 1.1  2002/04/11 18:48:48  haselbac                                            
! Initial revision                                                                       
!
! ******************************************************************************











