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
! Purpose: Specify deformation of patches so that moving-grid calculations 
!   can be carried out outside of GENX.
!
! Description: None.
!
! Input: 
!   region     Region data
!
! Output: None.
!
! Notes: 
!   1. This routine will have to be hard-coded for each case.
!   2. Only specify deformation for interior vertices.
!
! ******************************************************************************
!
! $Id: RFLU_USER_GetDeformation.F90,v 1.18 2008/12/06 08:44:30 mtcampbe Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RFLU_USER_GetDeformation(region)

  USE ModDataTypes
  USE ModDataStruct, ONLY: t_region
  USE ModBndPatch, ONLY: t_patch
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModMPI
  USE ModParameters

  IMPLICIT NONE

! *****************************************************************************
! Declarations and definitions
! *****************************************************************************

! =============================================================================
! Arguments
! =============================================================================

  TYPE(t_region) :: region

! =============================================================================
! Locals
! =============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: ibv,iPatch,ivg,speedSign
  REAL(RFREAL) :: nx,ny,nz,ps,ps1,ps2,x,x1,x2,y,z
  TYPE(t_global), POINTER :: global
  TYPE(t_patch), POINTER :: pPatch

! *****************************************************************************
! Start
! *****************************************************************************

  RCSIdentString = '$RCSfile: RFLU_USER_GetDeformation.F90,v $ $Revision: 1.18 $'

  global => region%global

  CALL RegisterFunction(global,'RFLU_USER_GetDeformation',&
  'RFLU_USER_GetDeformation.F90')

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME,'Getting deformation...'
  END IF ! global%myProcid

! *****************************************************************************
! Select case-dependent boundary patch deformation
! *****************************************************************************

  SELECT CASE ( TRIM(global%casename) )   

! =============================================================================
!   Simple box
! =============================================================================    
    
    CASE ( "box_hex2" )
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)
        
        IF ( pPatch%iPatchGlobal == 6 ) THEN 
          DO ibv = 1,pPatch%nBVert                  
            pPatch%dXyz(XCOORD,ibv) = 8.0_RFREAL*global%dtMin 
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL*global%dtMin             
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch    
      END DO ! iPatch      

! =============================================================================
!   Burning crack problem
! =============================================================================    
    
    CASE ( "burncrack" )
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)
        
        IF ( pPatch%iPatchGlobal == 1 ) THEN 
          DO ibv = 1,pPatch%nBVert
            nx = pPatch%bvn(XCOORD,ibv)
            ny = pPatch%bvn(YCOORD,ibv)
            nz = pPatch%bvn(ZCOORD,ibv)
    
            ivg = pPatch%bv(ibv)
            
            x = region%grid%xyz(XCOORD,ivg)
            y = region%grid%xyz(YCOORD,ivg)
            z = region%grid%xyz(ZCOORD,ivg)
               
            x1  =  0.04_RFREAL
            x2  =  0.07_RFREAL               
            ps1 = 10.00_RFREAL   
            ps2 =  5.00_RFREAL 
                          
            ps = (ps2 - ps1)/(x2 - x1)*x - (x1*ps2 - x2*ps1)/(x2 - x1) 
    
            pPatch%dXyz(XCOORD,ibv) = ps*global%dtMin*nx
            pPatch%dXyz(YCOORD,ibv) = ps*global%dtMin*ny
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL 
          END DO ! ibv  
        ELSE IF ( pPatch%iPatchGlobal == 4 ) THEN 
          DO ibv = 1,pPatch%nBVert
            nx = pPatch%bvn(XCOORD,ibv)
            ny = pPatch%bvn(YCOORD,ibv)
            nz = pPatch%bvn(ZCOORD,ibv)
    
            pPatch%dXyz(XCOORD,ibv) = 5.0_RFREAL*global%dtMin
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL 
          END DO ! ibv            
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch    
      END DO ! iPatch   
        
! =============================================================================
!   Deforming cube (check whether volume is conserved)
! =============================================================================

    CASE ( "cube_def" )
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)           

        IF ( pPatch%iPatchGlobal == 2 .OR. & 
             pPatch%iPatchGlobal == 5 .OR. & 
             pPatch%iPatchGlobal == 6 ) THEN 
          DO ibv = 1,pPatch%nBVert
            ivg = pPatch%bv(ibv)

            z = region%grid%xyz(ZCOORD,ivg) 
                                                  
            ps = 100.0_RFREAL*z/0.1_RFREAL

            pPatch%dXyz(XCOORD,ibv) = ps*global%dtMin
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL
          END DO ! ibv            
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv        
        END IF ! iPatch       
      END DO ! iPatch        
        
! =============================================================================
!   Cube with tetrahedra only - NOTE different patch numbering...
! =============================================================================    

    CASE ( "cube2pt","cube3pt","cube11pt","cube21pt" ) 
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)

        IF ( pPatch%iPatchGlobal == 5 ) THEN 
          DO ibv = 1,pPatch%nBVert                  
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 8.0_RFREAL*global%dtMin             
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch

    CASE ( "cube6pt" )    
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)

        IF ( pPatch%iPatchGlobal == 4 ) THEN 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = -100.0_RFREAL*global%dtMin                  
            pPatch%dXyz(YCOORD,ibv) =    0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) =    0.0_RFREAL       
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch        
        
! =============================================================================
!   Endburner problem 
! =============================================================================    

    CASE ( "endburner3pt","endburner5pt","endburner9pt" )      
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)                   

        IF ( pPatch%iPatchGlobal == 5 ) THEN                 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL 
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 8.0_RFREAL*global%dtMin            
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch
      
! =============================================================================
!   Endburner problem (new)
! =============================================================================    

    CASE ( "endburner3ptnew","endburner5ptnew","endburner9ptnew" )      
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)                   

        IF ( pPatch%iPatchGlobal == 6 ) THEN         
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL 
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL          
            pPatch%dXyz(ZCOORD,ibv) = 8.0_RFREAL*global%dtMin            
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch
      
! =============================================================================
!   Endburner problem (angled)
! =============================================================================    

    CASE ( "endburner3pt_angled" )      
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)                   

        IF ( pPatch%iPatchGlobal == 5 ) THEN                 
          DO ibv = 1,pPatch%nBVert
            nx = pPatch%bvn(XCOORD,ibv)
            ny = pPatch%bvn(YCOORD,ibv)
            nz = pPatch%bvn(ZCOORD,ibv)
          
            pPatch%dXyz(XCOORD,ibv) = 8.0_RFREAL*global%dtMin*nx 
            pPatch%dXyz(YCOORD,ibv) = 8.0_RFREAL*global%dtMin*ny  
            pPatch%dXyz(ZCOORD,ibv) = 8.0_RFREAL*global%dtMin*nz             
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch      
      
! =============================================================================
!   Piston problem
! =============================================================================    

    CASE ( "piston_exp","piston_comp" )
      IF ( TRIM(global%casename) == "piston_exp" ) THEN 
        speedSign = -1 
      ELSE 
        speedSign =  1
      END IF ! TRIM(global%casename)
          
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)                   

        IF ( pPatch%iPatchGlobal == 3 ) THEN                 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 10.0_RFREAL*global%dtMin*speedSign
            pPatch%dXyz(YCOORD,ibv) =  0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) =  0.0_RFREAL            
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch 
      END DO ! iPatch                  
      
! =============================================================================
!   Scalability problem
! =============================================================================    
    
    CASE ( "scalability" )
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)
        
        IF ( pPatch%iPatchGlobal == 1 ) THEN 
          DO ibv = 1,pPatch%nBVert
            nx = pPatch%bvn(XCOORD,ibv)
            ny = pPatch%bvn(YCOORD,ibv)
            nz = pPatch%bvn(ZCOORD,ibv)
    
            pPatch%dXyz(XCOORD,ibv) = 1.0_RFREAL*global%dtMin*nx
            pPatch%dXyz(YCOORD,ibv) = 1.0_RFREAL*global%dtMin*ny
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL 
          END DO ! ibv  
        ELSE 
          DO ibv = 1,pPatch%nBVert
            pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL      
          END DO ! ibv
        END IF ! iPatch    
      END DO ! iPatch   

! =============================================================================
!   Stargrain slice
! =============================================================================

    CASE ( "starslice" )
      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch)               

        IF ( pPatch%iPatchGlobal == 1 ) THEN 
          DO ibv = 1,pPatch%nBVert
            nx = pPatch%bvn(XCOORD,ibv)
            ny = pPatch%bvn(YCOORD,ibv)
            nz = pPatch%bvn(ZCOORD,ibv)

            pPatch%dXyz(XCOORD,ibv) = 10.0_RFREAL*global%dtMin*nx
            pPatch%dXyz(YCOORD,ibv) =  0.0_RFREAL
            pPatch%dXyz(ZCOORD,ibv) = 10.0_RFREAL*global%dtMin*nz                         
          END DO ! ibv
        END IF ! iPatch
      END DO ! iPatch

! =============================================================================
!   Default
! =============================================================================        
           
    CASE DEFAULT 
      global%warnCounter = global%warnCounter + 1     
    
      IF ( global%myProcid == MASTERPROC .AND. & 
           global%verbLevel > VERBOSE_LOW ) THEN        
        WRITE(STDOUT,'(A,3(1X,A))') SOLVER_NAME,'*** WARNING ***', & 
                                    'No displacements specified.', & 
                                    'Setting displacements to zero.'
      END IF ! global

      DO iPatch=1,region%grid%nPatches
        pPatch => region%patches(iPatch) 

        DO ibv = 1,pPatch%nBVert
          pPatch%dXyz(XCOORD,ibv) = 0.0_RFREAL
          pPatch%dXyz(YCOORD,ibv) = 0.0_RFREAL
          pPatch%dXyz(ZCOORD,ibv) = 0.0_RFREAL           
        END DO ! ibv
      END DO ! iPatch
 
    END SELECT ! global%casename
    
! *****************************************************************************
! Print diagnostic information on boundary patch movement
! *****************************************************************************
   
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW .AND. & 
       region%grid%nPatches > 0 ) THEN     
    WRITE(STDOUT,'(A,3X,A)') SOLVER_NAME,'Deformation extrema:'       
    WRITE(STDOUT,'(A,5X,A)') SOLVER_NAME,'Written only for patches '// &
                             'with non-zero actual vertices.'
     
    DO iPatch=1,region%grid%nPatches
      pPatch => region%patches(iPatch)       

      IF ( pPatch%nBVert > 0 ) THEN 
        WRITE(STDOUT,'(A,5X,A,1X,I3)') SOLVER_NAME,'Patch:',iPatch 

        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.x:', & 
              MINVAL(pPatch%dXyz(XCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(XCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.y:', & 
              MINVAL(pPatch%dXyz(YCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(YCOORD,1:pPatch%nBVert))
        WRITE(STDOUT,'(A,7X,A,2(1X,E15.8))') SOLVER_NAME,'dXyz.z:', & 
              MINVAL(pPatch%dXyz(ZCOORD,1:pPatch%nBVert)), & 
              MAXVAL(pPatch%dXyz(ZCOORD,1:pPatch%nBVert))
      END IF ! pPatch%nBVert
    END DO ! iPatch                                                                           
  END IF ! global%myProcid

! *****************************************************************************
! End
! *****************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
    WRITE(STDOUT,'(A,1X,A,1X,I3)') SOLVER_NAME,'Getting deformation done.'
  END IF ! global%myProcid
 
  CALL DeregisterFunction(global)


END SUBROUTINE RFLU_USER_GetDeformation

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_USER_GetDeformation.F90,v $
! Revision 1.18  2008/12/06 08:44:30  mtcampbe
! Updated license.
!
! Revision 1.17  2008/11/19 22:17:43  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.16  2005/05/26 22:06:19  haselbac
! Cosmetics only
!
! Revision 1.15  2004/03/15 21:04:51  haselbac
! Added deformations for piston problem
!
! Revision 1.14  2003/08/20 02:09:58  haselbac
! Changed verbosity conditions to reduce solver output in GENx runs
!
! Revision 1.13  2003/07/22 02:07:36  haselbac
! Added global%warnCounter
!
! Revision 1.12  2003/04/12 16:39:42  haselbac
! Added burning crack, reordered cases
!
! Revision 1.11  2003/03/31 16:18:47  haselbac
! Added CASE for endburner3pt_angled
!
! Revision 1.10  2003/03/25 19:17:48  haselbac
! Added deformation specification for case box_hex2
!
! Revision 1.9  2003/03/18 21:35:10  haselbac
! Added new endburner deformation, got lost in last merge...
!
! Revision 1.8  2003/03/15 19:03:40  haselbac
! Use iPatchGlobal, changes for || runs
!
! Revision 1.7  2003/02/20 19:48:50  haselbac
! Rewrote with different loop-select order, added starslice and default
!
! Revision 1.6  2003/02/01 00:26:35  haselbac
! Added deformation for new endburner geometry
!
! Revision 1.5  2003/01/28 14:53:43  haselbac
! Clean-up and added cube_def case
!
! Revision 1.4  2003/01/03 22:06:24  haselbac
! Added CASE for 6pt cube - different patch numbering
!
! Revision 1.3  2002/11/15 14:09:56  haselbac
! Added endburner patch deformation
!
! Revision 1.2  2002/11/08 21:36:56  haselbac
! Made deformation case-specific
!
! Revision 1.1  2002/10/27 19:20:28  haselbac
! Initial revision
!
! ******************************************************************************







