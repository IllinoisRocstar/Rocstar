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
! Purpose: Suite of routines related to bilinear patches.
!
! Description: None.
!
! Notes: None. 
!
! ******************************************************************************
!
! $Id: RFLU_ModBilinearPatch.F90,v 1.5 2008/12/06 08:44:20 mtcampbe Exp $
!
! Copyright: (c) 2005 by the University of Illinois
!
! ******************************************************************************

MODULE RFLU_ModBilinearPatch

  USE ModGlobal, ONLY: t_global 
  USE ModParameters
  USE ModDataTypes
  USE ModError
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_region
  USE ModGrid, ONLY: t_grid
  USE ModMPI

  IMPLICIT NONE
   
  PRIVATE
  PUBLIC :: RFLU_BLIN_ComputeNormal, &
            RFLU_BLIN_ComputeXSectLine, &  
            RFLU_BLIN_FindClosestPoint
      
! ******************************************************************************
! Declarations and definitions
! ******************************************************************************  
     
  CHARACTER(CHRLEN), PRIVATE :: & 
    RCSIdentString = '$RCSfile: RFLU_ModBilinearPatch.F90,v $ $Revision: 1.5 $'        
       
! ******************************************************************************
! Routines
! ******************************************************************************

  CONTAINS
  





! ******************************************************************************
!
! Purpose: Compute normal vector of bilinear patch at given (u,v) coordinates.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   ax,ay,az    Components of bilinear patch parametric representation
!   bx,by,bz    Components of bilinear patch parametric representation
!   cx,cy,cz    Components of bilinear patch parametric representation
!   dx,dy,dz    Components of bilinear patch parametric representation
!   u           u-coordinate
!   v           v-coordinate
!
! Output: 
!   nx          x-component of normal vector
!   ny          y-component of normal vector
!   nz          z-component of normal vector
!
! Notes: None.
!
! ******************************************************************************

SUBROUTINE RFLU_BLIN_ComputeNormal(global,ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz, &
                                   u,v,nx,ny,nz)

  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Parameters      
! ==============================================================================

  REAL(RFREAL), INTENT(IN) :: ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,u,v
  REAL(RFREAL), INTENT(OUT) :: nx,ny,nz
  TYPE(t_global), POINTER :: global   

! ==============================================================================
! Locals      
! ==============================================================================

  REAL(RFREAL) :: imag

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'RFLU_BLIN_ComputeNormal',&
  'RFLU_ModBilinearPatch.F90')

! ******************************************************************************
! Compute normal
! ******************************************************************************

  nx = v*(ay*cz - az*cy) - u*(ay*bz - az*by) + by*cz - bz*cy
  ny = v*(az*cx - ax*cz) - u*(az*bx - ax*bz) + bz*cx - bx*cz
  nz = v*(ax*cy - ay*cx) - u*(ax*by - ay*bx) + bx*cy - by*cx        

  imag = 1.0_RFREAL/SQRT(nx*nx + ny*ny + nz*nz)

  nx = imag*nx
  ny = imag*ny
  nz = imag*nz

! ******************************************************************************
! End
! ******************************************************************************  

  CALL DeregisterFunction(global)  

END SUBROUTINE RFLU_BLIN_ComputeNormal










! ******************************************************************************
!
! Purpose: Compute intersection distance of given location along given line 
!   vector and bilinear face. 
!
! Description: None.
!
! Input:
!   pRegion             Pointer to region data
!   xLoc,yLoc,zLoc      x-, y-, and z-coordinates of location in question
!   ex,ey,ez            x-, y-, and z-components of unit line vector
!   iPatch              Face location 
!   ifg                 Face index
!
! Output:
!   nt                  Number of intersections
!   t                   Array of intersection distances
!
! Notes: 
!   1. The line vector MUST be a unit line vector. If that is not correct, the
!      distance between the given point and the intersection of the line with 
!      the faces of the given cell will not be computed correctly.
!   2. For details on derivation, see: S.D. Ramsey, K. Potter, and C. Hansen, 
!      Ray Bilinear Patch Intersection, ACM J. Graphics Tools, Vol. 9, No. 3, 
!      pp. 41-47, 2004. 
!
! ******************************************************************************

SUBROUTINE RFLU_BLIN_ComputeXSectLine(pRegion,xLoc,yLoc,zLoc,ex,ey,ez,icg, & 
                                      iPatch,ifg,nt,t)
     
  USE ModTools, ONLY: SwapRFREALs   
     
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  INTEGER, INTENT(IN) :: icg,ifg,iPatch
  INTEGER, INTENT(INOUT) :: nt
  REAL(RFREAL), INTENT(IN) :: ex,ey,ez
  REAL(RFREAL), INTENT(IN) :: xLoc,yLoc,zLoc    
  REAL(RFREAL), INTENT(OUT) :: t(2)
  TYPE(t_region), POINTER :: pRegion

! ==============================================================================
! Locals
! ==============================================================================

  INTEGER :: iv,nv,v1g,v1l,v2g,v2l,v3g,v3l,v4g,v4l
  REAL(RFREAL) :: a,ae,ax,ay,az,A1,A2,b,be,bx,by,bz,B1,B2,c,ce,cx,cy,cz,C1, &
                  C2,d,de,denom1,denom2,disc,dotp,dx,dy,dz,D1,D2, &
                  fuzzyTolerance,px,py,pz,q,re,u,xCofg,x1,x2,x3,x4, &
                  yCofg,y1,y2,y3,y4,zCofg,z1,z2,z3,z4
  REAL(RFREAL), DIMENSION(2) :: v
  TYPE(t_global), POINTER :: global
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_patch), POINTER :: pPatch

! DEBUG
  INTEGER :: ivg
  REAL(RFREAL) :: nx,ny,nz,A3,B3,C3,D3
! END DEBUG

! ******************************************************************************
! Start
! ******************************************************************************

  global => pRegion%global

  CALL RegisterFunction(global,'RFLU_BLIN_ComputeXSectLine',&
  'RFLU_ModBilinearPatch.F90')

! ==============================================================================
! Set grid pointer and initialize variables
! ==============================================================================

  pGrid => pRegion%grid

  nt = 0
  
  fuzzyTolerance = -pRegion%mixtInput%tolerICT

  t(1:2) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)

! ******************************************************************************
! Get face vertices and coordinates
! ******************************************************************************

! DEBUG
  WRITE(0,*) '###100',pRegion%iRegionGlobal,iPatch,ifg
  WRITE(0,*) '###110',pRegion%iRegionGlobal,xLoc,yLoc,zLoc
  WRITE(0,*) '###120',pRegion%iRegionGlobal,ex,ey,ez
  
!  DO ivg = 1,pGrid%nVert
!    WRITE(*,'(1X,I2,3(1X,E13.6))') ivg,pGrid%xyz(:,ivg)
!  END DO ! ivg
! END DEBUG  

  IF ( iPatch == 0 ) THEN 
    v1g = pGrid%f2v(1,ifg)
    v2g = pGrid%f2v(2,ifg)
    v3g = pGrid%f2v(3,ifg)
    v4g = pGrid%f2v(4,ifg)  
    
    IF ( v4g == VERT_NONE ) THEN 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! v4g        
  ELSE IF ( iPatch > 0 ) THEN 
    pPatch => pRegion%patches(iPatch)
    
    v1l = pPatch%bf2v(1,ifg)
    v2l = pPatch%bf2v(2,ifg)
    v3l = pPatch%bf2v(3,ifg)
    v4l = pPatch%bf2v(4,ifg) 

    IF ( v4l == VERT_NONE ) THEN 
      CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
    END IF ! v4l
      
    v1g = pPatch%bv(v1l)
    v2g = pPatch%bv(v2l)
    v3g = pPatch%bv(v3l)
    v4g = pPatch%bv(v4l)                    
  ELSE
    CALL ErrorStop(global,ERR_REACHED_DEFAULT,__LINE__)
  END IF ! iPatch
  
! DEBUG
  WRITE(0,*) '###130',pRegion%iRegionGlobal,v1g,v2g,v3g,v4g
! END DEBUG  
  
  x1 = pGrid%xyz(XCOORD,v1g)
  x2 = pGrid%xyz(XCOORD,v2g)
  x3 = pGrid%xyz(XCOORD,v3g)
  x4 = pGrid%xyz(XCOORD,v4g)      

  y1 = pGrid%xyz(YCOORD,v1g)
  y2 = pGrid%xyz(YCOORD,v2g)
  y3 = pGrid%xyz(YCOORD,v3g)
  y4 = pGrid%xyz(YCOORD,v4g) 
  
  z1 = pGrid%xyz(ZCOORD,v1g)
  z2 = pGrid%xyz(ZCOORD,v2g)
  z3 = pGrid%xyz(ZCOORD,v3g)
  z4 = pGrid%xyz(ZCOORD,v4g)   

! DEBUG
  WRITE(0,*) '@@@131a',pRegion%iRegionGlobal,x1,x2,x3,x4
  WRITE(0,*) '@@@131b',pRegion%iRegionGlobal,y1,y2,y3,y4
  WRITE(0,*) '@@@131c',pRegion%iRegionGlobal,z1,z2,z3,z4
! END DEBUG

! ******************************************************************************
! Compute geometric terms
! ******************************************************************************

  re = xLoc*ex + yLoc*ey + zLoc*ez

  ax = x3 - x2 - x4 + x1
  ay = y3 - y2 - y4 + y1
  az = z3 - z2 - z4 + z1  
  ae = ax*ex + ay*ey + az*ez

  bx = x2 - x1
  by = y2 - y1
  bz = z2 - z1   
  be = bx*ex + by*ey + bz*ez 
  
  cx = x4 - x1
  cy = y4 - y1
  cz = z4 - z1    
  ce = cx*ex + cy*ey + cz*ez

  dx = x1
  dy = y1
  dz = z1        
  de = dx*ex + dy*ey + dz*ez

  A1 = az - ay - ae*(ez - ey)
  A2 = az - ax - ae*(ez - ex)
  A3 = A2 - A1
  
  B1 = bz - by - be*(ez - ey)
  B2 = bz - bx - be*(ez - ex)
  B3 = B2 - B1
  
  C1 = cz - cy - ce*(ez - ey)
  C2 = cz - cx - ce*(ez - ex)
  C3 = C2 - C1
  
  D1 = dz - dy - (zLoc - yLoc) - (de - re)*(ez - ey)
  D2 = dz - dx - (zLoc - xLoc) - (de - re)*(ez - ex)
  D3 = D2 - D1

  a = A1*C2 - A2*C1
  b = A1*D2 - A2*D1 + B1*C2 - B2*C1
  c = B1*D2 - B2*D1
  
! DEBUG
  WRITE(0,*) 'a:',ax,ay,az,ae
  WRITE(0,*) 'b:',bx,by,bz,be
  WRITE(0,*) 'c:',cx,cy,cz,ce
  WRITE(0,*) 'd:',dx,dy,dz,de
  WRITE(0,*) 'A-D1:',A1,B1,C1,D1
  WRITE(0,*) 'A-D2:',A2,B2,C2,D2
  WRITE(0,*) 'A-D3:',A3,B3,C3,D3
  WRITE(0,*) 'a-c:',a,b,c
! END DEBUG  
  
! ******************************************************************************
! Compute v 
! ******************************************************************************
  
  IF ( ABS(a) > 0.0_RFREAL ) THEN ! Quadratic equation   
    disc = b*b - 4.0_RFREAL*a*c

! DEBUG
  WRITE(0,*) 'disc=',disc
! END DEBUG  

    IF ( disc > 0.0_RFREAL ) THEN
      nv = 2
    
      IF ( ABS(c) > 0.0_RFREAL ) THEN 
        q = -0.5_RFREAL*(b + SIGN(1.0_RFREAL,b)*SQRT(disc))

        v(1) = q/a
        v(2) = c/q
      ELSE 
        v(1) = 0.0_RFREAL
        v(2) = -b/a 
      END IF ! nv
    ELSE IF ( disc < 0.0_RFREAL ) THEN 
      nv = 0

      v(1) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)  
      v(2) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)    
    ELSE 
      nv = 1

      v(1) = -0.5_RFREAL*b/a
      v(2) = REAL(CRAZY_VALUE_INT,KIND=RFREAL)
    END IF ! disc
  ELSE IF ( ABS(b) > 0.0_RFREAL ) THEN ! Linear equation
    nv = 1
    
    v(1) = -c/b
  ELSE 
    nv = 0
  END IF ! ABS(a)

! DEBUG
  WRITE(0,*) '###132',pRegion%iRegionGlobal,v(:)
! END DEBUG  
  

! ******************************************************************************
! Compute intersection distance
! ******************************************************************************

  DO iv = 1,nv

! ==============================================================================
!   Check value of v, proceed only if within range [0,1]
! ==============================================================================

    IF ( (v(iv) < 0.0_RFREAL) .OR. (v(iv) > 1.0_RFREAL) ) THEN 
      CYCLE
    END IF ! v(iv)
    
! ==============================================================================
!   Compute u, proceed only if within range [0,1]
! ==============================================================================

    denom1 = v(iv)*A2 + B2
    denom2 = denom1 - (v(iv)*A1 + B1)
    
    IF ( ABS(denom2) > ABS(denom1) ) THEN 
      u = (v(iv)*(C1 - C2) + D1 - D2)/denom2
    ELSE 
      u = -(v(iv)*C2 + D2)/denom1
    END IF ! ABS(denom1)
    
! DEBUG
  WRITE(0,*) '###134',pRegion%iRegionGlobal,iv,v(iv),u
! END DEBUG  
    
    IF ( (u < 0.0_RFREAL) .OR. (u > 1.0_RFREAL) ) THEN 
      CYCLE
    END IF ! u
    
! ==============================================================================
!   Compute intersection point 
! ==============================================================================
   
    px = u*v(iv)*ax + u*bx + v(iv)*cx + dx
    py = u*v(iv)*ay + u*by + v(iv)*cy + dy    
    pz = u*v(iv)*az + u*bz + v(iv)*cz + dz

! DEBUG
  WRITE(0,*) '###135',pRegion%iRegionGlobal,px,py,pz
! END DEBUG  
    
! ==============================================================================
!   Compute normal at intersection point and ensure it is an outward normal
! ==============================================================================
    
    CALL RFLU_BLIN_ComputeNormal(global,ax,ay,az,bx,by,bz,cx,cy,cz, &
                                 dx,dy,dz,u,v(iv),nx,ny,nz)
  
    xCofg = pGrid%cofg(XCOORD,icg)
    yCofg = pGrid%cofg(YCOORD,icg)
    zCofg = pGrid%cofg(ZCOORD,icg)  
  
    dotp = (px - xCofg)*nx + (py - yCofg)*ny + (pz - zCofg)*nz
  
! DEBUG
  WRITE(0,*) '###136',pRegion%iRegionGlobal,nx,ny,nz,dotp
! END DEBUG  
    
    IF ( dotp < 0.0_RFREAL ) THEN
      nx = -nx
      ny = -ny
      nz = -nz
    END IF ! dotp
    
! ==============================================================================
!   Compute intersection distance. NOTE expression for intersection
!   distance differs from that in the paper - it avoids IFs and divisions.
! ==============================================================================
  
    dotp = (px - xLoc)*nx + (py - yLoc)*ny + (pz - zLoc)*nz
  
! DEBUG
  WRITE(0,*) '###137',pRegion%iRegionGlobal,dotp
! END DEBUG  
    
    IF ( dotp > fuzzyTolerance ) THEN ! Inside cell 
! DEBUG
  WRITE(0,*) '###138',pRegion%iRegionGlobal,ex*nx + ey*ny + ez*nz
! END DEBUG  
    
      IF ( ex*nx + ey*ny + ez*nz > 0.0_RFREAL ) THEN ! In direction of traj
        nt = nt + 1
    
        t(nt) = (px - xLoc)*ex + (py - yLoc)*ey + (pz - zLoc)*ez
      END IF ! ex*nx
    END IF ! dotp
  END DO ! iv
  
! ******************************************************************************
! Swap so first entry is that with smaller intersection distance
! ******************************************************************************

  IF ( nt > 0 ) THEN 
    IF ( ABS(t(2)) < ABS(t(1)) ) THEN 
      CALL SwapRFREALs(t(1),t(2))
    END IF ! ABS(t(2)) 
  END IF ! nt

! DEBUG
  WRITE(0,*) '###140',pRegion%iRegionGlobal,nt,t(1:2)
! END DEBUG  
  
! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_BLIN_ComputeXSectLine








! ******************************************************************************
!
! Purpose: Given a point, compute closest point on bilinear patch.
!
! Description: None.
!
! Input:
!   global      Pointer to global data
!   ax,ay,az    Components of bilinear patch parametric representation
!   bx,by,bz    Components of bilinear patch parametric representation
!   cx,cy,cz    Components of bilinear patch parametric representation
!   dx,dy,dz    Components of bilinear patch parametric representation
!   xq,yq,zq    Coordinates of point
!
! Output: 
!   u,v         Parametric coordinates of closest point on bilinear patch
!   x,y,z       Coordinates of closest point on bilinear patch
!
! Notes:
!   1. Use Newton-Raphson method to solve non-linear system for (u,v). 
!   2. This routine does NOT check that the closest point lies in the patch, 
!      i.e., that the (u,v) coordinates are in the range [0,1].
!
! ******************************************************************************

  SUBROUTINE RFLU_BLIN_FindClosestPoint(global,ax,ay,az,bx,by,bz,cx,cy,cz, & 
                                        dx,dy,dz,xq,yq,zq,u,v,x,y,z)

    IMPLICIT NONE

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

! ==============================================================================
!   Parameters      
! ==============================================================================

    REAL(RFREAL), INTENT(OUT) :: u,v,x,y,z
    REAL(RFREAL), INTENT(IN) :: ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,xq,yq,zq
    TYPE(t_global), POINTER :: global   
    
! ==============================================================================
!   Locals      
! ==============================================================================

    INTEGER :: loopCounter
    REAL(RFREAL) :: dfdu,dfdv,dgdu,dgdv,du,dv,dxdu,dxdv,dydu,dydv,dzdu,dzdv, &
                    d2,d2xduv,d2yduv,d2zduv,f,g,idet,tol
    
! ******************************************************************************
!   Start
! ******************************************************************************

    CALL RegisterFunction(global,'RFLU_BLIN_FindClosestPoint',&
  'RFLU_ModBilinearPatch.F90')


! ******************************************************************************
!   Initialize
! ******************************************************************************

    loopCounter = 0
    
    tol = 1.0E-6_RFREAL

    u = 0.5_RFREAL ! Initial guess
    v = 0.5_RFREAL

! ******************************************************************************
!   Newton-Raphson iteration to compute closest point
! ******************************************************************************
        
    emptyLoop: DO 
      loopCounter = loopCounter + 1
    
      x = u*v*ax + u*bx + v*cx + dx 
      y = u*v*ay + u*by + v*cy + dy
      z = u*v*az + u*bz + v*cz + dz            
    
      dxdu = v*ax + bx 
      dydu = v*ay + by
      dzdu = v*az + bz  
      
      dxdv = u*ax + cx
      dydv = u*ay + cy
      dzdv = u*az + cz                  
    
      d2xduv = ax
      d2yduv = ay
      d2zduv = az
    
      d2 = (x-xq)*(x-xq) + (y-yq)*(y-yq) + (z-zq)*(z-zq) 
    
      f = (x-xq)*dxdu + (y-yq)*dydu + (z-zq)*dzdu
      g = (x-xq)*dxdv + (y-yq)*dydv + (z-zq)*dzdv
      
      dfdu = dxdu*dxdu + dydu*dydu + dzdu*dzdu
      dfdv = dxdu*dxdv + dydu*dydv + dzdu*dzdv & 
           + (x-xq)*d2xduv + (y-yq)*d2yduv + (z-zq)*d2zduv 

      dgdu = dfdv 
      dgdv = dxdv*dxdv + dydv*dydv + dzdv*dzdv
      
      idet = 1.0_RFREAL/(dfdu*dgdv - dfdv*dgdu)
      
      du = -idet*(dgdv*f - dgdu*g)
      dv = -idet*(dfdu*g - dfdv*f)

      IF ( ABS(du) < tol .AND. ABS(dv) < tol ) THEN 
        EXIT emptyLoop
      ELSE IF ( loopCounter > LIMIT_INFINITE_LOOP ) THEN 
        CALL ErrorStop(global,ERR_INFINITE_LOOP,__LINE__)
      END IF ! ABS(du) 

      u = u + du 
      v = v + dv      
    END DO emptyLoop 

! ******************************************************************************
!   End
! ******************************************************************************  

    CALL DeregisterFunction(global)  

  END SUBROUTINE RFLU_BLIN_FindClosestPoint






! ******************************************************************************
! End
! ******************************************************************************
  
END MODULE RFLU_ModBilinearPatch


! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_ModBilinearPatch.F90,v $
! Revision 1.5  2008/12/06 08:44:20  mtcampbe
! Updated license.
!
! Revision 1.4  2008/11/19 22:17:31  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.3  2006/04/07 15:19:18  haselbac
! Removed tabs
!
! Revision 1.2  2006/03/25 21:49:14  haselbac
! Renamed SwapFloats to SwapRFREALs
!
! Revision 1.1  2005/12/24 21:17:50  haselbac
! Initial revision
!
! ******************************************************************************
  









