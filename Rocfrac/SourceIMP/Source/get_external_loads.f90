!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
SUBROUTINE get_external_loads(coords,jdltyp, &
     adlmag,mdload,ndload,Rex,nodes,NumNp,&
     MapSElVolEl,NumElVol)
!
! *----------------------------------------------GET_EXTERNAL_LOADS----*
! |                                                                    |
! |   Add surface tractions and Add body forces                        |
! |                                                                    |
! |   Notes (for details, see report):                                 |
! |                                                                    |
! |   1- Non-Uniform surface fluxes are supported                      |
! |                                                                    |
! |   2- The face that the load is applied is designated by            |
! |      load type in the |DLOAD card. For example :                   |
! |          *DLOAD, OP=NEW                                            |
! |          SET_NAME, U1, 1.0                                         |
! |      describes a uniform pressure of magnitude 1.0                 |
! |      to be applied to the face 1 (number taken from U1)            |
! |      of all the elements in the set SET_NAME.                      |
! |                                                                    |
! |   3- Load types mean,                                              |
! |          U1,   pressure onto the 1-2-3-4 face                      |
! |          U2,   pressure onto the 5-8-7-6 face                      |
! |          U3,   pressure onto the 1-5-6-2 face                      |
! |          U4,   pressure onto the 2-6-7-3 face                      |
! |          U5,   pressure onto the 3-7-8-4 face                      |
! |          U6,   pressure onto the 4-8-5-1 face                      |
! |          U7,   body force in the global X-direction                |
! |          U8,   body force in the global Y-direction                |
! |          U9,   body force in the global Z-direction                |
! |                                                                    |
! *--------------------------------------------------------------------*

  IMPLICIT DOUBLE PRECISION (a-h, o-z)
      
!  Argument variables

  INTEGER :: NumNp
  INTEGER, DIMENSION(1:8,1:NumNp) :: nodes
  REAL*8, DIMENSION(1:3*NumNp) :: Rex
  REAL*8, DIMENSION(1:3,1:NumNp) :: Coords
  INTEGER, DIMENSION(1:NumElVol) :: MapSFElVolEl
  DIMENSION  jdltyp(mdload,1:2),adlmag(mdload)

  INTEGER :: IdVolEl

!  Local variables

  REAL*8, DIMENSION(1:3) :: veca, vecb, outwd_surface_normal, &
       half_normal_veca, half_normal_vecb

  INTEGER :: iload
  INTEGER :: ltype
  INTEGER :: nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8

  REAL*8, DIMENSION(1:24) :: dlvec
  INTEGER :: icnt

!  Data statements

  DATA eigth /0.125000000000000/

!  Compute external loads and add to the load vector

  DO iload = 1,         ! LOOP OVER LOADS

     dlvec(:) = 0.d0

     IdVolEl = MapSFElVolEl(i) ! Volumetric Element Id
     ltype = jdltyp(iload,2) ! get load type

     nd1 = nodes(1,IdVolEl)
     nd2 = nodes(2,IdVolEl)
     nd3 = nodes(3,IdVolEl)
     nd4 = nodes(4,IdVolEl)
     nd5 = nodes(5,IdVolEl)
     nd6 = nodes(6,IdVolEl)
     nd7 = nodes(7,IdVolEl)
     nd8 = nodes(8,IdVolEl)
         
     IF(ltype.LT.0) THEN    ! unsupported load type (UnNU)
        WRITE(6,1100)
        STOP
     END IF

! branch depending on load type
     SELECT CASE(ltype)

     CASE (100) ! pressure on Face 1

        veca(1:3) = coords(1:3,nd4) - coords(1:3,nd1)
        vecb(1:3) = coords(1:3,nd2) - coords(1:3,nd1)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd2) - coords(1:3,nd3)
        vecb(1:3) = coords(1:3,nd4) - coords(1:3,nd3)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb
        
        pmag = -eigth * adlmag(iload)   
        
        dlvec(1:3)   = pmag * outwd_surface_normal(1:3) ! node 1
        dlvec(4:6)   = pmag * outwd_surface_normal(1:3) ! node 2
        dlvec(7:9)   = pmag * outwd_surface_normal(1:3) ! node 3
        dlvec(10:12) = pmag * outwd_surface_normal(1:3) ! node 4
         
     CASE(200) ! pressure on Face 2
         
        veca(1:3) = coords(1:3,nd6) - coords(1:3,nd5)
        vecb(1:3) = coords(1:3,nd8) - coords(1:3,nd5)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd8) - coords(1:3,nd7)
        vecb(1:3) = coords(1:3,nd6) - coords(1:3,nd7)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb
        
        pmag = -eigth * adlmag(iload)
        
        dlvec(13:15) = pmag * outwd_surface_normal(1:3) ! node 5
        dlvec(16:18) = pmag * outwd_surface_normal(1:3) ! node 6
        dlvec(19:21) = pmag * outwd_surface_normal(1:3) ! node 7
        dlvec(22:24) = pmag * outwd_surface_normal(1:3) ! node 8
         
     CASE(300) ! pressure on Face 3
        
        veca(1:3) = coords(1:3,nd2) - coords(1:3,nd1)
        vecb(1:3) = coords(1:3,nd5) - coords(1:3,nd1)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd5) - coords(1:3,nd6)
        vecb(1:3) = coords(1:3,nd2) - coords(1:3,nd6)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb
        
        pmag = -eigth * adlmag(iload)
        
        dlvec(1:3)   = pmag * outwd_surface_normal(1:3) ! node 1
        dlvec(4:6)   = pmag * outwd_surface_normal(1:3) ! node 2
        dlvec(13:15) = pmag * outwd_surface_normal(1:3) ! node 5
        dlvec(16:18) = pmag * outwd_surface_normal(1:3) ! node 6
         
     CASE(400) ! pressure on Face 4
         
        veca(1:3) = coords(1:3,nd3) - coords(1:3,nd2)
        vecb(1:3) = coords(1:3,nd6) - coords(1:3,nd2)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd6) - coords(1:3,nd7)
        vecb(1:3) = coords(1:3,nd3) - coords(1:3,nd7)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb
        pmag = -eigth * adlmag(iload)
         
        dlvec(4:6)   = pmag * outwd_surface_normal(1:3) ! node 2
        dlvec(7:9)   = pmag * outwd_surface_normal(1:3) ! node 3
        dlvec(16:18) = pmag * outwd_surface_normal(1:3) ! node 6
        dlvec(19:21) = pmag * outwd_surface_normal(1:3) ! node 7
         

     CASE(500) ! pressure on Face 5
         
        veca(1:3) = coords(1:3,nd4) - coords(1:3,nd3)
        vecb(1:3) = coords(1:3,nd7) - coords(1:3,nd3)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd7) - coords(1:3,nd8)
        vecb(1:3) = coords(1:3,nd4) - coords(1:3,nd8)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb
        pmag = -eigth * adlmag(iload)
        
        dlvec(7:9)   = pmag * outwd_surface_normal(1:3) ! node 3
        dlvec(10:12) = pmag * outwd_surface_normal(1:3) ! node 4
        dlvec(19:21) = pmag * outwd_surface_normal(1:3) ! node 7
        dlvec(22:24) = pmag * outwd_surface_normal(1:3) ! node 8

     CASE(600) ! pressure on Face 6

        veca(1:3) = coords(1:3,nd5) - coords(1:3,nd1)
        vecb(1:3) = coords(1:3,nd4) - coords(1:3,nd1)
        
        half_normal_veca(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_veca(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_veca(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        veca(1:3) = coords(1:3,nd4) - coords(1:3,nd8)
        vecb(1:3) = coords(1:3,nd5) - coords(1:3,nd8)
        
        half_normal_vecb(1) = veca(2)*vecb(3) - vecb(2)*veca(3)
        half_normal_vecb(2) = veca(3)*vecb(1) - vecb(3)*veca(1)
        half_normal_vecb(3) = veca(1)*vecb(2) - vecb(1)*veca(2)
        
        outwd_surface_normal = half_normal_veca + half_normal_vecb

        pmag = -eigth * adlmag(iload)
        
        dlvec(1:3)   = pmag * outwd_surface_normal(1:3) ! node 1
        dlvec(10:12) = pmag * outwd_surface_normal(1:3) ! node 4
        dlvec(13:15) = pmag * outwd_surface_normal(1:3) ! node 5
        dlvec(22:24) = pmag * outwd_surface_normal(1:3) ! node 8
         
     CASE(700) ! body force in X-direction

     CASE(800) ! body force in Y-direction

     CASE(900)! body force in Z-direction



     END SELECT 

     icnt = 1
     DO i = 1, 8
        nd = nodes(i,IdVolEl)
        Rex(3*nd-2) = Rex(3*nd-2) + dlvec(icnt)
        Rex(3*nd-1) = Rex(3*nd-1) + dlvec(icnt+1)
        Rex(3*nd  ) = Rex(3*nd  ) + dlvec(icnt+2)
        icnt = icnt + 3
     ENDDO

  END DO                    ! END LOOP OVER LOADS

1100 FORMAT (/,2x,'>>> Non-uniform surface tractions are not', &
          /,2x,'    supported (i.e. UnNU). Aborting !')

END SUBROUTINE get_external_loads

