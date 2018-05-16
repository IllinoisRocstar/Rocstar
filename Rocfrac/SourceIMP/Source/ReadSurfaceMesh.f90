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
subroutine ReadSufaceMesh(glb,NumProcs,iunit,iElType2D,myid)


  integer :: iunit

  integer :: iaux,k,kk,id2d
  
  READ(15,*) iaux, glb%iElType2D
  DO i = 1, NumProcs
     READ(iunit,*) id2d ! processor id
     IF(id2d.EQ.myid+1)THEN
!
! Read No Solid/Fluid mesh
!
        READ(iunit,*) glb%InterfaceSNumNodes, glb%InterfaceSNumElems
        ALLOCATE(glb%InterfaceSNodalCoors(1:3,1:glb%InterfaceSNumNodes))
        ALLOCATE(glb%MapNodeS(1:glb%InterfaceSNumNodes))
          
          DO j = 1, glb%InterfaceSNumNodes
             READ(iunit,*) glb%InterfaceSNodalCoors(1,j), &
                        glb%InterfaceSNodalCoors(2,j), &
                        glb%InterfaceSNodalCoors(3,j), glb%MapNodeS(j)
          ENDDO
          
          ALLOCATE(glb%InterfaceSElemConn(1:glb%iElType2D,1:glb%InterfaceSNumElems))
          IF(glb%EnforceTractionS)THEN
             ALLOCATE(glb%MapSElVolEl(1:glb%InterfaceSNumElems))
             DO j = 1, glb%InterfaceSNumElems
                IF ( glb%iElType2D==3) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j),&
                        glb%MapSElVolEl(j)
                ELSE IF ( glb%iElType2D==6) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j),&
                        glb%InterfaceSElemConn(4,j),&
                        glb%InterfaceSElemConn(5,j),&
                        glb%InterfaceSElemConn(6,j),&
                        glb%MapSElVolEl(j)
                ELSE IF ( glb%iElType2D==4) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j),&
                        glb%InterfaceSElemConn(4,j),&
                        glb%MapSElVolEl(j)
                END IF
             ENDDO
             EXIT
          ELSE
             DO j = 1, glb%InterfaceSNumElems
                IF ( glb%iElType2D==3) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j)
                ELSE IF ( glb%iElType2D==6) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j),&
                        glb%InterfaceSElemConn(4,j),&
                        glb%InterfaceSElemConn(5,j),&
                        glb%InterfaceSElemConn(6,j)
                ELSE IF ( glb%iElType2D==4) THEN
                   READ(iunit,*) glb%InterfaceSElemConn(1,j),&
                        glb%InterfaceSElemConn(2,j),&
                        glb%InterfaceSElemConn(3,j),&
                        glb%InterfaceSElemConn(4,j)
                END IF
             ENDDO
             EXIT
          ENDIF
       ELSE ! Not the processor's surface mesh
          READ(iunit,*) k, kk
          DO j = 1, k
             READ(iunit,'()')
          ENDDO
          DO j = 1, kk
             READ(iunit,'()')
          ENDDO
       ENDIF
    ENDDO

