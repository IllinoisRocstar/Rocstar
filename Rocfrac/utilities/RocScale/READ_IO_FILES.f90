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
SUBROUTINE READ_IO_FILES
  use data_declarations
  
  IMPLICIT NONE
  
  CHARACTER*50 :: ichr50
  
  integer :: IDpacket

  OPEN(UNIT = 14, FILE=PREFIX(1:LENGTH(PREFIX))//".inp", STATUS="OLD")
  
  
  ichr50 = 'mkdir '//PREFIX(1:LENGTH(PREFIX))
  CALL system(ichr50)

      ! -----READ LOAD DATA
!      READ(14,*) numnp, numel          ! number of global nodes and elements              
!
!      READ(14,*) delta
!      READ(14,*) numprop
!      ALLOCATE(ts_proportion(1:numprop),proportion(1:numprop))
!      DO l = 1,numprop
!         READ(14,*) ts_proportion(l), proportion(l)
!      ENDDO

  ALLOCATE(frontnode); NULLIFY(frontnode%next)
  first_front_node => frontnode; numscale_np = 0
  
  ALLOCATE(backnode); NULLIFY(backnode%next)
  first_back_node => backnode; numscale_np = 0
  
  numscale_bk = 0
  numscale_ft = 0
  maxZ = -10000
  minZ = 10000

  DO
     read(14,*) IDpacket

     print*,IDpacket
     
     select case (IDpacket)
     
     case(1)
        read(14,'()')
        
     case(2)
     
      !-----READ NODAL COORDINATES.
        READ(14,*) numnp  !, iaux, iaux, iaux, iaux
        
        ALLOCATE(node(numnp))
      
      ! This initially sets all nodes to have a boundary flag for internal
        node(1:numnp)%marker = 0
        node(1:numnp)%meshmark = 0
        
        ALLOCATE(meshcoor(1:3,1:numnp))
        ALLOCATE(nboundtype(1:numnp))
        DO l=1,numnp
           READ(14,*) j,node(l)%coord(1:3) ! ,iaux
           
           ! Defines the maximum and minimum Z coord values for matching nodes
           maxZ = MAX(maxZ, node(l)%coord(3))
           minZ = MIN(minZ, node(l)%coord(3))
           
           !  new ! for ALE formulation
           meshcoor(1:3,l) = node(l)%coord(1:3)
           
        END DO
        
        WRITE(6,*) "MinZ = ",minZ, "MaxZ = ", maxZ
        ! Reassigns maxZ & minZ to offset possible discrepencies in node coords
        maxZ = maxZ - .00001*maxZ
        minZ = minZ + .00001*minZ
        
        
        
        DO l = 1, numnp
           
           IF (node(l)%coord(3) .LE. minZ) THEN
              backnode%ID = l
              backnode%coord(1:3) = node(l)%coord(1:3)
              ALLOCATE(backnode%next); backnode => backnode%next; NULLIFY(backnode%next)
              numscale_np = numscale_np + 1
              numscale_bk = numscale_bk + 1
              
           ELSE IF (node(l)%coord(3) .GE. maxZ) THEN
              frontnode%ID = l
              frontnode%coord(1:3) = node(l)%coord(1:3)
              ALLOCATE(frontnode%next); frontnode => frontnode%next; NULLIFY(frontnode%next)
              numscale_np = numscale_np + 1
              numscale_ft = numscale_ft + 1
              
              !  new ! for ALE formulation
              meshcoor(1:3,l) = node(l)%coord(1:3)
           END IF
        END DO
        
     case(3)
        
        !---------------------------------!
        !--- Read Nodal Boundary Flags ---!
        !---------------------------------!
        
        frontnode => first_front_node
        backnode => first_back_node
        
        READ(14,*) numbound ! , iaux
        ALLOCATE(id(1:numbound))
        
        numboundmid = numbound
        numboundend = numbound
        numboundfirst = numbound
        
        DO l = 1, numbound
           READ(14,*) id(l)%ID,id(l)%NdBCflag ! , iaux
           ! Apply marker information to node data struct
           node(id(l)%ID)%marker = id(l)%NdBCflag
           IF((id(l)%NdBCflag .GE. 10) .AND. (id(l)%NdBCflag .LT. 100)) THEN
              numboundmid = numboundmid - 1
              numboundend = numboundend - 1
              DO WHILE(associated(backnode%next))
                 IF (backnode%ID .EQ. id(l)%ID) THEN
                    backnode%marker = id(l)%NdBCflag
                 ENDIF
                 backnode => backnode%next
              ENDDO
              
           ELSE IF (id(l)%NdBCflag .GE. 100) THEN
              numboundmid = numboundmid - 1
              numboundfirst = numboundfirst - 1
              DO WHILE(associated(frontnode%next))
                 IF (frontnode%ID .EQ. id(l)%ID) THEN
                    frontnode%marker = id(l)%NdBCflag
                 ENDIF
                 frontnode => frontnode%next
              ENDDO
           END IF
           frontnode => first_front_node
           backnode => first_back_node
        ENDDO
        
     case(4)

        !---------------------------------------------!
        !--- Read Nodal Mesh Motion Boundary Flags ---!
        !---------------------------------------------!
        
        READ(14,*) numboundmesh !, iaux
        ALLOCATE(idmesh(1:numboundmesh))  !,rmesh(1:3,1:numboundmesh))
        
        numbmeshmid = numboundmesh
        numbmeshend = numboundmesh
        numbmeshfirst = numboundmesh
        
        DO l = 1, numboundmesh
           READ(14,*) idmesh(l)%ID, idmesh(l)%NdBCflag !, iaux
           ! Apply marker information to node data struct
           node(idmesh(l)%ID)%meshmark = idmesh(l)%NdBCflag
           IF((idmesh(l)%NdBCflag .GE. 10) .AND. (idmesh(l)%NdBCflag .LT. 100)) THEN
              numbmeshmid = numbmeshmid - 1
              numbmeshend = numbmeshend - 1
              DO WHILE(associated(backnode%next))
                 IF (backnode%ID .EQ. idmesh(l)%ID) THEN
                    backnode%mesh = idmesh(l)%NdBCflag
                 ENDIF
                 backnode => backnode%next
              ENDDO
              
           ELSE IF (idmesh(l)%NdBCflag .GE. 100) THEN
              numbmeshmid = numbmeshmid - 1
              numbmeshfirst = numbmeshfirst - 1
              DO WHILE(associated(frontnode%next))
                 IF (frontnode%ID .EQ. idmesh(l)%ID) THEN
                    frontnode%mesh = idmesh(l)%NdBCflag
                 ENDIF
                 frontnode => frontnode%next
              ENDDO
           END IF
           frontnode => first_front_node
           backnode => first_back_node
        ENDDO
        
        
      
     !-----Read cohesive element data, determine the length and angle
!      READ(14,*) !numco, numclst, numcohshared, nproc_neigh, num_border_coh
!      READ(14,*)

     case(5)

!      READ(14,*) numcstet,l,num_border_cst, priv1, priv2
        READ(14,*) numcstet,num_border_cst
        
        ALLOCATE(lmcstet(1:4,1:numcstet),matcstet(1:numcstet))
        
     !
     ! -- Read element data.
        ALLOCATE(elemlist); NULLIFY(elemlist%next); NULLIFY(elemlist%previous)
        ALLOCATE(eboundlist); NULLIFY(eboundlist%next); NULLIFY(eboundlist%previous)
        firstelem => elemlist
        firstebound => eboundlist
        
        NBoundryEl3D = 0
        
        DO l = 1,numcstet
           READ(14,*) matcstet(l), lmcstet(1:4,l)
           
           DO m = 1, numnp
              
              ! If there is a node on the boundary...
              IF((node(m)%coord(3) .GE. maxZ) .OR. (node(m)%coord(3) .LE. 0.000000)) THEN
                 
                 DO n = 1, 4
                    IF (m .EQ. lmcstet(n,l)) THEN
                       eboundlist%ID = l; eboundlist%conn(1:4) = lmcstet(1:4,l)
                       ALLOCATE(eboundlist%next); eboundlist => eboundlist%next; NULLIFY(eboundlist%next)
                       NBoundryEl3D = NBoundryEl3D + 1
                       
                       GOTO 50
                    END IF
                 END DO
                 
              END IF
              
           END DO
           elemlist%ID = l; elemlist%conn(1:4) = lmcstet(1:4,l)
           ALLOCATE(elemlist%next); elemlist => elemlist%next; NULLIFY(elemlist%next)
50      END DO
        PRINT*,'Finished Case 5'
     case (99)
        exit   
     end select
  enddo
  
!    DEBUG ROUTINE...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      OPEN(UNIT = 1000, FILE = "elemlist.dat", STATUS = "UNKNOWN")
!
!      elemlist => firstelem
!      DO WHILE(associated(elemlist%next))
!         write(1000,*) elemlist%ID, elemlist%conn(1:4)
!         elemlist => elemlist%next
!      END DO
!      DO WHILE(associated(elemlist%previous))
!         write(1000,*) elemlist%ID, elemlist%conn(1:4)
!        elemlist => elemlist%previous
!      END DO
!     
!      CLOSE(1000)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !
     ! -- Read coh nodes that need to be sent to other processors
     !
!      ALLOCATE(neigh(0:scale-1))
!      ALLOCATE(neigh_proc(1:nproc_neigh))

!      DO i = 0,scale-1
!         READ(14,*) iaux ! should always be zero
!      ENDDO

 !     READ(14,*) iaux ! should always be zero

     !
     ! -- Read lst nodes that need to be sent to other processors
     !

      CLOSE(14)

    END SUBROUTINE READ_IO_FILES

