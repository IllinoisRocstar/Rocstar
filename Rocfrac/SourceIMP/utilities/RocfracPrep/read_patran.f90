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
SUBROUTINE read_patran(numvertx2d,dhmin,nprocs)
  
  USE linked_list
  USE meshdata
  USE Linked_List2
  
  IMPLICIT NONE
  
  INTEGER :: numvertx2d
  REAL*8 :: dhmin

! Local
  REAL*8 :: xx,yy,zz,size1,size2,size3,size4,size5,size6
  REAL*8 :: dt_courant
  REAL*8, DIMENSION(1:6) :: DisFlagValue, TmpVal
  INTEGER :: numDisFlag
  
  INTEGER :: id
  INTEGER :: i, j
  INTEGER :: itype
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER :: iface
  REAL*8 :: value
  REAL*8 :: press
!-- Tempory holding array of partitioned metis  arrays
  INTEGER, ALLOCATABLE, DIMENSION(:) :: npart
  INTEGER :: edgecut
  INTEGER :: nprocs
  INTEGER :: ip
!
  INTEGER :: iaux,iaux1,iaux2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: imin,nmin,imax,nmax,ninc
  INTEGER, ALLOCATABLE, DIMENSION(:) :: imap
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lm4node
  INTEGER :: nnd1, nnd2, nnd3, nnd4, nnd5, nnd6, nnd7, nnd8, nnd9, nnd10
  INTEGER :: numnp4node
  integer, dimension(1:6) :: DisFlag
  INTEGER :: iFaceId, itmp
  
  
  INTEGER :: ProcId
  
  
  LOGICAL :: Hex8, Tet4, Tet10, fileExist
  LOGICAL :: fileExist1, fileExist2, fileExist3, fileExist4, fileExist5, fileExist6
  integer :: NumGeomBodies,iunit,ibody

  INTEGER, POINTER, DIMENSION(:) :: NumNP_glb,NumEL_glb,NumMat_glb

  INTEGER :: iElemTest, icnt, iNodeTest, MinEl, MaxEl,maxnd, minnd, iEL10

  INTEGER, POINTER, DIMENSION(:) :: epart_loc 

  INTEGER, POINTER, dimension(:) :: ProcDist
  
  TYPE(Link_Ptr_Type)  :: Link

  INTEGER :: AccumEl, AccumNd, AccumNd_old, AccumEl_old

  INTEGER :: ierror, iv, kc

  size1 = 0.d0
  size2 = 0.d0
  size3 = 0.d0
  size4 = 0.d0
  size5 = 0.d0
  size6 = 0.d0

  ierror = 0 
  
  Hex8  = .FALSE.
  Tet4  = .FALSE.
  Tet10 = .FALSE.
  OverlayMesh = .FALSE.
  InteractMesh = .FALSE.
  
  PRINT*,'MESH OPTION:'
  PRINT*,'  READING PATRAN MESH'
  PRINT*,' '
  
  itmp   = 0
  MaxNumBC_str = 0
  MaxNumBC_mm  = 0
  MaxNumBC_th  = 0

  NumNdHistory = 0
  
  ALLOCATE(NumElHex2D(1:5,1:nprocs))
  ALLOCATE(NumEltet2D(1:5,1:nprocs))
  
  NumElHex2D(:,:) = 0
  NumEltet2D(:,:) = 0

  fileExist1 = .FALSE.
  fileExist2 = .FALSE.
  fileExist3 = .FALSE.
  fileExist4 = .FALSE.
  fileExist5 = .FALSE.
  fileExist6 = .FALSE.

  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out', EXIST=fileExist1)
  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out', EXIST=fileExist2)
  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.3.out', EXIST=fileExist3)
  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.4.out', EXIST=fileExist4)
  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.5.out', EXIST=fileExist5)
  INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.6.out', EXIST=fileExist6)

  NumGeomBodies = 1

  IF(fileExist6)THEN
     NumGeomBodies = 6
     OPEN(106,FILE='Modin/'//prefx(1:prefx_lngth)//'.6.out')
     OPEN(105,FILE='Modin/'//prefx(1:prefx_lngth)//'.5.out')
     OPEN(104,FILE='Modin/'//prefx(1:prefx_lngth)//'.4.out')
     OPEN(103,FILE='Modin/'//prefx(1:prefx_lngth)//'.3.out')
     OPEN(102,FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out')
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE IF (fileExist5)THEN
     NumGeomBodies = 5
     OPEN(105,FILE='Modin/'//prefx(1:prefx_lngth)//'.5.out')
     OPEN(104,FILE='Modin/'//prefx(1:prefx_lngth)//'.4.out')
     OPEN(103,FILE='Modin/'//prefx(1:prefx_lngth)//'.3.out')
     OPEN(102,FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out')
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE IF (fileExist4)THEN
     NumGeomBodies = 4
     OPEN(104,FILE='Modin/'//prefx(1:prefx_lngth)//'.4.out')
     OPEN(103,FILE='Modin/'//prefx(1:prefx_lngth)//'.3.out')
     OPEN(102,FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out')
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE IF (fileExist3)THEN
     NumGeomBodies = 3
     OPEN(103,FILE='Modin/'//prefx(1:prefx_lngth)//'.3.out')
     OPEN(102,FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out')
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE IF (fileExist2) THEN
     NumGeomBodies = 2
     OPEN(102,FILE='Modin/'//prefx(1:prefx_lngth)//'.2.out')
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE IF (fileExist1) THEN
     OPEN(101,FILE='Modin/'//prefx(1:prefx_lngth)//'.1.out')
  ELSE
     fileExist = .FALSE.
     
     INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.pat', EXIST=fileExist)
     
     IF (.NOT. fileExist)THEN
        PRINT*,'File '//'Modin/'//prefx(1:prefx_lngth)//'.pat NOT found'
        PRINT*,'...Trying '//'Modin/'//prefx(1:prefx_lngth)//'.out'
        
        fileExist = .FALSE.
        INQUIRE (FILE='Modin/'//prefx(1:prefx_lngth)//'.out', EXIST=fileExist)
     
        IF (.NOT. fileExist)THEN
           PRINT*,'ERROR: File '//'Modin/'//prefx(1:prefx_lngth)//'.out NOT found'
           PRINT*,'ERROR: No Patran Neutral file Found, stopping'
           STOP
        ELSE
           PRINT*,'File '//'Modin/'//prefx(1:prefx_lngth)//'.out found'

           OPEN(101,file='Modin/'//prefx(1:prefx_lngth)//'.out', FORM='formatted')
        ENDIF
     ELSE
        OPEN(101,file='Modin/'//prefx(1:prefx_lngth)//'.pat', FORM='formatted')
     ENDIF

  ENDIF

  PRINT*,'No. of Mesh Segment Files Found =', NumGeomBodies    

  ALLOCATE(NumNP_glb(1:NumGeomBodies),NumEL_glb(1:NumGeomBodies),NumMat_glb(1:NumGeomBodies))

  DO ibody = 1, NumGeomBodies
     iunit = 100 + ibody


!
! - Packet Type 25: Title Card, Format(i2,8i8)
     READ(iunit,'()')
     READ(iunit,'()')

! -  Packet Type 26: Summary Data, Format(i2,8i8)
! -      26 ID IV KC N1 N2 N3 N4 N5
! -  N1 = number of nodes
! -  N2 = number of elements
! -  N3 = number of materials
! -  N4 = number of Element Properties
! -  N5 = number of Coordinate Frames

     READ(iunit,*) i,i,i,i,NumNP_glb(ibody),NumEL_glb(ibody),NumMat_glb(ibody)
  
     IF(NumMat_glb(ibody).LE.0) NumMat_glb(ibody) = 1

! - part of title card
     READ(iunit,'()')

  ENDDO

  numnp_prmry = SUM(NumNP_glb)
  numelv_prmry = SUM(NumEL_glb)
  NumMat = MAXVAL(NumMat_glb)
  
  PRINT*,' No. of elements = ', numelv_prmry
  PRINT*,' No. of nodes = ', numnp_prmry
  PRINT*,' No. of Materials = ', NumMat
  PRINT*,' No. of vertex=',numvertx

  iNodeTest = 0
  iElemTest = 0

  iaux = 0

  ALLOCATE(epart(1:numelv_prmry))

  AccumNd = 0

  AccumEl = 0

  DO ibody = 1, NumGeomBodies
     iunit = 100 + ibody
     numnp4node = 0
!
! --  Read Nodal coordinates

     IF(ibody.EQ.1)THEN
        ALLOCATE(coor(1:3,1:numnp_prmry),lmelv_prmry(1:numvertx,1:numelv_prmry))
        
!!$        IF(numvertx.NE.4)THEN
!!$           ALLOCATE(lm4node( 1:4,1:numelv_prmry))
!!$           ALLOCATE(imap(1:numnp_prmry))
!!$           imap(:) = 0
!!$           numnp4node = 0
!!$        ELSE
!!$           numnp4node = numnp_prmry
!!$        ENDIF
     ENDIF

     iEL10 = 0

!      ConvertUnit = .0254d0 ! for titan, units in engles
      
     DO i = 1, NumNP_glb(ibody)

! - Type 01: Node Data,  format(i2,8i8)
!        1 ID IV KC
!    ID = node id
!    IV = 0 n/a
!    KC = 2
        READ(iunit,*) j,id
        id = id + AccumNd
! - format(3e16.9)
! - caresian coordinate of Nodes x, y, z
        READ(iunit,'(3e16.9)') coor(1:3,id)
        coor(1:3,id) = coor(1:3,id)*ConvertUnit
! - ICF GTYPE NDF CONFIG CID PSPC
        READ(iunit,'()')

        iNodeTest = iNodeTest + 1
     ENDDO
     AccumNd_old = AccumNd
     AccumNd =  AccumNd + NumNP_glb(ibody)
!
! --  Read element connectivity array
     IF(ibody.EQ.1)THEN
        ALLOCATE(MatId(1:numelv_prmry))
        ALLOCATE(ElTypeId(1:numelv_prmry))
     ENDIF
! flages for element type id
!
!    4 = 4 node tetrahedral
!    8 = 8 node hexahedral
!   10 = 10 node tetrahedral
!   12 = 12 node wedge* 
!   15 = 15 node wedge*
!   20 = 20 node hexahdral*
!
!    * Not yet implemented

     MinEL = 1e9
     MaxEL = 0

     Element: DO i = 1, NumEL_glb(ibody)

! - Type 02: Element Data, format(i2,8i8)
!        2 ID IV KC N1
!     ID = element ID
!     IV = Shape (5 = tet)
        READ(iunit,*) j, id
        iElemTest = iElemTest + 1
        AccumEl_old = AccumEl
        id = AccumEl + id


!     NODES CONFIG CID CEID 
        READ(iunit,'(3i8)') numvertx,j, MatId(id)
        IF(numvertx.NE.4)THEN
           IF(numvertx.NE.10)THEN
              IF(numvertx.NE.8)THEN
                 PRINT*,'ERROR: FOUND UNSUPPORTED ELEMENT TYPE'
                 PRINT*,'Found element with ', numvertx,' nodes '
                 PRINT*,'Stopping'
                 STOP
              ENDIF
           ENDIF
        ENDIF

        ElTypeId(id) = numvertx
     
        IF(MatId(id).LE.0) MatId(id) = 1
!     LNODES
!        LNODES = Element corner nodes followed by additional nodes
        READ(iunit,*) lmelv_prmry(1:numvertx,id) 
        lmelv_prmry(1:numvertx,id) = lmelv_prmry(1:numvertx,id) + AccumNd_old
        
        iEL10 = iEL10 + 1
        IF(numvertx.EQ.10.AND.i.EQ.1.AND.nprocs.GT.1)THEN
           ALLOCATE(lm4node( 1:4,1:NumEL_glb(ibody)))
           IF(ibody.EQ.1) ALLOCATE(imap(1:numnp_prmry))
           imap(:) = 0
           PRINT*,'Re-setting imap'
        ENDIF
        IF(numvertx.EQ.4.AND.i.EQ.1.AND.nprocs.GT.1)THEN
           ALLOCATE(lm4node( 1:4,1:NumEL_glb(ibody)))
           IF(ibody.EQ.1) ALLOCATE(imap(1:numnp_prmry))
           imap(:) = 0
           PRINT*,'Re-setting imap'
        ENDIF
        IF(numvertx.EQ.10.AND.nprocs.GT.1)THEN
           DO j = 1, 4
              IF(imap(lmelv_prmry(j,id)).EQ.0)THEN
                 numnp4node = numnp4node + 1
                 imap(lmelv_prmry(j,id)) = numnp4node
              ENDIF
              iaux = iaux + 1
              lm4node(j,iEL10) = imap(lmelv_prmry(j,id))
           ENDDO

        ELSE IF(numvertx.EQ.4.AND.nprocs.GT.1)THEN
           DO j = 1, 4
              IF(imap(lmelv_prmry(j,id)).EQ.0)THEN
                 numnp4node = numnp4node + 1
                 imap(lmelv_prmry(j,id)) = numnp4node
              ENDIF
              iaux = iaux + 1
              lm4node(j,iEL10) = imap(lmelv_prmry(j,id))
           ENDDO
        ENDIF
        MaxEl = SUM(NumEL_glb(1:ibody)) ! SUM(          MAX(MaxEl,id)
        MinEl = AccumEl_old + 1 !MIN(,id)


!
! -- Find the size of the smallest element
!
        xx = coor(1,lmelv_prmry(1,id)) - coor(1,lmelv_prmry(2,id))
        yy = coor(2,lmelv_prmry(1,id)) - coor(2,lmelv_prmry(2,id))
        zz = coor(3,lmelv_prmry(1,id)) - coor(3,lmelv_prmry(2,id))
        size1 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = coor(1,lmelv_prmry(2,id)) - coor(1,lmelv_prmry(3,id))
        yy = coor(2,lmelv_prmry(2,id)) - coor(2,lmelv_prmry(3,id))
        zz = coor(3,lmelv_prmry(2,id)) - coor(3,lmelv_prmry(3,id))
        size2 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = coor(1,lmelv_prmry(3,id)) - coor(1,lmelv_prmry(1,id))
        yy = coor(2,lmelv_prmry(3,id)) - coor(2,lmelv_prmry(1,id))
        zz = coor(3,lmelv_prmry(3,id)) - coor(3,lmelv_prmry(1,id))
        size3 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(1,id))
        yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(1,id))
        zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(1,id))
        size4 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(2,id))
        yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(2,id))
        zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(2,id))
        size5 = SQRT(xx*xx+yy*yy+zz*zz)
        xx = coor(1,lmelv_prmry(4,id)) - coor(1,lmelv_prmry(3,id))
        yy = coor(2,lmelv_prmry(4,id)) - coor(2,lmelv_prmry(3,id))
        zz = coor(3,lmelv_prmry(4,id)) - coor(3,lmelv_prmry(3,id))
        size6 = SQRT(xx*xx+yy*yy+zz*zz)
        dhmin = MIN(size1,size2,size3,size4,size5,size6,dhmin)
     ENDDO Element

     PRINT*,'Done reading Element'
     AccumEl = AccumEl + NumEL_glb(ibody)

     IF(nprocs.GT.1)THEN
        IF(numvertx.EQ.10)THEN
           ALLOCATE(npart(1:numnp4node),stat=ierror)
           IF(ierror .NE.0)THEN
              ! Space for npart could not be allcoated
              PRINT*,'Program could not allocate space for npart'
              PRINT*,'Dimensions 1 to ', numnp4node
              STOP
           ENDIF
              
!!$           IF(ibody.EQ.4) THEN
!!$              WRITE(44,*) NumEL_glb(ibody),numnp4node,nprocs
!!$              MaxNd = MAXVAL(lm4node)
!!$              MinNd = MINVAL(lm4node)
!!$              PRINT*,'MaxNd =', MaxNd
!!$              PRINT*,'MinND =', MinNd
!!$              PRINT*,'MinNDlocation',MINLOC(lm4node)
!!$              PRINT*,'Maximum renumbered node =', numnp4node
!!$              PRINT*,'Number of elements in section',NumEL_glb(ibody)
!!$              PRINT*,'MinEl:MaxEl',MinEl,MaxEl
!!$              DO i = 1, NumEL_glb(ibody)
!!$                 WRITE(44,*) lm4node(:,i)
!!$                 IF(ANY(lm4node(:,i) .EQ. 0)) THEN
!!$                    PRINT*,'zero connectivity'
!!$                    STOP
!!$                 ENDIF
!!$              ENDDO
!!$              CLOSE(44)
!!$           ENDIF
!           PRINT*,'METIS for elements ',MinEL,':',MaxEl


!!$           CALL METIS_PartMeshDual(NumEL_glb(ibody),numnp4node, &
!!$                lm4node,2,1,nprocs,edgecut,epart(MinEl:MaxEl),npart)

           PRINT*,'METIS for elements ',MinEL,':',MaxEl,NumEL_glb(ibody)

           ALLOCATE(epart_loc(1:NumEL_glb(ibody)),stat=ierror)

            IF(ierror .NE.0)THEN
              ! Space for npart could not be allcoated
              PRINT*,'Program could not allocate space for epart_loc'
              PRINT*,'Dimensions 1 to ', NumEL_glb(ibody)
              STOP
           ENDIF

           PRINT*,'CALLING METIS'

           CALL METIS_PartMeshNodal(NumEL_glb(ibody),numnp4node, &
                lm4node,2,1,nprocs,edgecut,epart_loc,npart)

           epart(MinEl:MaxEl) = epart_loc(1:NumEL_glb(ibody))
           print*,MinEl,MaxEl,1,NumEL_glb(ibody)
           deallocate(epart_loc)

           PRINT*,'Finished Metis'
           DEALLOCATE(npart)
           DEALLOCATE(lm4node)

        ELSE IF(numvertx.EQ.4)THEN
           ALLOCATE(npart(1:NumNP_glb(ibody)),stat=ierror)
           IF(ierror .NE.0)THEN
              ! Space for npart could not be allcoated
              PRINT*,'Program could not allocate space for npart'
              PRINT*,'Dimensions 1 to ', numnp4node
              STOP
           ENDIF

!!$           CALL METIS_PartMeshDual(NumEL_glb(ibody),numnp4node, &
!!$                lm4node,2,1,nprocs,edgecut,epart(MinEl:MaxEl),npart)

           PRINT*,'METIS for elements ',MinEL,':',MaxEl,NumEL_glb(ibody)

           ALLOCATE(epart_loc(1:NumEL_glb(ibody)),stat=ierror)

            IF(ierror .NE.0)THEN
              ! Space for npart could not be allcoated
              PRINT*,'Program could not allocate space for epart_loc'
              PRINT*,'Dimensions 1 to ', NumEL_glb(ibody)
              STOP
           ENDIF

           PRINT*,'CALLING METIS',NumEL_glb(ibody),NumNP_glb(ibody)

           CALL METIS_PartMeshNodal(NumEL_glb(ibody),NumNP_glb(ibody), &
                lm4node,2,1,nprocs,edgecut,epart_loc,npart)

           epart(MinEl:MaxEl) = epart_loc(1:NumEL_glb(ibody))
           print*,MinEl,MaxEl,1,NumEL_glb(ibody)
           deallocate(epart_loc)

           PRINT*,'Finished Metis'
           DEALLOCATE(npart)
           DEALLOCATE(lm4node)


        ENDIF
     ELSE
        epart = 1
     ENDIF
  ENDDO

!!$
!!$  OPEN(45,file='LoadBalanceStats.out')
!!$  ALLOCATE(ProcDist(1:nprocs))
!!$  ProcDist(:) = 0
!!$  DO i = 1, numelv_prmry
!!$     IF(epart(i).EQ.0)THEN
!!$        print*,'epart = 0',i
!!$        stop
!!$     endif
!!$     ProcDist(epart(i)) = ProcDist(epart(i)) + 1
!!$  ENDDO
!!$
!!$  DO i=1, nprocs
!!$     write(45,*) i, ProcDist(i)
!!$  enddo
!!$  PRINT*,'*** Load Balancing Stats ***'
!!$  PRINT*,'Maximum number of elements on a processor =', MAXVAL(ProcDist(:))
!!$  PRINT*,'Minimum number of elements on a processor =', MINVAL(ProcDist(:))
!!$  close(45)
!!$  DEALLOCATE(ProcDist)



  IF(iNodeTest.NE.numnp_prmry)THEN
     PRINT*,'ERROR number of combined nodes not equal total nodes'
     PRINT*,'Sum of Nodes =', numnp_prmry
     PRINT*,'Actual Number Read =', iNodeTest
     STOP
  ENDIF

  IF(iElemTest.NE.numelv_prmry)THEN
     PRINT*,'ERROR number of combined elements not equal total elements'
     PRINT*,'Sum of Elements =', numelv_prmry
     PRINT*,'Actual Number Read =', iElemTest
     STOP
  ENDIF

! -- Partition the finite element mesh

  PRINT*,' REQUESTED NUMBER OF PARTITIONS =',nprocs

!
!     Call METIS using the partition the mesh
!


  PRINT*,'CALLING METIS'

  IF(nprocs.GT.1)THEN
!-ALT       CALL METIS_PartMeshNodal(numel_old_z,nn,
!-ALT$          elmnts,1,1,nprocs,edgecut,epart_p,npart_p)
     IF(numvertx.EQ.8)THEN
        ALLOCATE(npart(1:numnp_prmry))
        CALL METIS_PartMeshDual(numelv_prmry,numnp_prmry, &
             lmelv_prmry,3,1,nprocs,edgecut,epart,npart)
!        DEALLOCATE(npart)
     ENDIF
!!$     ELSE IF(numvertx.EQ.4)THEN
!!$        ALLOCATE(npart(1:numnp_prmry))
!!$        CALL METIS_PartMeshDual(numelv_prmry,numnp_prmry, &
!!$             lmelv_prmry,2,1,nprocs,edgecut,epart,npart)
!!$        DEALLOCATE(npart)
!!$     ENDIF
  ELSE
     epart(:) = 1
  ENDIF

  PRINT*,' Called METIS'
  OPEN(45,file='LoadBalanceStats.out')
  ALLOCATE(ProcDist(1:nprocs))
  ProcDist(:) = 0
  DO i = 1, numelv_prmry
     IF(epart(i).EQ.0)THEN
        print*,'epart = 0',i
        stop
     endif
     ProcDist(epart(i)) = ProcDist(epart(i)) + 1
  ENDDO

  DO i=1, nprocs
     write(45,*) i, ProcDist(i)
  enddo
  PRINT*,'*** Load Balancing Stats ***'
  PRINT*,'Maximum number of elements on a processor =', MAXVAL(ProcDist(:))
  PRINT*,'Minimum number of elements on a processor =', MINVAL(ProcDist(:))
  close(45)
  DEALLOCATE(ProcDist)


! create communication file ! Note: should be able to remove restriction

  ALLOCATE(ProcNdList(1:numnp_prmry,1:MaxNumberOfProcsToShareNode) )
  ALLOCATE(NumProcPerNd(1:numnp_prmry))
  ALLOCATE(NumElPerProc(1:nprocs))
  ALLOCATE(NumNdPerProc(1:nprocs))
  
  NumNdPerProc(:) = 0
  NumElPerProc(:) = 0

  PRINT*,'Calling NewCommList'
  CALL NewCommList(numelv_prmry, numnp_prmry, nprocs, NumVertx,  &
       lmelv_prmry, epart,NumProcPerNd,ProcNdList,MaxNumberOfProcsToShareNode,NumElPerProc, &
       NumNdPerProc)

  PRINT*,'  .. completed'

!
! Set up the numbering of nodes for the implicit solver
!

  IF ( IMP ) THEN

     ALLOCATE(StartNumNP_loc_implicit(1:nprocs))
     ALLOCATE(NumNP_loc_implicit(1:nprocs))
     ALLOCATE(MapNodeImp(1:numnp_prmry))
     ALLOCATE(NodeProcImpGlobal(1:numnp_prmry))
     NumNP_loc_implicit(:) = 0
 
     DO i = 1, numnp_prmry
        NumNP_loc_implicit(npart(i)) = NumNP_loc_implicit(npart(i)) + 1
        MapNodeImp(i) = NumNP_loc_implicit(npart(i))
        NodeProcImpGlobal(i) = npart(i) - 1
     END DO
     
     StartNumNP_Loc_implicit(1) = 0
     DO i = 2, nprocs
        StartNumNP_Loc_implicit(i) = SUM(NumNP_loc_implicit(1:i-1))
     ENDDO

     DO i = 1, numnp_prmry
        MapNodeImp(i) = MapNodeImp(i) + StartNumNP_Loc_implicit(npart(i))
     END DO

  END IF

  IF(nprocs.GT.1)THEN
     IF(numvertx.EQ.8)THEN
        DEALLOCATE(npart)
     END IF
  END IF


! -- create a view of the partitioned mesh
! -- PMVIS software
! -- Command: pmvis.sgi.bin -n pmvis.nod -c pmvis.ele -p pmvis.part -o 1

  iopmvis = 0
  
  IF(iopmvis.EQ.1)THEN
     
     OPEN(13,FILE='pmvis.nod')
     DO i = 1, numnp_prmry
        WRITE(13,'(3e16.9)') coor(1:3,i)
     END DO

     CLOSE(13)
     OPEN(13,FILE='pmvis.ele')
     DO i = 1, numelv_prmry
        WRITE(13,'(4i10)') lmelv_prmry(1:4,i)
     ENDDO
     CLOSE(13)
     
     OPEN(13,FILE='pmvis.part')
     DO i = 1, numelv_prmry
        WRITE(13,'(1i10)') epart(i)
     ENDDO
     CLOSE(13)
  ENDIF

!
!
! -- Continue to read in the rest of the input parameters

  numbc_prmry = 0
  numbc_prmry_mm = 0
  numbc_prmry_ht = 0
  NumSerBC = 0
  
!!$  ALLOCATE(numel_2d(1:nprocs,1:2))
!!$  numel_2d(:,:) = 0
  

  ALLOCATE(iNdsBurnFlg(1:numnp_prmry))
  
  iNdsBurnFlg(:) = 0

! Initialize list for boundary conditions

  NULLIFY(BC_structural_head,BC_structural_tail)
  NULLIFY(BC_meshmotion_head,BC_meshmotion_tail)
  NULLIFY(BC_thermal_head,BC_thermal_tail)
  
  NULLIFY(SurfMesh_tri3_Ov1_head,SurfMesh_tri3_Ov1_tail)
  NULLIFY(SurfMesh_tri3_Ov2_head,SurfMesh_tri3_Ov2_tail)
  NULLIFY(SurfMesh_tri3_S_head,SurfMesh_tri3_S_tail)
  NULLIFY(SurfMesh_tri3_SF_head,SurfMesh_tri3_SF_tail)
  NULLIFY(SurfMesh_tri6_S_head,SurfMesh_tri6_S_tail)
  NULLIFY(SurfMesh_tri6_SF_head,SurfMesh_tri6_SF_tail)
  NULLIFY(SurfMesh_hex8_S_head,SurfMesh_hex8_S_tail)
  NULLIFY(SurfMesh_hex8_SF_head,SurfMesh_hex8_SF_tail)

! number of boundary conditions per processors

  ALLOCATE(NumBC_structural(1:nprocs))
  ALLOCATE(NumBC_meshmotion(1:nprocs))
  ALLOCATE(NumBC_thermal(1:nprocs))
  ALLOCATE(NumBC_Flag(1:nprocs))
  NumBC_structural(:) = 0
  NumBC_meshmotion(:) = 0
  NumBC_thermal(:)    = 0 
  NumBC_Flag(:)       = 0
  ALLOCATE(BC_Flag(1:3,1:numnp_prmry))
 
  BC_Flag(:,:) = 0
  

  PRINT*,'FINISH READING PATRAN FILE'
  PRINT*,'Reading Boundary Conditions'
  AccumEl = 0
  AccumNd = 0
  DO ibody = 1, NumGeomBodies
     iunit = 100 + ibody


     DO
        READ(IUNIT,*) itype,id,iv,kc
        IF(itype.EQ.99) EXIT
        
        IF(itype.EQ.4)THEN
           READ(IUNIT,'()')

! Marks of IntFaceFlag :
!
!   SolidFluid Interface = 1
!   NotASolidFluid Interface = 0
         
        ELSE IF(itype.EQ. 3) THEN
           DO i = 1, 20
              READ(IUNIT,'()')
           ENDDO
        ELSE IF(itype .EQ. 6)THEN      ! pressure loading
           id = id + AccumEl
           ip = epart(id)
           READ(IUNIT,'(i1,i1,i1,i6,8i1)')i,i,i,i,n1,n2,n3,n4,n5,n6,n7,n8
           READ(IUNIT,*) press
           IntFaceFlag = ABS(INT(press))

           IF(.NOT.(InteractMesh))THEN
              IF(IntFaceFlag.GE.0.AND.IntFaceFlag.LE.2)  InteractMesh = .TRUE.
           ENDIF

! tet elements
           
           IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.1   &
                .AND.n4.EQ.0.AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.4)THEN
                 nnd1 = lmelv_prmry(1,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(2,id)
                 !itmp = itmp + 1
                 iFaceID = 1
              ELSE
                 nnd1 = lmelv_prmry(1,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(2,id)
                 nnd4 = lmelv_prmry(7,id)
                 nnd5 = lmelv_prmry(6,id)
                 nnd6 = lmelv_prmry(5,id)
                 iFaceID = 1
              ENDIF
           ELSE IF(n1.EQ.1.AND.n2.EQ.1.AND.n4.EQ.1  &
                .AND.n3.EQ.0.AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.4)THEN
                 nnd1 = lmelv_prmry(1,id)
                 nnd2 = lmelv_prmry(2,id)
                 nnd3 = lmelv_prmry(4,id)
                 !itmp = itmp + 1
                 iFaceID = 2
              ELSE
                 nnd1 = lmelv_prmry(1,id)
                 nnd2 = lmelv_prmry(2,id)
                 nnd3 = lmelv_prmry(4,id)
                 nnd4 = lmelv_prmry(5,id)
                 nnd5 = lmelv_prmry(9,id)
                 nnd6 = lmelv_prmry(8,id)
                 iFaceID = 2
              ENDIF
           ELSE IF(n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1  &
                .AND.n1.EQ.0.AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.4)THEN
                 nnd1 = lmelv_prmry(2,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(4,id)
                 !itmp = itmp + 1
                 iFaceID = 3
              ELSE
                 nnd1 = lmelv_prmry(2,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(4,id)
                 nnd4 = lmelv_prmry(6,id)
                 nnd5 = lmelv_prmry(10,id)
                 nnd6  = lmelv_prmry(9,id)
                 iFaceID = 3
              ENDIF
           ELSE IF(n1.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1  &
                .AND.n2.EQ.0.AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.4)THEN
                 nnd1 = lmelv_prmry(4,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(1,id)
                 itmp = itmp + 1
                 iFaceID = 4
              ELSE
                 nnd1 = lmelv_prmry(4,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(1,id)
                 nnd4 = lmelv_prmry(10,id)
                 nnd5 = lmelv_prmry(7,id)
                 nnd6 = lmelv_prmry(8,id)
                 iFaceID = 4
              ENDIF
! Hex elements
           ELSE IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.1  &
                .AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(4,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(2,id)
                 nnd4 = lmelv_prmry(1,id)
                 itmp = itmp + 1
                 iFaceID = 1
              ENDIF
           ELSE IF(n1.EQ.0.AND.n2.EQ.0.AND.n3.EQ.0.AND.n4.EQ.0  &
                .AND.n5.EQ.1.AND.n6.EQ.1.AND.n7.EQ.1.AND.n8.EQ.1)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(5,id)
                 nnd2 = lmelv_prmry(6,id)
                 nnd3 = lmelv_prmry(7,id)
                 nnd4 = lmelv_prmry(8,id)
                 itmp = itmp + 1
                 iFaceID = 2
              ENDIF
           ELSE IF(n1.EQ.1.AND.n2.EQ.1.AND.n3.EQ.0.AND.n4.EQ.0  &
                .AND.n5.EQ.1.AND.n6.EQ.1.AND.n7.EQ.0.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(1,id)
                 nnd2 = lmelv_prmry(2,id)
                 nnd3 = lmelv_prmry(6,id)
                 nnd4 = lmelv_prmry(5,id)
                 itmp = itmp + 1
                 iFaceID = 3
              ENDIF
           ELSE IF(n1.EQ.0.AND.n2.EQ.1.AND.n3.EQ.1.AND.n4.EQ.0  &
                .AND.n5.EQ.0.AND.n6.EQ.1.AND.n7.EQ.1.AND.n8.EQ.0)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(2,id)
                 nnd2 = lmelv_prmry(3,id)
                 nnd3 = lmelv_prmry(7,id)
                 nnd4 = lmelv_prmry(6,id)
                 itmp = itmp + 1
                 iFaceID = 4
              ENDIF
           ELSE IF(n1.EQ.0.AND.n2.EQ.0.AND.n3.EQ.1.AND.n4.EQ.1  &
                .AND.n5.EQ.0.AND.n6.EQ.0.AND.n7.EQ.1.AND.n8.EQ.1)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(3,id)
                 nnd2 = lmelv_prmry(4,id)
                 nnd3 = lmelv_prmry(8,id)
                 nnd4 = lmelv_prmry(7,id)
                 itmp = itmp + 1
                 iFaceID = 5
              ENDIF
           ELSE IF(n1.EQ.1.AND.n2.EQ.0.AND.n3.EQ.0.AND.n4.EQ.1  &
                .AND.n5.EQ.1.AND.n6.EQ.0.AND.n7.EQ.0.AND.n8.EQ.1)THEN
              IF(ElTypeId(id).EQ.8)THEN
                 nnd1 = lmelv_prmry(5,id)
                 nnd2 = lmelv_prmry(8,id)
                 nnd3 = lmelv_prmry(4,id)
                 nnd4 = lmelv_prmry(1,id)
                 itmp = itmp + 1
                 iFaceID = 6
              ENDIF
           ELSE
              PRINT*,'Error in pressure face numbering'
              STOP
           ENDIF
        
! BC flags
! ignitable (=1)
! non-ignitable (=0)
! noninteracting (=2)


! 3-node triangle, (via 4-node tetrahedral)
        
           IF(ElTypeId(id).EQ.4)THEN

              ALLOCATE(SurfMesh_tri3_item)
              
              SurfMesh_tri3_item%ElemData(1) = nnd1
              SurfMesh_tri3_item%ElemData(2) = nnd2
              SurfMesh_tri3_item%ElemData(3) = nnd3
              SurfMesh_tri3_item%ElemData(4) = id
              
              IF(IntFaceFlag.EQ.1)THEN ! Fluid Structure interaction, with burning
                 CALL add_SurfMesh_tri3(SurfMesh_tri3_item, SurfMesh_tri3_SF_head, SurfMesh_tri3_SF_tail)
                 NumEltet2D(1,ip) = NumEltet2D(1,ip) + 1
                 
              ELSE IF(IntFaceFlag.EQ.2)THEN ! fluid Structure interaction, without burning
                 CALL add_SurfMesh_tri3(SurfMesh_tri3_item, SurfMesh_tri3_SF_NonIgnt_head, SurfMesh_tri3_SF_NonIgnt_tail)
                 NumEltet2D(3,ip) = NumEltet2D(3,ip) + 1
              ELSE IF(IntFaceFlag.EQ.0)THEN ! No fluid structure interaction
                 CALL add_SurfMesh_tri3(SurfMesh_tri3_item, SurfMesh_tri3_S_head, SurfMesh_tri3_S_tail)
                 NumEltet2D(2,ip) = NumEltet2D(2,ip) + 1
              ELSE IF(IntFaceFlag.EQ.10)THEN ! No fluid structure interaction for non-matching meshes
                 OverlayMesh = .TRUE.
                 SurfMesh_tri3_item%ElemData(5) =  iFaceID
                 CALL add_SurfMesh_tri3(SurfMesh_tri3_item, SurfMesh_tri3_Ov1_head, SurfMesh_tri3_Ov1_tail)
                 NumEltet2D(4,ip) = NumEltet2D(4,ip) + 1
              ELSE IF(IntFaceFlag.EQ.100)THEN ! No fluid structure interaction for non-matching meshes    
                 OverlayMesh = .TRUE.    
                 SurfMesh_tri3_item%ElemData(5) =  iFaceID         
                 CALL add_SurfMesh_tri3(SurfMesh_tri3_item, SurfMesh_tri3_Ov2_head, SurfMesh_tri3_Ov2_tail)
                 NumEltet2D(5,ip) = NumEltet2D(5,ip) + 1

              ENDIF
              Tet4 = .TRUE.

! 6-node triangle, (via 10-node tetrahedral)

           ELSE IF(ElTypeId(id).EQ.10)THEN
              
              ALLOCATE(SurfMesh_tri6_item)
              
              SurfMesh_tri6_item%ElemData(1) = nnd1
              SurfMesh_tri6_item%ElemData(2) = nnd2
              SurfMesh_tri6_item%ElemData(3) = nnd3
              SurfMesh_tri6_item%ElemData(4) = nnd4
              SurfMesh_tri6_item%ElemData(5) = nnd5
              SurfMesh_tri6_item%ElemData(6) = nnd6
              SurfMesh_tri6_item%ElemData(7) = id
              
              IF(IntFaceFlag.EQ.1)THEN
                 CALL add_SurfMesh_tri6(SurfMesh_tri6_item, SurfMesh_tri6_SF_head, SurfMesh_tri6_SF_tail)
                 NumEltet2D(1,ip) = NumEltet2D(1,ip) + 1
              ELSEIF(IntFaceFlag.EQ.2)THEN
                 CALL add_SurfMesh_tri6(SurfMesh_tri6_item, SurfMesh_tri6_SF_NonIgnt_head, SurfMesh_tri6_SF_NonIgnt_tail)
                 NumEltet2D(3,ip) = NumEltet2D(3,ip) + 1
              ELSE IF(IntFaceFlag.EQ.0) THEN
                 CALL add_SurfMesh_tri6(SurfMesh_tri6_item, SurfMesh_tri6_S_head, SurfMesh_tri6_S_tail)
                 NumEltet2D(2,ip) = NumEltet2D(2,ip) + 1
              ENDIF
              tet10 = .TRUE.

! 4-node quad, (via 8-node hexahedral)
          
           ELSE IF(ElTypeId(id).EQ.8)THEN
              
              ALLOCATE(SurfMesh_hex8_item)
              
              SurfMesh_hex8_item%ElemData(1) = nnd1
              SurfMesh_hex8_item%ElemData(2) = nnd2
              SurfMesh_hex8_item%ElemData(3) = nnd3
              SurfMesh_hex8_item%ElemData(4) = nnd4
              SurfMesh_hex8_item%ElemData(5) = id
              IF(IntFaceFlag.EQ.1)THEN
                 CALL add_SurfMesh_hex8(SurfMesh_hex8_item, SurfMesh_hex8_SF_head, SurfMesh_hex8_SF_tail)
                 NumElhex2D(1,ip) = NumElhex2D(1,ip) + 1
              ELSE IF(IntFaceFlag.EQ.2)THEN
                 CALL add_SurfMesh_hex8(SurfMesh_hex8_item, SurfMesh_hex8_SF_NonIgnt_head, SurfMesh_hex8_SF_NonIgnt_tail)
                 NumElhex2D(3,ip) = NumElhex2D(3,ip) + 1
              ELSE IF(IntFaceFlag.EQ.0) THEN
                 CALL add_SurfMesh_hex8(SurfMesh_hex8_item, SurfMesh_hex8_S_head, SurfMesh_hex8_S_tail)
                 NumElhex2D(2,ip) = NumElhex2D(2,ip) + 1
              ELSE
                 PRINT*,'ERROR, interface flag type not found',IntFaceFlag
                 stop
              ENDIF
              
              hex8 = .TRUE.
           ENDIF
        

!        numel_2d(ip,IntFaceFlag) = numel_2d(ip,IntFaceFlag) + 1
        

        ELSE IF(itype .EQ. 8)THEN ! boundary conditions
           id = id + AccumNd
! 
!     CID = Coordinate frame ID
!     ICOMP = 6 displacement component flags (0 or 1)
!     note the flag is passed to node number one (i.e. fixed displacement)

           READ(IUNIT,'(I8,6I1)') iaux, DisFlag(1:6)
           NumDisFlag = SUM(DisFlag(1:6))
           

           NumSerBC = NumSerBC + NumDisFlag

           READ(IUNIT,*) tmpVal(1:NumDisFlag)
           
           icnt = 0
           do i = 1, 6
              IF(DisFlag(i).eq.1)THEN
                 icnt = icnt + 1
                 DisFlagValue(i) = tmpVal(icnt)
              endif
           enddo

           IF(icnt.NE.NumDisFlag)THEN
              PRINT*,'ERROR in processing boundary conditions'
              STOP
           ENDIF
        
! add to the number of boundary conditions per processor  

           IF(DisFlag(1).eq.1 )THEN ! structural boundary condition
              DO i = 1,NumProcPerNd(id)
                 ProcId = ProcNdList(id,i)
                 NumBC_structural(ProcId) = NumBC_structural(ProcId) + 1
                 IF(BC_Flag(1,id).EQ.0) NumBC_Flag(ProcId) = NumBC_Flag(ProcId) + 1 
              ENDDO

              BC_Flag(1,id) = BC_Flag(1,id) + INT(DisFlagValue(1))
              
              MaxNumBC_str = MAX(MaxNumBC_str,INT(DisFlagValue(1)))
           endif

           IF(DisFlag(2).eq.1 )THEN ! temperature boundary condition
              DO i = 1,NumProcPerNd(id)
                 ProcId = ProcNdList(id,i)
                 NumBC_thermal(ProcId) = NumBC_thermal(ProcId) + 1
                 IF(BC_Flag(2,id).EQ.0) NumBC_Flag(ProcId) = NumBC_Flag(ProcId) + 1 
              enddo
              
              BC_Flag(2,id) = BC_Flag(2,id) + 100*INT(DisFlagValue(2))
              
              MaxNumBC_th = MAX(MaxNumBC_th,INT(DisFlagValue(2)))
           endif


           IF(DisFlag(3).eq.1 )THEN ! mesh motion boundary condition
              DO i = 1,NumProcPerNd(id)
                 ProcId = ProcNdList(id,i)
                 NumBC_meshmotion(ProcId) = NumBC_meshmotion(ProcId) + 1
                 IF(BC_Flag(3,id).EQ.0) NumBC_Flag(ProcId) = NumBC_Flag(ProcId) + 1 
              ENDDO
              BC_Flag(3,id) = BC_Flag(3,id) + 10000*INT(DisFlagValue(3))
              
              MaxNumBC_mm = MAX(MaxNumBC_mm,INT(DisFlagValue(3)))
           endif


        
        ELSE IF(itype .EQ. 7)THEN ! mesh motion boundary conditions
           id = id + AccumNd
! 
!     CID = Coordinate frame ID
!     ICOMP = 6 displacement component flags (0 or 1)
!     note the flag is passed to node number one (i.e. fixed displacement)

           READ(IUNIT,'()')
           READ(IUNIT,*) value
           
!           ALLOCATE(BC_meshmotion_item)

!           BC_meshmotion_item%BC_nodeGlb = id
!           BC_meshmotion_item%BC_flagGlb = INT(value)
           
!           CALL add_BC(BC_meshmotion_item, BC_meshmotion_head, BC_meshmotion_tail)
           

! add to the number of boundary conditions per processor         

           DO i = 1,NumProcPerNd(id)
              ProcId = ProcNdList(id,i)
              NumBC_meshmotion(ProcId) = NumBC_meshmotion(ProcId) + 1
              IF(BC_Flag(3,id).EQ.0) NumBC_Flag(ProcId) = NumBC_Flag(ProcId) + 1 
           ENDDO

           BC_Flag(3,id) = BC_Flag(3,id) + 10000*INT(value)

           MaxNumBC_mm = MAX(MaxNumBC_mm,INT(value))

        ELSE IF(itype .EQ. 5)THEN ! coordinate frames
           READ(IUNIT,'()')
           READ(IUNIT,'()')
           READ(IUNIT,'()')
           READ(IUNIT,'()')
        ELSE IF(itype .EQ. 10)THEN ! nodal temperature (corners)
           id = id + AccumNd
           READ(IUNIT,*) value
           
!           ALLOCATE(BC_thermal_item)
           
!           BC_thermal_item%BC_nodeGlb = id
!           BC_thermal_item%BC_flagGlb = INT(value)
           
!           CALL add_BC(BC_thermal_item, BC_thermal_head, BC_thermal_tail)

! add to the number of boundary conditions per processor
  
           DO i = 1,NumProcPerNd(id)
              ProcId = ProcNdList(id,i)
              NumBC_thermal(ProcId) = NumBC_thermal(ProcId) + 1
              IF(BC_Flag(2,id).EQ.0) NumBC_Flag(ProcId) = NumBC_Flag(ProcId) + 1 
           ENDDO

           BC_Flag(2,id) = BC_Flag(2,id) + 100*INT(value)
           
           MaxNumBC_th = MAX(MaxNumBC_th,INT(value))

        ELSE
           DO i = 1, kc
              READ(iunit,'()')
           ENDDO
        ENDIF
     ENDDO

     AccumEl =  SUM(NumEL_glb(1:ibody))

     AccumNd =  SUM(NumNP_glb(1:ibody))
  ENDDO

! if negative number then higher order element


  IF(tet4.AND..NOT.(tet10).AND..NOT.(hex8))THEN
     MeshType2D = 3
  ELSE IF (.NOT.(tet4).AND.tet10.AND..NOT.(hex8))THEN
     MeshType2D = 6
  ELSE IF (.NOT.(tet4).AND..NOT.(tet10).AND.hex8)THEN
     MeshType2D = 4
  ELSE
     MeshType2D = 5
  ENDIF

!  PRINT*,'Number of Surface elements, Interacting with burning =', NumEltet2D(1,:)
!  PRINT*,'Number of Surface elements, Interacting with Non-burning =', NumEltet2D(3,:)
!  PRINT*,'Number of Surface elements, Non-Interacting =', NumEltet2D(2,:)

  PRINT*,'BOUNDARY CONDITIONS'
  PRINT*,' Total of all types = ',NumSerBC
  print*,'  Number of structural = ',MaxNumBC_str
  print*,'  Number of meshmotion = ',MaxNumBC_mm
  print*,'  Number of thermal = ',MaxNumBC_th

  IF(numvertx.EQ.4)THEN
     numvertx2d = 3
  ELSE IF(numvertx.EQ.10)THEN
     numvertx2d = 6
  ELSE IF(numvertx.EQ.8)THEN
     numvertx2d = 4
  ENDIF

  DO ibody = 1, NumGeomBodies
     iunit = 100 +  ibody 
     CLOSE(iunit)
  enddo
  
  RETURN
  
END SUBROUTINE read_patran
    

