      program hdf_check
!
! Checks HDF files written for Rocketeer v1.0 and v1.1.  
! Reports on array dimensions, ranges, etc.
!
! Prompts for an HDF file name.  If a blank is input, reads file 
! named test.hdf
!
! Compile on CSAR Suns using:
!
! % f90 reader.f -L/projects1/jnorris/lib -ldf -lz -ljpeg
!
! Compile on NCSA's Origin 2000 using:
!
! % f90 -n32 reader.f -L/afs/ncsa/packages/hdf/4.1r4_irix64-n32/lib -ldf -lz -ljpeg
!
! Compile on linux with the Absoft compiler using:
!
! % f90 reader.f -L/home/rfiedler/HDF/Absoft -ldf -lz -ljpeg
! 
! Compile on frost.llnl.gov using:
!
! xlf90_r -qfixed reader_hep.f -L/g/g20/mdbrandy/HDF/lib -ldf -lz -ljpeg -o reader_hep.x
!
! note: replace -L with path to your copy of HDF4 libraries.
!
! Written by Robert A. Fiedler, 4/21/2001
!
      implicit none
      integer ne,nn                ! Numbers of elements and nodes
      integer n1,n2,n3             ! structured mesh dimensions
!
! Arrays to store data sets, both single and double precision
!
      real, dimension (:), allocatable :: data
      real (kind=8), dimension (:), allocatable :: ddata
      integer, dimension (:,:), allocatable :: conn
      real, dimension (:), allocatable :: x1,x2,x3
      real (kind=8), dimension(:), allocatable :: dx1,dx2,dx3
      real, dimension (:,:,:), allocatable :: c1,c2,c3
      real (kind=8), dimension(:,:,:), allocatable :: dc1,dc2,dc3
      real, dimension (:,:,:), allocatable :: data3d
      real (kind=8), dimension (:,:,:), allocatable :: ddata3d

      integer ret, rank, maxrank
      parameter(maxrank=3)
      integer sizes(maxrank)
      integer num_type
      integer l_label,l_unit,l_format,l_coord,l_mat,l_block
      real version
      real fmin,fmax
      real (kind=8) dfmin,dfmax
      integer idata(2)
      integer l,m

      integer nchar
      parameter(nchar=64)
      character cdata(nchar)
      character*32 filename,filegeom
      character*64 string
      character*64 coordsys
      character*64 blockname
      character*64 citime, material
      character*64 label, units, format
      character*64 cheader
      character*1  bk
      logical skip_rank
      integer it
      integer expi
!
! HDF library routines
!
      integer  dsgdims, dsgdata, dsgrang, dsgnt, dsgdast, dsgdaln
      external dsgdims, dsgdata, dsgrang, dsgnt, dsgdast, dsgdaln

      bk = " "

!......................................................................
!
! Save head end pressure to a file for block_0001
!
      open(8,file='hep.txt',position='append')
!
! Get HDF file name
!
      filename = bk
      print *,'Enter file name: '
      read *,filename
      if (filename .eq. bk) filename = 'test.hdf'  ! Default name
      print *,bk
      print *,'Reading ',filename

!-------------------------------+
!                               |
!     Loop over data blocks     |
!                               |
!-------------------------------+
!
      do m=1,1000000
!
! Look for block headers.
!
! A block header consists of a mesh description array and a set of
! attributes.
!
! For version 1.0, the mesh description array has two integer values:
!
!    1) grid type:
!       0 -- uniform
!       1 -- rectilinear
!       2 -- structured
!       3 -- unstructured
!
!    2) number of layers of ghost cells (ignored for unstructured)
!
! For version 1.1, the mesh is described by a character string:
!
!    "grid_type|number_of_ghost_cells|file_with_mesh_data".
!
! For example,
!
!    cdata = '3' // '|' // '0' // '|' // 'unstr0.hdf' // char(0)
!
! ..indicates an unstructured mesh, no ghost zones, coordinate
! arrays and connectivity data located in file 'unstr0.hdf'.
! Note the terminating null character "char(0)".  The HDF
! library treats this as an array of 15 characters.  The
! file name is optional; by default, the coordinates are assumed
! to be in the current data file.
!
! The attributes of the mesh description array are set via:
!
!   call DFSDsetdatastrs(label, units, format, coord_sys)
!
!   block name in "label" attribute, e.g. 'block_001'
!   time level in "units" attribute, e.g. '5', '2.5ms' (v1.1)
!   the string 'block header' in "format" attribute
!   material type in "coord_sys" attribute
!
        print *,bk
        print *,'Looking for block header'
        if (m .eq. 1) then
!
! For the first block only; successive blocks will
! already have been detected, so skip this.
!
          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read header'
            stop
          endif

          if (rank .eq. 1) then
            print *,'rank of mesh type array is ',rank
          else
            print *,'mesh type array should be of rank 1, was ',rank
            stop
          endif

! FOR EXTRACTION OF HEAD-END PRESSURE
        else
          stop
        endif
!
! Do for all blocks
!
        skip_rank = .false.
!
!   Check lengths of attribute strings; there's no limit in
!   Rocketeer, but this program is limited by the value of nchar.
!
        ret = dsgdaln(l_block,l_unit,l_format,l_mat)
        if (ret .ne. 0) then
          print *,'Could not get header string lengths'
          stop
        endif

        if (max(l_block,l_unit,l_format,l_mat) >= nchar) then
          print *,'One of the header string labels is > ',nchar,' chars'
          stop
        endif

        ret = dsgdast(blockname,citime,cheader,material)
        if (ret .ne. 0) then
          print *,'Failed to read block header strings'
          stop
        endif
!
! Check the header attributes
!
        if (l_block .le. 1 .or. blockname(1:1) .eq. char(0)) then
          print *,'Invalid blockname: ',blockname,l_block,' characters'
          stop
        else
          print *,'blockname is ',blockname(1:l_block)
        endif

        if (l_unit .le. 1 .or. citime(1:1) .eq. char(0)) then
          print *,'Invalid time level: ',citime,l_unit,' characters'
          stop
        else
          print *,'time level is ',citime(1:l_unit)
        endif

        if (cheader(1:l_format) .ne. 'block header' .and.
     1      cheader(1:l_format) .ne. 'block header' // char(0)) then
          print *,'this string should be block header: ',  
     1             cheader(1:l_format)
          stop
        endif

        if (l_mat .le. 1 .or.  material(1:1) .eq. char(0)) then
          print *,'Invalid material type'
          stop
        else
          print *,'material type is ',material(1:l_mat)
        endif
!
! Check the mesh type array
!
        print *,bk
        print *,'Reading mesh type data'
        ret = dsgnt(num_type)
        if (ret .ne. 0) then
          print *,'Failed to get data type'
          stop
        endif

        if (num_type .eq. 24) then
!
! Valid for version 1.0 and later -- 2 integers.
!
          print *,'mesh type array is 32 bit signed integers'
          print *,'valid for Rocketeer version 1.0 and later'
          version = 1.00
          ret = dsgdata(filename,rank,sizes,idata)
        else if (num_type .eq. 4) then
!
! Valid for version 1.1 -- characters.
!
          print *,'mesh type array is character'
          print *,'valid for Rocketeer version 1.1 and later'
          version = 1.1
          ret = dsgdata(filename,rank,sizes,cdata)
        else
          print *,'mesh type array has wrong data type: ',num_type
          stop
        endif

        if (ret .ne. 0) then
          print *,'Failed to read mesh type data'
          stop
        endif

        if (num_type .eq. 24) then
!
! Version 1.0.
!
          if (sizes(1) .ne. 2) then
            print *,'Mesh type array should have 2 elements, has ',
     1              sizes(1)
            stop
          endif
          print *,'mesh type is ',idata(1)
          print *,'number of ghost zones is ',idata(2)
          print *,bk
          filegeom = filename  ! Geometrical data always in same file
        else
!
! Version 1.1
!
          if (sizes(1).lt.3) then
            print *,'Mesh type string should be of the form'
            print *,'type|ghost cells|filename (filename optional)'
            print *,'You need to include a | character,' 
            print *,'the number of ghost zones (0 for unstructured),'
            print *,'and terminate with a null character char(0)'
            stop
          endif

          print *,'mesh type is ',cdata(1)

          if (cdata(2) .ne. '|') then
            print *,'A | character is needed in the mesh type string'
            stop
          endif

          print *,'number of ghost zone layers is ',cdata(3)

          if (sizes(1).gt.4) then
            print *,'file containing geometry is ',cdata(5:sizes(1))

            if (cdata(sizes(1)) .ne. char(0)) then
              print *,'Geometry file name must end with null char'
              stop
            endif
!
! Write character array cdata to string filegeom
!
!            if (cdata(1)(sizes(1):sizes(1)) .eq. char(0)) then
!              filegeom = cdata(1)(5:sizes(1)-1)
!            else
!              filegeom = cdata(1)(5:sizes(1))
!            endif
!
          else
!
! This file has the geometry data
!
            filegeom = filename
          endif  ! length of mesh type string
!
! Check mesh type
!
          print *,bk
          if (cdata(1).eq.'0') then
            idata(1) = 0
          else if (cdata(1).eq.'1') then
            idata(1) = 1
          else if (cdata(1).eq.'2') then
            idata(1) = 2
          else if (cdata(1).eq.'3') then
            idata(1) = 3
          else if (cdata(1).eq.'4') then
            idata(1) = 3  ! Treat discontinuous like unstructured
            print *,'Discontinous Galerkin FE mesh, treated like'
          else if (cdata(1).eq.'5') then
            idata(1) = 5
          else
            print *,'Invalid mesh type: ',idata(1)
            stop
          endif
!
! Check number of ghost cell layers
!
          if (cdata(3).eq.'0') then
            idata(2) = 0
          else if (cdata(3).eq.'1') then
            idata(2) = 1
          else if (cdata(3).eq.'2') then
            idata(2) = 2
          else if (cdata(3).eq.'3') then
            idata(2) = 3
          else if (cdata(3).eq.'4') then
            idata(2) = 4
          else
            print *,'Too many ghost cell layers'
            stop
          endif

        endif  ! integers or characters in mesh type
!
        if (idata(1).eq.0) then
!
!------------------------------+
!                              |
!         Uniform mesh         |
!                              |
!------------------------------+
!
! No coordinate arrays to read

          print *,'Uniform mesh'
          print *,bk

! Prepare to read in field variables

          print *,'Reading field variables'
          print *,bk
          do l=1,100

            ret = dsgdims(filename,rank,sizes,maxrank)
            if (ret .ne. 0) then
!              print *,'Failed to read field variable dimensions'
              stop
            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (rank.eq.1 .and. (num_type.eq.24 .or. num_type.eq.4))   
     1        exit  ! This could be a block header
            if (rank .ne. 3) then
              print *,'Field variable arrays must be of rank 3, read ',
     1                rank
              stop
            endif

            print *,'Dimensions are ',sizes(1),sizes(2),sizes(3)

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get attribute string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'One of the attribute string labels is > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(string,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read field variable attributes'
              stop
            endif

            print *,'label is ',string(1:l_label)
            print *,'units are ',units(1:l_unit)
            print *,'short name is ',format(1:l_format)
            print *,'component is ',coordsys(1:l_coord)
            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range'
              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is from ',fmin,' to ',fmax
              allocate(data3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,data3d)
              deallocate(data3d)
            else
              print *,'Range is from ',dfmin,' to ',dfmax
              allocate(ddata3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,ddata3d)
              deallocate(ddata3d)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read field variable'
              stop
            endif

            print *,'Field variable read in'
            print *,bk

          end do  ! loop on field variables

        else if (idata(1).eq.1) then
!
!---------------------------------+
!                                 |
!         Rectilinear mesh        |
!                                 |
!---------------------------------+
!
! Read nodal coordinates as 3 separate 1-D arrays.

          print *,'Rectilinear mesh'
          print *,'Reading coordinate arrays'
          print *,bk
!
! X
!
          print *,'Reading X coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 1st coord array'
            stop
          endif

          if (rank .ne. 1) then
            print *,'Coordinates must be stored as 3 separate ',
     1              '1-D arrays'
            stop
          endif

          print *,'Size of 1st dimension is ',sizes(1)

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     1              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 1st coord array attributes'
            stop
          endif

          print *,'Label of 1st coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)
          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 1st coord array'
c            stop
          endif

          n1 = sizes(1)
          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(x1(n1))
            ret = dsgdata(filename,rank,sizes,x1)
            deallocate(x1)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dx1(n1))
            ret = dsgdata(filename,rank,sizes,dx1)
            deallocate(dx1)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 1st coord array'
            stop
          endif

          print *,'X coordinates read in'
          print *,bk 
!
! Y
!
          print *,'Reading Y coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 2nd coord array'
            stop
          endif

          if (rank .ne. 1) then
            print *,'Coordinates must be stored as 3 separate ',
     1              '1-D arrays'
            stop
          endif

          n2 = sizes(1)
          print *,'Size of 2nd dimension is ',n2

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     1              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 2nd coord array attributes'
            stop
          endif

          print *,'Label of 2nd coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)
          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 2nd coord array'
c            stop
          endif

          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(x2(n2))
            ret = dsgdata(filename,rank,sizes,x2)
            deallocate(x2)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dx2(n2))
            ret = dsgdata(filename,rank,sizes,dx2)
            deallocate(dx2)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 2nd coord array'
            stop
          endif

          print *,'Y coordinates read in'
          print *,bk 
!
! Z
!
          print *,'Reading Z coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 3rd coord array'
            stop
          endif

          if (rank .ne. 1) then
            print *,'Coordinates must be stored as 3 separate ',
     1              '1-D arrays'
            stop
          endif

          n3 = sizes(1)
          print *,'Size of 3rd dimension is ',n3

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     !              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 3rd coord array attributes'
            stop
          endif

          print *,'Label of 3rd coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)

          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 3rd coord array'
c            stop
          endif

          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(x3(n3))
            ret = dsgdata(filename,rank,sizes,x3)
            deallocate(x3)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dx3(n3))
            ret = dsgdata(filename,rank,sizes,dx3)
            deallocate(dx3)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 3rd coord array'
            stop
          endif

          print *,'Z coordinates read in'
          print *,bk 

! Prepare to read in field variables

          print *,'Reading field variables'
          print *,bk
          do l=1,100

            ret = dsgdims(filename,rank,sizes,maxrank)
            if (ret .ne. 0) then
!              print *,'Failed to read field variable dimensions'
              stop
            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (rank.eq.1 .and. (num_type.eq.24 .or. num_type.eq.4))   
     1        exit  ! This could be a block header

            if (rank .ne. 3) then
              print *,'Field variable arrays must be of rank 3, read ',
     1                rank
              stop
            endif

            if (sizes(1)-2*idata(2).eq.n1-1) then
              print *,'Element centered data'
            else if (sizes(1)-2*idata(2).eq.n1) then
              print *,'Node centered data'
            else if (sizes(1).eq.n1) then
              print *,'Node centered data (no ghost zones)'
            else
              print *,'WARNING: Array is the wrong size'
!              stop
            endif

            print *,'Dimensions are ',sizes(1),sizes(2),sizes(3)

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get attribute string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'One of the attribute string labels is > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(string,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read field variable attributes'
              stop
            endif

            print *,'label is ',string(1:l_label)
            print *,'units are ',units(1:l_unit)
            print *,'short name is ',format(1:l_format)
            print *,'component is ',coordsys(1:l_coord)

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range'
              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is from ',fmin,' to ',fmax
              allocate(data3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,data3d)
              deallocate(data3d)
            else
              print *,'Range is from ',dfmin,' to ',dfmax
              allocate(ddata3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,ddata3d)
              deallocate(ddata3d)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read field variable'
              stop
            endif

            print *,'Field variable read in'
            print *,bk

          end do  ! loop on field variables

        else if (idata(1).eq.2) then
!
!---------------------------------+
!                                 |
!         Structured mesh         |
!                                 |
!---------------------------------+
!
!         Read coordinates as 3 separate 3-D arrays.

          print *,'Structured mesh'
          print *,'Reading coordinate arrays'
          print *,bk
!
! X
!
          print *,'Reading X coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 1st coord array'
            stop
          endif

          if (rank .ne. 3) then
            print *,'Nodal coords must be stored as 3 separate ',
     1              '3-D arrays'
            stop
          endif

          print *,'Coordinate dimensions are ',sizes(1),sizes(2),
     1            sizes(3)

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     1              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 1st coord array attributes'
            stop
          endif

          print *,'Label of 1st coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)

          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 1st coord array'
c            stop
          endif

          n1 = sizes(1)
          n2 = sizes(2)
          n3 = sizes(3)
          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(c1(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,c1)
            deallocate(c1)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dc1(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,dc1)
            deallocate(dc1)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 1st coord array'
            stop
          endif

          print *,'X coordinates read in'
          print *,bk 
!
! Y
!
          print *,'Reading Y coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 2nd coord array'
            stop
          endif

          if (rank .ne. 3) then
            print *,'Nodal coords must be stored as 3 separate ',
     1              '3-D arrays'
            stop
          endif

          if (sizes(1).ne.n1 .or. sizes(2).ne.n2 .or. sizes(3).ne.n3) 
     1      then
            print *,'Dimensions ',sizes(1),sizes(2),sizes(3),
     1              ' should be ',n1,n2,n3
            stop
          endif

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     1              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 2nd coord array attributes'
            stop
          endif

          print *,'Label of 2nd coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)

          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 2nd coord array'
c            stop
          endif

          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(c2(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,c2)
            deallocate(c2)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dc2(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,dc2)
            deallocate(dc2)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 2nd coord array'
            stop
          endif

          print *,'Y coordinates read in'
          print *,bk 
!
! Z
!
          print *,'Reading Z coordinate values'

          ret = dsgdims(filename,rank,sizes,maxrank) 
          if (ret .ne. 0) then
            print *,'Failed to read dims of 3rd coord array'
            stop
          endif

          if (rank .ne. 3) then
            print *,'Nodal coords must be stored as 3 separate ',
     1              '3-D arrays'
            stop
          endif

          if (sizes(1).ne.n1 .or. sizes(2).ne.n2 .or. sizes(3).ne.n3) 
     1      then
            print *,'Dimensions ',sizes(1),sizes(2),sizes(3),
     1              ' should be ',n1,n2,n3
            stop
          endif

          ret = dsgdaln(l_label,l_unit,l_format,l_coord)
          if (ret .ne. 0) then
            print *,'Could not get coordinate string lengths'
            stop
          endif

          if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
            print *,'One of the coordinate string labels is > ',nchar,
     1              ' chars'
            stop
          endif

          ret = dsgdast(label,units,format,coordsys)
          if (ret .ne. 0) then
            print *,'Failed to read 3rd coord array attributes'
            stop
          endif

          print *,'Label of 3rd coord array is ',label(1:l_label)
          print *,'Units are ',units(1:l_unit)

          if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
            print *,'Format should be a null string'
            stop
          endif

          if (version .eq. 1.00) then
            print *,'Coord sys is ',coordsys(1:l_coord)
          else

            if (coordsys(1:l_coord) .eq. 
     1          material(1:l_mat) // '|' // blockname(1:l_block)) then
              print *,'material is ',coordsys(1:l_coord)
            else
              print *,'material should be of the form gas|block'
              print *,'material is ',coordsys(1:l_coord)
            endif

          endif

          ret = dsgnt(num_type)
          if (ret .ne. 0) then
            print *,'Failed to get data type'
            stop
          endif

          if (num_type .lt. 5 .or. num_type .gt. 6) then
            print *,'Wrong data type: ',num_type
            stop
          endif

          if (num_type .eq. 5) then
            ret = dsgrang(fmax,fmin)
          else
            ret = dsgrang(dfmax,dfmin)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read range of 3rd coord array'
c            stop
          endif

          if (num_type .eq. 5) then
            print *,'Range is ',fmin,' to ',fmax
            allocate(c3(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,c3)
            deallocate(c3)
          else
            print *,'Range is ',dfmin,' to ',dfmax
            allocate(dc3(n1,n2,n3))
            ret = dsgdata(filename,rank,sizes,dc3)
            deallocate(dc3)
          endif

          if (ret .ne. 0) then
            print *,'Failed to read 3rd coord array'
            stop
          endif

          print *,'Z coordinates read in'
          print *,bk 

! Prepare to read in field variables

          print *,'Reading field variables'
          print *,bk

          do l=1,100

            ret = dsgdims(filename,rank,sizes,maxrank)
            if (ret .ne. 0) then
!              print *,'Failed to read field variable dimensions'
              stop
            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (rank.eq.1 .and. (num_type.eq.24 .or. num_type.eq.4))   
     1        exit  ! This could be a block header

            if (rank .ne. 3) then
              print *,'Field variable arrays must be of rank 3, read ',
     1                rank
              stop
            endif

            if (sizes(1)-2*idata(2).eq.n1-1) then
              print *,'Element centered data'
            else if (sizes(1)-2*idata(2).eq.n1) then
              print *,'Node centered data'
            else
              print *,'WARNING: Array is the wrong size'
!              stop
            endif

            print *,'Dimensions are ',sizes(1),sizes(2),sizes(3)

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get attribute string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'One of the attribute string labels is > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(string,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read field variable attributes'
              stop
            endif

            print *,'label is ',string(1:l_label)
            print *,'units are ',units(1:l_unit)
            print *,'short name is ',format(1:l_format)
            print *,'component is ',coordsys(1:l_coord)

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range'
              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is from ',fmin,' to ',fmax
              allocate(data3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,data3d)
              deallocate(data3d)
            else
              print *,'Range is from ',dfmin,' to ',dfmax
              allocate(ddata3d(sizes(1),sizes(2),sizes(3)))
              ret = dsgdata(filename,rank,sizes,ddata3d)
              deallocate(ddata3d)
!
! Extract head end pressure for block_0001 of lab scale rocket
!
              if ((TRIM(blockname) == 'block_0001' .or. 
     1             TRIM(blockname) == '0001') .and.
     2            (TRIM(format) == 'pressure' .or.
     3             TRIM(format) == 'press')) then
!                it = INDEX(string,'=')
!                write(8,*) TRIM(string(it+1:)),dfmax
!                write(0,*) TRIM(blockname),TRIM(string(it+1:)),dfmax
!                it = INDEX(filename,'_')
!                write(8,*) TRIM(filename(it+1:it+9)),'. ',dfmax
!                write(0,*) TRIM(blockname),' ',TRIM(filename(it+1:it+9)),
!     1                     ,'. ',dfmax
                it = INDEX(filename,'_')
! Convert text format time to integer microseconds
                select case (citime(1:2))
                case ('01')
                  expi = 1
                case ('02')
                  expi = 2
                case ('03')
                  expi = 3
                case ('04')
                  expi = 4
                case ('05')
                  expi = 5
                case ('06')
                  expi = 6
                case ('07')
                  expi = 7
                case ('08')
                  expi = 8
                case ('09')
                  expi = 9
                case ('10')
                  expi = 10
                case ('11')
                  expi = 11
                case ('12')
                  expi = 12
                end select
                write(8,181) '0',citime(3:9),'E',expi-9,' ',dfmax
181             Format (a1,a7,a1,SP,i3.2,a,SS,e15.8)
                write(0,*) TRIM(blockname),' ',TRIM(filename(it+1:it+9)),
     1                     ,'. ',dfmax
              endif
            endif

            if (ret .ne. 0) then
              print *,'Failed to read field variable'
              stop
            endif

            print *,'Field variable read in'
            print *,bk

          end do  ! loop on field variables

        else if (idata(1).ge.3 .and. idata(1).le.5) then
!
!----------------------------------+
!                                  |
!         Unstructured mesh        |
!                                  |
!----------------------------------+
!
! Read nodal coordinates as 3 separate 1-D arrays.

          print *,'Unstructured mesh'

          if (filename .eq. filegeom) then
            print *,'Reading nodal coordinate arrays'
            print *,bk
!
! X
!
            print *,'Reading X coordinate values'

            ret = dsgdims(filename,rank,sizes,maxrank) 
            if (ret .ne. 0) then
              print *,'Failed to read dims of 1st coord array'
              stop
            endif

            if (rank .ne. 1) then
            print *,'Nodal coords must be stored as 3 separate ',
     1                '1-D arrays'
              stop
            endif

            print *,'Number of nodes is ',sizes(1)

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get coordinate string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'A the coordinate string label has > ',nchar,
     1                'chars'
              stop
            endif

            ret = dsgdast(label,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read 1st coord array attributes'
              stop
            endif

            print *,'Label of 1st coord array is ',label(1:l_label)
            print *,'Units are ',units(1:l_unit)

            if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
              print *,'Format should be a null string'
              stop
            endif

            if (version .eq. 1.00) then
              print *,'Coord sys is ',coordsys(1:l_coord)
            else

              if (coordsys(1:l_coord) .eq. 
     1            material(1:l_mat) // '|' // blockname(1:l_block)) then
                print *,'material is ',coordsys(1:l_coord)
              else
                print *,'material should be of the form gas|block'
                print *,'material is ',coordsys(1:l_coord)
              endif

            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range of 1st coord array'
c              stop
            endif

            nn = sizes(1)
            if (num_type .eq. 5) then
              print *,'Range is ',fmin,' to ',fmax
              allocate(x1(nn))
              ret = dsgdata(filename,rank,sizes,x1)
              deallocate(x1)
            else
              print *,'Range is ',dfmin,' to ',dfmax
              allocate(dx1(nn))
              ret = dsgdata(filename,rank,sizes,dx1)
              deallocate(dx1)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read 1st coord array'
              stop
            endif

            print *,'X coordinates read in'
            print *,bk 
!
! Y
!
            print *,'Reading Y coordinate values'

            ret = dsgdims(filename,rank,sizes,maxrank) 
            if (ret .ne. 0) then
              print *,'Failed to read dims of 2nd coord array'
              stop
            endif

            if (rank .ne. 1) then
              print *,'Nodal coords must be stored as 3 separate ',
     1                '1-D arrays'
              stop
            endif

            if (sizes(1) .ne. nn) then
              print *,'Number of nodes is ',sizes(1), ' should be ',nn
              stop
            endif

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get coordinate string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'A coordinate string label has > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(label,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read 2nd coord array attributes'
              stop
            endif

            print *,'Label of 2nd coord array is ',label(1:l_label)
            print *,'Units are ',units(1:l_unit)

            if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
              print *,'Format should be a null string'
              stop
            endif

            if (version .eq. 1.00) then
              print *,'Coord sys is ',coordsys(1:l_coord)
            else

              if (coordsys(1:l_coord) .eq. 
     1            material(1:l_mat) // '|' // blockname(1:l_block)) then
                print *,'material is ',coordsys(1:l_coord)
              else
                print *,'material should be of the form gas|block'
                print *,'material is ',coordsys(1:l_coord)
              endif

            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range of 2nd coord array'
c              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is ',fmin,' to ',fmax
              allocate(x2(nn))
              ret = dsgdata(filename,rank,sizes,x2)
              deallocate(x2)
            else
              print *,'Range is ',dfmin,' to ',dfmax
              allocate(dx2(nn))
              ret = dsgdata(filename,rank,sizes,dx2)
              deallocate(dx2)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read 2nd coord array'
              stop
            endif

            print *,'Y coordinates read in'
            print *,bk 
!
! Z
!
            print *,'Reading Z coordinate values'

            ret = dsgdims(filename,rank,sizes,maxrank) 
            if (ret .ne. 0) then
              print *,'Failed to read dims of 3rd coord array'
              stop
            endif

            if (rank .ne. 1) then
              print *,'Nodal coords must be stored as 3 separate ',
     1                '1-D arrays'
              stop
            endif

            if (sizes(1) .ne. nn) then
              print *,'Number of nodes is ',sizes(1), ' should be ',nn
              stop
            endif

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get coordinate string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'A coordinate string label has > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(label,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read 3rd coord array attributes'
              stop
            endif

            print *,'Label of 3rd coord array is ',label(1:l_label)
            print *,'Units are ',units(1:l_unit)

            if (format(1:1) .ne. bk .and. format(1:1) .ne. ' ') then
              print *,'Format should be a null string'
              stop
            endif

            if (version .eq. 1.00) then
              print *,'Coord sys is ',coordsys(1:l_coord)
            else

              if (coordsys(1:l_coord) .eq. 
     1            material(1:l_mat) // '|' // blockname(1:l_block)) then
                print *,'material is ',coordsys(1:l_coord)
              else
                print *,'material should be of the form gas|block'
                print *,'material is ',coordsys(1:l_coord)
              endif

            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range of 3rd coord array'
c              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is ',fmin,' to ',fmax
              allocate(x3(nn))
              ret = dsgdata(filename,rank,sizes,x3)
              deallocate(x3)
            else
              print *,'Range is ',dfmin,' to ',dfmax
              allocate(dx3(nn))
              ret = dsgdata(filename,rank,sizes,dx3)
              deallocate(dx3)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read 3rd coord array'
              stop
            endif

            print *,'Z coordinates read in'
            print *,bk 

! Read connectivity data

            print *,'Looking for connectivity data'

            ret = dsgdims(filename,rank,sizes,maxrank)
            if (ret .ne. 0) then
              print *,'Failed to read dimensions for connectivity data'
              stop
            endif

            if (rank .eq. 2) then

              if (sizes(2) .gt. 8 .or. sizes(2) .lt. 3) then
                print *,'Connectivity should be (# elements,# vertices)'
                stop
              endif

              print *,'number of elements is ',sizes(1)

              ne = sizes(1)
              print *,'max number of vertices is ',sizes(2)

              ret = dsgdaln(l_label,l_unit,l_format,l_coord)
              if (ret .ne. 0) then
                print *,'Could not get attribute string lengths'
                stop
              endif

              if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
                print *,'One of the attribute string labels is > ',
     1                   nchar,' chars'
                stop
              endif

              ret = dsgnt(num_type)
              if (ret .ne. 0) then
                print *,'Failed to get data type'
                stop
              endif

              if (num_type .ne. 24) then
                print *,'Connectivity array must be 32 bit integers'
                stop
              endif

              ret = dsgdast(label,units,format,coordsys)
              if (ret .ne. 0) then
                print *,'Failed to read connectivity data attributes'
                stop
              endif

              print *,'label is ',label(1:l_label)
              print *,'units are ',units(1:l_unit)
              print *,'format is ',format(1:l_format)
              if (version .eq. 1.00) then
                print *,'coordsys is ',coordsys(1:l_coord)
              else

                if (coordsys(1:l_coord) .eq. 
     1              material(1:l_mat) // '|' // blockname(1:l_block)) 
     2              then
                  print *,'material is ',coordsys(1:l_coord)
                else
                  print *,'material should be of the form gas|block'
                  print *,'material is ',coordsys(1:l_coord)
                endif

              endif

              allocate(conn(ne,sizes(2)))

              ret = dsgdata(filename,rank,sizes,conn)
              if (ret .ne. 0) then
                print *,'Failed to read connectivity array'
                stop
              endif

              print *,'Connectivity data read in'
              print *,bk

              deallocate(conn)

            else if (rank.eq.1) then
              print *,'Found no connectivity data; assuming point data'
              print *,bk
              skip_rank = .true.
            else
              print *,'Connectivity data must be rank 2, but I read ',
     1                rank
              stop
            endif

          else
            print *,'Skip reading nodal coordinates and connectivity'
            print *,bk
          endif

! Prepare to read in field variables

          print *,'Reading field variables'
          print *,bk

          do l=1,100

            if (skip_rank) then
              skip_rank = .false.
            else
              ret = dsgdims(filename,rank,sizes,maxrank)
              if (ret .ne. 0) then
!                print *,'Failed to read field variable dimensions'
                stop
              endif
            endif

            ret = dsgnt(num_type)
            if (ret .ne. 0) then
              print *,'Failed to get data type'
              stop
            endif

            if (rank.eq.1 .and. (num_type.eq.24 .or. num_type.eq.4))   
     1        exit  ! This could be a block header

            if (rank .ne. 1) then
              print *,'Field variable arrays must be of rank 1, read ',
     1                rank
              stop
            endif

            if (filename .eq. filegeom) then

              if (sizes(1) .eq. ne) then
                print *,'Element centered data'
              else if (sizes(1) .eq. nn) then
                print *,'Node centered data'
              else
                print *,'WARNING: Array is the wrong size'
!                stop
              endif

            else
              print *,'Field variable array length is ',sizes(1)
            endif

            ret = dsgdaln(l_label,l_unit,l_format,l_coord)
            if (ret .ne. 0) then
              print *,'Could not get attribute string lengths'
              stop
            endif

            if (max(l_label,l_unit,l_format,l_coord) >= nchar) then
              print *,'One of the attribute string labels is > ',nchar,
     1                ' chars'
              stop
            endif

            ret = dsgdast(string,units,format,coordsys)
            if (ret .ne. 0) then
              print *,'Failed to read field variable attributes'
              stop
            endif

            print *,'label is ',string(1:l_label)
            print *,'units are ',units(1:l_unit)
            print *,'short name is ',format(1:l_format)
            print *,'component is ',coordsys(1:l_coord)

            if (num_type .lt. 5 .or. num_type .gt. 6) then
              print *,'Wrong data type: ',num_type
              stop
            endif

            if (num_type .eq. 5) then
              ret = dsgrang(fmax,fmin)
            else
              ret = dsgrang(dfmax,dfmin)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read range'
              stop
            endif

            if (num_type .eq. 5) then
              print *,'Range is from ',fmin,' to ',fmax
              allocate(data(sizes(1)))
              ret = dsgdata(filename,rank,sizes,data)
              deallocate(data)
            else
              print *,'Range is from ',dfmin,' to ',dfmax
              allocate(ddata(sizes(1)))
              ret = dsgdata(filename,rank,sizes,ddata)
              deallocate(ddata)
            endif

            if (ret .ne. 0) then
              print *,'Failed to read field variable'
              stop
            endif

            print *,'Field variable read in'
            print *,bk

          end do  ! loop on field variables

        else
          print *,'Invalid mesh type: ',idata(1)
          stop
        endif  ! mesh type

      end do  ! loop over mesh blocks

      end

