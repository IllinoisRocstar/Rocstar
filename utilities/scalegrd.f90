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
      IMPLICIT NONE
      INTEGER, PARAMETER ::  double = SELECTED_REAL_KIND(P=14,R=30)
      real(kind=double), allocatable :: x(:),y(:),z(:)
      INTEGER, allocatable :: blockDims(:,:)
      character*50 infile,outfile
      real(kind=double) scaleFactor
      INTEGER numBlocks, numNodes, blocks, i, j, k
      
      Write(*,*) "**********************************************************"
      Write(*,*) "This program is designed to take a GRIDGEN-generated text"
      Write(*,*) "grid file (.grd), and multiply every x, y, and z in the "
      Write(*,*) "file by a uniform factor. It was written to translate a "
      Write(*,*) "grid in millimeters to one in meters. "      
      Write(*,*) "**********************************************************"

      write(*,*)'Enter Input File Name:'
      read(*,*)infile
      write(*,*)'Enter Output File Name:'
      read(*,*)outfile
      write(*,*)'Enter Multiplicative Scale Factor for all x,y,z:'
      read(*,*)scaleFactor

      open(unit=10,file=infile,form='formatted',action="read")
      open(unit=11,file=outfile,form='formatted',action="write")
	  
	  WRITE(*,*) "Ready to read numBlocks"
      READ(10,*) numBlocks
	  WRITE(*,*) "numblocks = ", numBlocks
      allocate (blockDims(numBlocks, 3))
	  
	  WRITE(*,*) "Read block dimensions."
      do blocks = 1, numBlocks
         READ(10,*) blockDims(blocks, 1), blockDims(blocks, 2), blockDims(blocks, 3)
      enddo
	  
	  WRITE(*,*) "Write new block header in new file."
      WRITE(11,102) numBlocks
      do blocks = 1, numBlocks
         WRITE(11,102) blockDims(blocks, 1), blockDims(blocks, 2), blockDims(blocks, 3)
      enddo
	  
	  WRITE(*,*) "Read and write block data."
      do blocks = 1,numBlocks
	    
         WRITE(*,*) "Processing block:", blocks, " out of ", numBlocks, " block(s)"
	    
	 numNodes = blockDims(blocks, 1)*blockDims(blocks, 2)*blockDims(blocks, 3)
	    
	 allocate (x(numNodes), y(numNodes), z(numNodes))
	    
	 READ(10,*) (x(i), i=1, numNodes)
	 READ(10,*) (y(j), j=1, numNodes)
	 READ(10,*) (z(k), k=1, numNodes)
	   
	 WRITE(11,101) (x(i)*scaleFactor, i=1, numNodes)
	 WRITE(11,101) (y(j)*scaleFactor, j=1, numNodes)
	 WRITE(11,101) (z(k)*scaleFactor, k=1, numNodes)
	    
	 DEALLOCATE(x, y, z)

      enddo

      close (10)
      close (11)   

101   format(4(ES19.11))
102   format(I8, 1X, I8, 1X, I8)

      stop
      end






