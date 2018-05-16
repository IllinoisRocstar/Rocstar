!******************************************************************************
!
! Purpose: calculation of pointers to grid cells/nodes.
!
! Description: - IndIJ     = pointer to grid cell/node based on i,j indices
!              - IndIJK    = pointer to grid cell/node based on i,j,k indices
!
! Notes: subroutines RFLO_GetCellOffset or RFLO_GetNodeOffset must be called
!        BEFORE the functions are used in order to set the correct offsets.
!
!******************************************************************************
!
! $Id: Indexing.h,v 1.1 2003/05/15 03:10:38 jblazek Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

#define IndIJK(i,j,k,iOff,ijOff) ((i)+((j)-1)*(iOff)+((k)-1)*(ijOff))
#define IndIJ(i,j,iOff) ((i)+((j)-1)*(iOff))
