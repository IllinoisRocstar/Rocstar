!!! $Id: rocmanf90.h,v 1.3 2004/05/05 14:15:29 brtnfld Exp $

!!! **********************************************************************
!!! *                                                                    *
!!! * Author:      X. Jiao, A. Haselbacher, and M. Breitenfeld           *
!!! * Created on:  Nov. 14, 2002                                         *
!!! * Purpose:     Defines parameters in GENx that are used in more than *
!!! *         one physics component. Added initially for computing       *
!!! *         integrals in sanity checking (such as mass conservation).  *
!!! **********************************************************************


!!! Name of the variable containing the integrals to be registered with Roccom.
CHARACTER(*), PARAMETER :: MAN_INTEG_NAME      = "integrals"

!!! Name of the function that computes the integrals to be registered with Roccom.
CHARACTER(*), PARAMETER :: MAN_COMP_INTEG_NAME = "compute_integrals"

!!! Length of the array for storing integrals
INTEGER, PARAMETER :: MAN_INTEG_SIZE=9

!!! Indices of the entries in the array
INTEGER, PARAMETER :: MAN_INTEG_VOL=1, MAN_INTEG_MASS=2    ! Volume and mass
INTEGER, PARAMETER :: MAN_INTEG_XMOM=3, MAN_INTEG_YMOM=4, MAN_INTEG_ZMOM=5 ! Momemtum
INTEGER, PARAMETER :: MAN_INTEG_ENER=6    ! Energy
INTEGER, PARAMETER :: MAN_INTEG_IBAREA=7  ! Area of burning f-s interface
INTEGER, PARAMETER :: MAN_INTEG_INBAREA=8 ! Area of non-burning f-s interface
INTEGER, PARAMETER :: MAN_INTEG_VOL_UND=9 ! Undeformed volume

!!! Double precision data type
INTEGER, PARAMETER :: MAN_DBL = SELECTED_REAL_KIND(P=14,R=30)
