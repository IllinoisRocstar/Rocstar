################################################################################
#
# $Id: Makefile,v 1.14 2006/01/10 21:53:07 wasistho Exp $
#
# Purpose: Makefile for Rocflo/Rocflu-MP pre- and postprocessing utilities.
#
# Description: none.
#
# Notes: none.
#
# Copyright: (c) 2001-2005 by the University of Illinois
#
################################################################################

all:
      ifdef RFLO
	$(MAKE) -C rocflo/prep UTIL=yes
	$(MAKE) -C rocflo/to3d UTIL=yes
	$(MAKE) -C rocflo/split UTIL=yes
	$(MAKE) -C rocflo/surfgrid UTIL=yes
	$(MAKE) -C rocflo/blocks UTIL=yes
	$(MAKE) -C rocflo/post UTIL=yes
        ifdef PLAG
	  $(MAKE) -C rocpart/post UTIL=yes
        endif
        ifdef PEUL
	  $(MAKE) -C rocsmoke/post UTIL=yes
        endif
        ifdef TOFLU
	  $(MAKE) -C rocflo/toflu UTIL=yes
        endif
	$(MAKE) -C vav
      endif
      ifdef RFLU
	$(MAKE) -C rocflu/conv UTIL=yes
	$(MAKE) -C rocflu/init UTIL=yes
	$(MAKE) -C rocflu/map UTIL=yes
	$(MAKE) -C rocflu/part UTIL=yes
	$(MAKE) -C rocflu/pick UTIL=yes
	$(MAKE) -C rocflu/post UTIL=yes
      endif

clean install:
      ifdef RFLO
	$(MAKE) -C rocflo/prep $@
	$(MAKE) -C rocflo/post $@
        ifdef PLAG
	  $(MAKE) -C rocpart/post $@
        endif
        ifdef PEUL
	  $(MAKE) -C rocsmoke/post $@
        endif
        ifdef TOFLU
	  $(MAKE) -C rocflo/toflu $@
        endif 
	$(MAKE) -C rocflo/to3d $@
	$(MAKE) -C rocflo/split $@
	$(MAKE) -C rocflo/surfgrid $@
	$(MAKE) -C rocflo/blocks $@
	$(MAKE) -C vav $@
      endif 
      ifdef RFLU
	$(MAKE) -C rocflu/conv $@
	$(MAKE) -C rocflu/init $@
	$(MAKE) -C rocflu/map $@
	$(MAKE) -C rocflu/part $@
	$(MAKE) -C rocflu/pick $@
	$(MAKE) -C rocflu/post $@
      endif 

################################################################################
#
# RCS Revision history:
#
# $Log: Makefile,v $
# Revision 1.14  2006/01/10 21:53:07  wasistho
# make rflopost as default
#
# Revision 1.13  2005/04/15 15:05:48  haselbac
# Removed Charm/FEM stuff, cosmetics
#                                                             
# Revision 1.12  2004/10/19 19:22:24  haselbac                                   
# Changed so rflumap gets built and cleaned out if CHARM == 1                    
#
# Revision 1.11  2004/08/17 00:59:05  wasistho                                   
# prepared for utilities/rocflo/toflu                                            
#
# Revision 1.10  2004/03/19 21:16:56  haselbac                                   
# Added rinit for RFLU                                                           
#
# Revision 1.9  2003/09/25 15:40:22  jferry                                      
# Implented Rocsmoke post-processing                                             
#
# Revision 1.8  2003/06/04 21:53:15  haselbac                                    
# Removed cells and added pick                                                   
#
# Revision 1.7  2003/05/06 16:20:36  fnajjar                                     
# Included calls for Rocpart postprocessing tool                                 
#
# Revision 1.6  2003/04/11 15:15:18  jiao                                        
# Changed RFLU_CHARM back to CHARM.                                              
#
# Revision 1.5  2003/04/01 16:35:08  haselbac                                    
# Added rflucells, restriction for serial utils                                  
#
# Revision 1.4  2003/03/25 19:09:13  haselbac                                    
# Removed POST for Rocflu                                                        
#
# Revision 1.3  2003/03/21 15:36:07  haselbac                                    
# Added install target                                                           
#
# Revision 1.2  2003/03/20 23:45:48  jiao                                        
# ACH: Added rocflu/conv.                                                        
#
# Revision 1.1.1.1  2003/03/20 22:48:04  haselbac                                
# Initial revision                                                               
#
# Revision 1.14  2002/12/20 19:38:02  jblazek                                    
# Added tool to count blocks in a box.                                           
#
# Revision 1.13  2002/10/19 00:40:31  jblazek                                    
# Added utility (rflosurf) to write out surface grids for GenX.                  
#
# Revision 1.12  2002/10/16 21:20:30  haselbac                                   
# Added rflump                                                                   
#
# Revision 1.11  2002/09/17 13:43:00  jferry                                     
# Added Time-dependent boundary conditions                                       
#
# Revision 1.10  2002/09/10 00:00:52  f-najjar                                   
# Reactivate RocVAV compliation flag                                             
#
# Revision 1.9  2002/09/05 17:40:22  jblazek                                     
# Variable global moved into regions().                                          
#
# Revision 1.8  2002/07/12 21:50:07  jblazek                                     
# Added tool to split single grid into multiple regions.                         
#
# Revision 1.7  2002/06/22 01:13:38  jblazek                                     
# Modified interfaces to BC routines.                                            
#
# Revision 1.6  2002/06/10 22:30:13  haselbac                                    
# Fixed spacing problem, moved vav inside ifdef RFLO                             
#
# Revision 1.5  2002/06/10 18:03:55  jblazek                                     
# Added V&V to utilities/Makefile. Corrected path to libraries in vav/Makefile.  
#
# Revision 1.4  2002/01/12 00:02:49  jblazek                                     
# Added postprocessor.                                                           
#
# Revision 1.3  2001/12/22 00:09:39  jblazek                                     
# Added routines to store grid and solution.                                     
#
# Revision 1.2  2001/12/21 23:56:51  jblazek                                     
# Added utility to convert 2D grids to 3D.                                       
#
# Revision 1.1.1.1  2001/12/03 21:44:05  jblazek                                 
# Import of RocfluidMP                                                           
#
################################################################################
