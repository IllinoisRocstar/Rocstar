#!/bin/sh

################################################################################
#
# Purpose: Remove carriage return (ctrl M) from files. 
#
# Description: None.
#
# Notes: 
#   1. IMPORTANT: This file may not be displayed properly in CVSWeb because the
#      ctrl-M character is interpreted. Nevertheless, the file is stored 
#      correctly in CVS, so DO NOT CHANGE IT.
# 
# Author: Andreas Haselbacher
#
################################################################################

TEMPFILE=rmcm.tmp

sed 's///g' $1 > $TEMPFILE
mv $TEMPFILE $1

################################################################################
#
# RCS Revision history:
#
# $Log: rmcm.sh,v $
# Revision 1.2  2005/09/19 18:27:38  haselbac
# Added note
#
# Revision 1.1  2005/09/19 16:30:48  haselbac
# Initial revision
#
################################################################################
