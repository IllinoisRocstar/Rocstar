#!/bin/tcsh

set OutFile=$1
set TmpOut=${OutFile}_tmp.txt
set BINDIR=${3}

# Test the serial example program output
echo "In this random serial_tests.csh script"

${BINDIR}/sep -o test_sep_out.txt Makefile 

if (! -e "test_sep_out.txt") then
   echo "Couldn't find test_sep_out.txt to diff"
   exit 1
endif

if (! -e "Makefile") then
   echo "Couldn't find Makefile to diff"
   exit 1
endif

set STEST=`diff Makefile test_sep_out.txt`
rm -f test_sep_out.txt
printf "ExampleProgram:Works=" >> ${TmpOut}
if( "$STEST" == "") then
  printf "1\n" >> ${TmpOut}
else
  printf "0\n" >> ${TmpOut}
endif
cat ${TmpOut} >> ${OutFile}
rm -f ${TmpOut}
exit 0
