#!/bin/tcsh

set RESULTSFILE = ${1}
set SRCDIR = ${2}
set BINDIR = ${3}

mpirun -np 4 ${BINDIR}/pepi -o tmpresults_2.txt 1000000
printf "PEPI:Runs=" >> ${RESULTSFILE}
@ err = 0
if ( -e tmpresults_2.txt ) then
   printf "1\n" >> ${RESULTSFILE}
else
   printf "0\n" >> ${RESULTSFILE}
   @ err += 1
endif
set RESULTS=`cat tmpresults_2.txt | grep 3.141592653589`
printf "PEPI:Works=" >> ${RESULTSFILE}
if ( "$RESULTS" == "") then
   printf "0\n" >> ${RESULTSFILE}
   @ err += 1
else
   printf "1\n" >> ${RESULTSFILE}
endif  
rm -f tmpresults_2.txt
exit ${err}
