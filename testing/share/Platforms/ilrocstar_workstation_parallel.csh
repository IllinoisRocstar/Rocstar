#!/bin/tcsh

set RESULTSFILE = ${1}
set SRCDIR = ${2}
set BINDIR = ${3}

rm -f tmpresults_1.txt
mpirun -np 4 ${BINDIR}/rocstar_parallel_test -o tmpresults_1.txt
@ i = 1
while($i <= 240)
    @ i += 1
    if( -e tmpresults_1.txt ) then
        @ i += 241;
    else
        sleep 30;
    endif
end
mv tmpresults_1.txt ${RESULTSFILE}
exit 0
