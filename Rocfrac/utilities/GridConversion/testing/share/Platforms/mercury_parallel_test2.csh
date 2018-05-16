#!/bin/tcsh

rm -f tmpresults_2.txt
cat <<EOF > ./gridconversion_parallel_test_batch.csh
#!/bin/tcsh
#
#PBS -l nodes=2:ppn=8
#PBS -l walltime=00:30:00 
#PBS -j oe
#PBS -o gridconversion_parallel_test_batch_output
#PBS -A IRRD

cd \${PBS_O_WORKDIR}
mpiexec -np 16 ./bin/pepi -o tmpresults_2.txt 1000000
EOF
qsub gridconversion_parallel_test_batch.csh
@ i = 1
while($i <= 720)
    @ i += 1
    if( -e tmpresults_2.txt ) then
        @ i += 721;
    else
        sleep 10;
    endif
end
printf "PEPI:Runs=" >> $1
@ err = 0
if ( -e tmpresults_2.txt ) then
   printf "1\n" >> $1
else
   printf "0\n" >> $1
   @ err += 1
endif
set RESULTS=`cat tmpresults_2.txt | grep 3.141592653589`
printf "PEPI:Works=" >> $1
if ( "$RESULTS" == "") then
   printf "0\n" >> $1
   @ err += 1
else
   printf "1\n" >> $1
endif  
rm -f tmpresults_2.txt
rm -f ./gridconversion_parallel_test_batch.csh
rm -f gridconversion_parallel_test_batch_output
exit ${err}
