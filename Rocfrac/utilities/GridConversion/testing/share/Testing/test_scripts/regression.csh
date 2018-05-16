#!/bin/tcsh

#Enter necessary filename variables here
set OutFile=$1
set TmpOut=${OutFile}_tmp.txt 
set InputDir=#ENTER THE NAME OF THE INPUT DATA DIRECTORY HERE
set Outputs=#ENTER EXEC OUTPUT FILE NAME(S) HERE
            #IF ONE USE FORMAT Outputs=Filename
            #IF MULTIPLE USE FORMAT Outputs=(Filename1 Filename2 Filename3 ...)
set OutputsCheck=#ENTER SAVED OUTPUT FILE TO CHECK AGAINST HERE
            #IF ONE USE FORMAT OutputsCheck=CheckFilename
            #IF MULTIPLE USE FORMAT OutputsCheck=(CheckFilename1 CheckFilename2 
            # CheckFilename3 ...)
            #THERE MUST BE ONE CHECK FILE FOR EACH OUTPUT FILE in Outputs AND
            #THEY MUST CORRESPOND! (i.e., Filname1 is diffed with CheckFilename1).
set TestName=#ENTER TEST NAME FOR RUNTEST TO PARSE HERE

#Remove old test InputDir if present
if( -d  ${InputDir}) then
  echo "removing ${InputDir} directory"
  rm -r ${InputDir}
endif

#Make InputDir directory to run test in
mkdir ${InputDir}
cd ${InputDir}

#Copy input data into InputDir
cp -r $2/share/Testing/test_data/${InputDir}/* .

#Run executable to generate output data
#ENTER YOUR COMMAND FOR RUNNING YOUR EXECUTABLE HERE
#(HINT: $3 IS THE PATH TO THE BIN DIRECTORY GIVEN WHEN CALLING RUNTEST)

#Make sure the necesary output was generated
foreach file (${Outputs})
  if( ! -e ${file} ) then
    echo "No ${file} results file from run!"
    exit 1
  endif
end

#variable for test passing
@ result = 1

#diff the new output file with the saved one to check
#(This uses our own speical diff (diffdatafiles) that
#can compare numbers with a tolerance. See documentation
#for more information on how to use it.)
@ i = 1
foreach file (${Outputs})
  $4/diffdatafiles ${file} $OutputsCheck[$i] -t 1.0e-10 -n
  if($? != 0) then
    echo "${file} differs from $OutputsCheck[$i]"
    @ result = 0
  endif
  @ i += 1
end

#print test results to OutFile
printf "${TestName}=" >> ${TmpOut}
printf "$result\n" >> ${TmpOut}
cat ${TmpOut} >> ../${OutFile}
cd ..
rm -r ${InputDir}

exit 0
