#!/bin/sh
# Bourne Shell script to replace %...%-style boundary conditions
# in a .flo file with actual hard numbers.

echo "Fix_bc.sh: postprocesses .flo file to insert actual boundary"
echo "  conditions instead of %...% placeholders."

if [ $# -ne 2 ]
then
	echo "Usage: fix_bc.sh <in .flo file> <out .flo file>"
	exit 1
fi

in="$1"
out="$2"
tmp="tmp1.flo"
tmp2="tmp2.flo"

rm -fr $tmp $tmp2 2> /dev/null > /dev/null
cp $in $tmp

while [ true ]
do
	echo
	echo "Please enter the boundary condition type to replace,"
	echo " or 'q' to exit.  Ex: unspecified, generic1, custom"
	bc="q"
	read bc
	[ "$bc" = "q" ] && break
	echo "Enter the rocflo boundary condition number for $bc:"
	num="0"
	read num
	echo "Enter a human-readable description for $bc:"
	desc=""
	read desc
	echo "Enter any parameters needed by $bc, finishing with Control-D:"
	cat > tmp_params
	params=`awk -F% '{print $1"\\\\"}' tmp_params`
	rm tmp_params
	
# Use sed to splice the input data into the .flo file
	cat > sedScript <<END
s/%bc_$bc%/$num/g
s/%bc_${bc}_desc%/$desc/g
s/%bc_${bc}_params%/$params
/g
END
	echo "Replacing $bc by ($num/$desc/$params) in $in file..."
	sed -f sedScript $tmp > $tmp2
	mv $tmp2 $tmp
done

cp $tmp $out
