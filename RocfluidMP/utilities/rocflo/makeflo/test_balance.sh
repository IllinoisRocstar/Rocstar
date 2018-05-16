#!/bin/sh

prog=$1

for blocks in 142 145 150 155 160 165 170 180 190 200 225 250 275 300 325 350 400 450 500 550 600 650 700 750 800 900 1000 1500 2000 2500 3000
do
	$prog  test/test.grd $blocks /dev/null /dev/null > tmp.$$
	echo $blocks '->' `grep "Heaviest" tmp.$$` 
	rm tmp.$$
done
