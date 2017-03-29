#!/bin/bash

# To convert in a directory
# for file in ./data1/*.tri; do a=${file##*/}; a=`echo $a | rev | cut -c5- | rev` ; ./convert.sh data1/$a; done
input=$1

triangles=$1.tri
vertices=$1.vert

output=$1.off

echo "OFF" > $output
a=`cat $vertices | wc -l`
a=$[$a+1]
b=`cat $triangles | wc -l`
echo "$a $b 0" >> $output

echo "0.0 0.0 0.0" >> $output
echo "$(cat $vertices)" >> $output
while read -r line
do 
	name="$line"
	echo "3 $name" >> $output
done < "$triangles"
# echo "$(cat $triangles)" >> $output

