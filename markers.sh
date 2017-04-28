cut -d ' ' --complement -f -1,-2,-3,-4,-5 $1.picked > /tmp/$1
cp /tmp/$1 $1.picked
rm /tmp/$1	