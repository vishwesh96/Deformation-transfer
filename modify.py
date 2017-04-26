import os
import sys
b = True
w = open(sys.argv[1]+"1",'a')
lineno =0
with open(sys.argv[1]) as f:
	for line in f:
		l = line.split();
		lineno +=1;
		if lineno == 2:
			w.write(str(int(l[0])-1)+" "+str(int(l[1]))+" "+str(int(l[2]))+"\n");
		elif len(l) == 4:
			w.write(str(l[0])+" "+str(int(l[1])-1)+" "+str(int(l[2])-1)+" "+str(int(l[3])-1)+"\n");
		elif l[0]=='0.0' and l[1]=='0.0' and l[2]=='0.0' and b:
			print("ignore")
			b=False
		else:
			w.write(line)
