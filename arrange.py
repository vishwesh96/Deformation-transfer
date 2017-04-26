import sys
import os
import random
import string

from optparse import OptionParser

if __name__ == '__main__':
	parser = OptionParser()

	parser.add_option("-i", "--ifile", dest="i_file", help="obj.off", metavar="FILE", default="a.off")
	parser.add_option("-a", "--afile", dest="a_file", help="ans.txt", metavar="FILE", default="ans.txt")
	parser.add_option("-o", "--ofile", dest="o_file", help="output", metavar="FILE", default="out.off")
	(options, args) = parser.parse_args()

	i_file = options.i_file
	a_file = options.a_file
	o_file = options.o_file
	#o_file = *.off
	#ans = *.txt
	
	fi = open(i_file, 'r')
	fa = open(a_file, 'r')
	fo = open(o_file, 'w')

	l = fi.readline();
	fo.write(l)

	l = fi.readline();
	fo.write(l)
	
	l = l.replace('\n','')
	n = int(l.split(' ')[0])

	i = 0
	while True:
		l1 = fi.readline()
		l2 = fa.readline()

		l = l1.replace('\n','')
		if l == '' or l is None:
			break

		if i<n:
			fo.write(l2)
			i = i + 1
		else:
			fo.write(l1)

	fi.close()
	fa.close()
	fo.close()