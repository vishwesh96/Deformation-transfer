import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg.dsolve import linsolve
from scipy.sparse.linalg import *
from scipy import io

import sys
import os
import random
import string

from optparse import OptionParser

if __name__ == '__main__':
	parser = OptionParser()

	parser.add_option("-a", "--afile", dest="a_file", help="a.txt", metavar="FILE", default="A.txt")
	parser.add_option("-c", "--cfile", dest="c_file", help="c.txt", metavar="FILE", default="c.txt")

	(options, args) = parser.parse_args()

	a_file = options.a_file
	c_file = options.c_file

	fa = open(a_file, 'r')
	fc = open(c_file, 'r')

	dim_a = fa.readline().replace('\n','')
	dim_a = dim_a.split(' ')

	row_a = int(dim_a[0])
	col_a = int(dim_a[1])

	#initialize sparse matrix A
	mat_a = sps.lil_matrix((row_a, col_a), dtype=np.float64)

	while True:
		l = fa.readline().replace('\n','')
		if l is None or l == '':
			break;

		l = l.split(' ')
		mat_a[int(l[0]), int(l[1])] = np.float64(np.float64(l[2]) * 1.0)


	fa.close()

	#initilaize matrix c
	c_size = fc.readline().replace('\n','')
	c_size = c_size.split(' ')
	c_size = int(c_size[0])

	mat_c = np.zeros(c_size, dtype=np.float64)

	while True:
		l = fc.readline().replace('\n','')
		if l is None or l == '':
			break

		l = l.split(' ')
		if int(l[0])%2 == 0:
			mat_c[int(l[0])] = np.float64(np.float64(l[1]) * 1.0)

	fc.close()

	# print(mat_c)
	# print(mat_a)

	at = mat_a.transpose()

	# print(mat_a)
	# print(mat_a.transpose())

	at  = at.tocsr()
	mat_a = mat_a.tocsr()

	ata = at.dot(mat_a)
	atc = at.dot(mat_c)
	# print(mat_c)
	# print(ata)

	# print(mat_c)
	# print(at)
	# print(atc)

	# print(ata)
	io.mmwrite("AtA.mtx", ata)
	
	f = open("Atc.mtx", 'w')
	f.write(str(atc.size))
	f.write('\n')
	for a in atc:
		f.write(str(a))
		f.write('\n')

	f.close()

	ata = ata.tocsr()
	x = spsolve(ata, atc, use_umfpack = True)

# rand = np.random.rand
# mtx = sps.lil_matrix((10000000, 10000000), dtype=np.float64)
# mtx[0, :1000000] = rand(1000000)
# mtx[1, 1000000:2000000] = mtx[0, 1000000:2000000]
# mtx.setdiag(rand(10000000))
# mtx=mtx.tocsr()
# rhs=rand(100000)
# x = linsolve.spsolve(mtx, rhs)
