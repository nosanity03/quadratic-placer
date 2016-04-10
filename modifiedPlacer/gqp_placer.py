#! /usr/bin/python3

import sys

from copy import deepcopy
import numpy as np
from numpy import array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import math
from itertools import repeat

class MyError(Exception):
	def __init__(self, value):
		self.value = value

class mothercore:
	def __init__(self, side):
		self.gate = dict([])
		self.pad = dict([])
		self.numG = 0
		self.numP = 0
		self.numN = 0
		self.side = side
		self.gateX = dict([])
		self.gateY = dict([])
		self.nets = dict([])
		self.sortedgate = [] # [None]*numG

	def get_numG(self):
		return deepcopy(self.numG)

	def get_numP(self):
		return deepcopy(self.numP)

	def get_numN(self):
		return deepcopy(self.numN)

	def set_side(self, side):
		self.side = deepcopy(side)
	
	def get_gate(self):
		return deepcopy(self.gate)

	def get_gateconnections(self, gatenum):
		return deepcopy(self.gate[gatenum])

	def get_gatelocation(self, gatenum):
		return deepcopy([self.gateX[gatenum], self.gateY[gatenum]])

	def get_pad(self):
		return deepcopy(self.pad)

	def get_padlocation(self, padnum):
		return deepcopy(self.pad[padnum])

	def add_gate(self, gatenum, listofconnections):
		self.gate[gatenum] = deepcopy(listofconnections)
		self.numG += 1
		return 1

	def remove_gate(self, gatenum):
		del self.gate[gatenum]
		del self.gateX[gatenum]
		del self.gateY[gatenum]
		self.numG -= 1
		return self.numG

	def add_pad(self, padnum, listofconnections):
		self.pad[padnum] = deepcopy(listofconnections)
		self.numP += 1
		return 1

	def add_net(self, netnum, connection, gateorpad):
		# 0 for gate, 1 for pad

		## if netnum doesn't already exist in dictionary
		if netnum not in self.nets:
			self.nets[netnum] = [[],[]] # first list for gates, second for pads
			self.numN += 1

		if (gateorpad == 1):
			self.nets[netnum][1].append(deepcopy(connection))
		else:
			self.nets[netnum][0].append(deepcopy(connection))
		return 1

	def add_netconns(self, netnum, listofconnections):
		# 0 for gate, 1 for pad

		## if netnum doesn't already exist in dictionary
		if netnum not in self.nets:
			self.nets[netnum] = [[],[]] # first list for gates, second for pads
			self.numN += 1

		self.nets[netnum] = deepcopy(listofconnections)
		return 1

	def get_otherconns(self, netnum, gatenum):
		# returns the other entities connected to Net "netnum", other than "gatenum"
		# returns list of 2 lists: first list for gates, second for pads
		netconns = deepcopy(self.nets[netnum])

		# gatenum is not in netconnections
		if (gatenum not in netconns[0]): return [[],[]]

		# gatenum is in netconnections
		netconns[0].remove(gatenum)
		return netconns

	def get_nets(self):
		return deepcopy(self.nets)


	def add_location(self, x, y, data):	# data is required to know the keys of x,y values
		if (len(data) == len(x)):
			for l in range(0,len(data)):
				xloc = x[l]
				yloc = y[l]
				if x[l] > self.side[1]:
					xloc = self.side[1]
				elif x[l] < self.side[0]:
					xloc = self.side[0]
				if y[l] > self.side[3]:
					yloc = self.side[3]
				elif y[l] < self.side[2]:
					yloc = self.side[2]
				self.gateX[data[l]] = xloc
				self.gateY[data[l]] = yloc
			return 1
		else:
			return 0

	def add_sorted(self, sort):
		self.sortedgate = deepcopy(sort)
		return 1

	def get_sorted(self):
		return deepcopy(self.sortedgate)

	def get_location(self):
		l = [None]*2
		l[0] = deepcopy(self.gateX)
		l[1] = deepcopy(self.gateY)
		return l

def create(filename):
	print('Creating data structure from file ...')
	cf = open(filename,'r')

	## Find the number of gates, nets and pads for creating object
	# for i, line in enumerate(cf):
	# 	line = line.strip()
	# 	if line == '': continue
	# 	v = line.split()
	# 	if i == 0:
	# 		numG = int(v[0])
	# 		numN = int(v[1])
	# 		numP = 1
	# 	elif i == numG+1:
	# 		numP = int(v[0])
	# 		break
	# 	else:
	# 		continue
	#cf.close()

	## New object
	returncore = mothercore([0,100,0,100])

	## Add gates, nets and pads to object from file
	#cf = open(filename,'r')
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()

		if i == 0:
			numG = int(v[0])
			numN = int(v[1])
			numP = 1
		elif i == numG+1:
			numP = int(v[0])

		## Add gates and nets
		elif i in range(1,numG+1):
			nets = [int(v[1])] # number of nets
			for j in range(2,int(v[1])+2):
				nets.append(int(v[j]))
				returncore.add_net(int(v[j]), int(v[0]), 0)
			returncore.add_gate(int(v[0]), nets)

		## Add pads and nets
		elif i in range(numG+2,numG+numP+2):
			nets = []
			for j in range(1,4):
				nets.append(int(v[j]))
			returncore.add_net(int(v[1]), int(v[0]), 1)
			returncore.add_pad(int(v[0]), nets)
		else:
			continue

	cf.close()
	print('Done. Added '+ str(returncore.get_numG()) + ' gates, '+ str(returncore.get_numP()) + ' pads, '+ str(returncore.get_numN()) + ' nets.')
	return returncore

def solveforx(core):
	print('Solving for locations ...')

	###########################################################################################
	### Initializations
	gate = core.get_gate()
	G = core.get_numG()
	key = list(gate.keys())
	pad = core.get_pad()
	P = core.get_numP()
	keyP = list(pad.keys())
	if (G == 0):
		return 1

	C = [0]*G
	for i in range(0,G):
		C[i] = [0]*G

	A = [0]*G
	for i in range(0,G):
		A[i] = [0]*G

	bx = [0]*G
	for i in range(0,G):
		bx[i] = 0

	by = [0]*G
	for i in range(0,G):
		by[i] = 0

	nets = core.get_nets()
	weights = dict([])

	###########################################################################################
	### Calculating weights
	print('Calculating weights ...')
	for val in nets:
		k = len(nets[val][0])+len(nets[val][1])
		#print(val, k)
		#print(len(nets[val][0]), len(nets[val][1]))
		weights[val] = 1/(k-1)
	
	#print('nets = ', nets)
	#print('weights = ', weights)
	
	###########################################################################################
	### Gate numbers may vary, this dictionary keeps them in order
	gateorder = dict([])
	for i in range(0,G):
		gateorder[key[i]] = i

	###########################################################################################
	### C and A matrices, bx and by vectors
	print('Calculating valid contributions to the cost function ...')

	for netval in nets.keys():
		weight = weights[netval]
		gates = nets[netval][0]
		numgates = len(gates)
		pads = nets[netval][1]
		numpads = len(pads)
		#print(netval, weight, gates, pads)
		if (numgates > 1):
			i = 0
			while (i < numgates-1):
				j = i+1
				while (j < numgates):
					C[gateorder[gates[i]]][gateorder[gates[j]]] = weight
					A[gateorder[gates[i]]][gateorder[gates[j]]] = -weight
					C[gateorder[gates[j]]][gateorder[gates[i]]] = weight
					A[gateorder[gates[j]]][gateorder[gates[i]]] = -weight
					A[gateorder[gates[i]]][gateorder[gates[i]]] += weight
					A[gateorder[gates[j]]][gateorder[gates[j]]] += weight
					j += 1
				i += 1

		if (numpads > 0):
			i = 0
			while (i < numgates):
				j = 0
				while (j < numpads):
					A[gateorder[gates[i]]][gateorder[gates[i]]] += weight
					bx[gateorder[gates[i]]] += weight*pad[pads[j]][1]
					by[gateorder[gates[i]]] += weight*pad[pads[j]][2]
					j += 1
				i += 1

	#for i in range(0,G):
	#	print('C ', i+1, ': ',  C[i])

	#for i in range(0,G):
	#	print('A ', i+1, ': ',  A[i])

	#print('bx : ',  bx)
	#print('by : ',  by)

	###########################################################################################
	### R, C, V matrices for sparse matrix generation
	print('Forming R, C, V matrices ...')
	R = []
	C = []
	V = []
	for i in range(0,G):
		#if (i%10 == 0): print('4 -> Gate ',i+1,' of ',G)
		for j in range(0,G):
			if (A[i][j] != 0) & (math.isnan(A[i][j]) is False) & (math.isinf(A[i][j]) is False):
				R.append(i)
				C.append(j)
				V.append(A[i][j])
		
	###########################################################################################
	### Solve for x and y vectors
	print('Solving for x and y ...')
	R = np.asarray(R)
	C = np.asarray(C)
	V = np.asarray(V)
	A = coo_matrix((V, (R, C)), shape=(G, G))
	bx = np.asarray(bx)
	by = np.asarray(by)
	#print('A = ', A)
	x = spsolve(A.tocsr(), bx)
	y = spsolve(A.tocsr(),by)
	x = x.tolist()
	y = y.tolist()
	#print('x = ', x, 'y = ', y)

	###########################################################################################
	update_location(core, x, y,key)
	print('Done.')	
	return 1
	###########################################################################################

def writeback(core,filename):
	G = core.get_location()
	wf = open(filename,'w')
	for l in range(1,core.numG+1):
		print(l," ",G[0][l]," ", G[1][l], file=wf)
	wf.close()

def assign(core, G, hORv, lORr):	# size of square, 1 horizontal, then 2 vertical
	print('Assignment ...')

	### horizontal sort
	x = 0
	var = deepcopy(G[x])           # 0 for horizontal (x) and 1 for vertical (y)
	var2 = deepcopy(G[1-x])

	# merge 2 sort keys into 1 using (100000*x+y). This will work as long
	# as there are less than 100k gates.
	for i in var.keys():
		var2[i] = var[i]*100000+var2[i]
	length = len(var.keys())
	midgate = math.floor(length/2)
	balance = length - midgate

	# sort var using the modified values in var2
	var_sorted = sorted(var, key=var2.__getitem__)
	if (len(var_sorted) == 0):
		return 1

	### vertical sort (if necessary)
	if (hORv == 1):
		sortdata = deepcopy(var_sorted)
		gateslORr = sortdata[0+(midgate*lORr):midgate+(balance*lORr)]
		
		x = 1
		var = deepcopy(G[x])
		var2 = deepcopy(G[1-x])
		keys = deepcopy(var)
		for i in keys:
			if i not in gateslORr:
				del var[i]
				del var2[i]
		for i in var.keys():
			var2[i] = var[i]*100000+var2[i]
		var_sorted = sorted(var, key=var2.__getitem__)
		if (len(var_sorted) == 0):
			return 1	
	core.add_sorted(deepcopy(var_sorted))		
	return 1

def update_coordinates(coordinate, side):	
	## Update the coordinate according to bounding box
	coordinatetemp = coordinate
	# X less than Xmin
	if coordinate[0] < side[0]:
		coordinatetemp[0] = side[0]
	# X more than Xmax
	if coordinate[0] > side[1]:
		coordinatetemp[0] = side[1]
	# Y less than Ymin
	if coordinate[1] < side[2]:
		coordinatetemp[1] = side[2]
	# Y more than Ymax
	if coordinate[1] > side[3]:
		coordinatetemp[1] = side[3]

	return coordinatetemp

def update_location(core, x, y,key):
	xcoord = deepcopy(x)
	ycoord = deepcopy(y)
	gatekeys = deepcopy(key)
	try:
		#print(gatekeys)
		#print(xcoord)
		#print(ycoord)
		if (core.add_location(xcoord, ycoord, gatekeys) is False):
			raise MyError('failed to add location')
	except MyError as e:
		print('My exception occurred, value:', e.value)
	return 1

def containNrun(core, side, hORv, lORr):
	print('Containment ...')
	new = mothercore(side)

	sortgate = core.get_sorted()
	#print(sortgate)
	midgate = math.floor(len(sortgate)/2)
	balance = len(sortgate) - midgate
	# print('sorteddata =', sortgate)

	## new has gates on the given side
	presentgates = []
	for gate in sortgate[0+(midgate*lORr):midgate+(balance*lORr)]:
		presentgates.append(gate)
		new.add_gate(gate, core.get_gateconnections(gate))

	donepads = []
	donegates = []
	donenets = []
	# i in gates on given side
	for i in range(0, len(presentgates)):
		pgatenum = presentgates[i]
		pgateconns = new.get_gateconnections(pgatenum)[1:] # first one is discarded, it has num connections
		# k recurses through nets in each presentgate
		for k in range(0, len(pgateconns)):
			pnet = pgateconns[k]
			if pnet in donenets: continue
			else: donenets.append(pnet)

			netconn = core.get_otherconns(pnet, pgatenum)
			#new.add_netconns(pnet, netconn)
			new.add_net(pnet, pgatenum, 0)

			# process connected pads as new pads
			for l in range(0, len(netconn[1])):
				newconn = netconn[1][l]	# pad number
				new.add_net(pnet, newconn, 1)
				if newconn in donepads: continue
				else: donepads.append(newconn)

				padtemp = [None]*3
				padtemp[0] = pnet # net

				padtemp = core.get_padlocation(newconn)
				if (padtemp[0] != pnet):
					raise MyError('returned pad is not connected to the net')
				else:
					## Correct coordinates which are outside the bounding box
					padtemp[1:3] = update_coordinates(padtemp[1:3], side)
					## Add new pad
					#print('Pads')
					#print(newconn, padtemp)
					new.add_pad(newconn, padtemp)

			# process connected gates as new pads
			for l in range(0, len(netconn[0])):
				newconn = netconn[0][l]
				if newconn in presentgates: 
					new.add_net(pnet, newconn, 0)
					continue
				new.add_net(pnet, newconn+100000, 1)
				if newconn in donegates: continue
				else: donegates.append(newconn)

				padtemp = [None]*3
				padtemp[0] = pnet # net

				#print(pnet, pgatenum, newconn)
				padtemp[1:3] = core.get_gatelocation(newconn)
				## Correct coordinates which are outside the bounding box
				padtemp[1:3] = update_coordinates(padtemp[1:3], side)
				## Correct coordinates which may be inside the bounding box
				if (hORv == 0):
					if (lORr == 0):
						padtemp[1] = side[1]
					else:
						padtemp[1] = side[0]
				else:
					if (lORr == 0):
						padtemp[2] = side[3]
					else:
						padtemp[2] = side[2]
				## Add new pad
				#print('Gates')
				#print(newconn+100000, padtemp)
				new.add_pad(newconn+100000, padtemp) # Offset to remove clash between padnum and gatenum

	#print('gates = ',len(list(new.get_gate().keys())))
	print('Done. Added '+ str(new.get_numG()) + ' gates, '+ str(new.get_numP()) + ' pads, '+ str(new.get_numN()) + ' nets.')
	#print(new.get_nets())
	#print('gates = ',new.get_gate())
	#print('pads = ',new.get_pad())
	#print('nets = ',new.get_nets())

	## if solveforx returns successfully, add all the new locations to core
	if (solveforx(new)):
		xcoord = deepcopy(list(new.gateX.values()))
		ycoord = deepcopy(list(new.gateY.values()))
		gatekeys = deepcopy(list(new.get_gate().keys()))
		update_location(core, xcoord, ycoord, gatekeys)
	else:
		print('Error')
	
	return 1

def place(core,side,n):
	if n >= 8:
		return 1
	else:
		side = deepcopy(side)
		n = n*2
		print()
		print('n = ', n)

		G = deepcopy(core.get_location())
		###
		for l in range(1,core.numG+1):	
			if G[0][l] < side[0]:
				del G[0][l]
				del G[1][l]
			elif G[0][l] > side[1]:
				del G[0][l]
				del G[1][l]
			elif G[1][l] < side[2]:
				del G[1][l]
				del G[0][l]
			elif G[1][l] > side[3]:
				del G[1][l]
				del G[0][l] 
		print('G = ',len(list(G[0].keys())))
		## horizontal sort
		assign(core,G, 0, 0)	
		midx = (side[1]-side[0])/2
			## left containment
		print('#####hor_left')	
		containNrun(core,[side[0],side[1]-midx,side[2],side[3]], 0, 0)
			## right containment
		print('#####hor_right')	
		containNrun(core,[side[0]+midx,side[1],side[2],side[3]], 0, 1)

		## vertical sort
		assign(core,G,1,0)
		midy = (side[3]-side[2])/2
			## left below
		print('#####left_ver_below')		
		containNrun(core,[side[0],side[1]-midx,side[2],side[3]-midy], 1, 0)
				## left above
		print('#####left_ver_above')		
		containNrun(core,[side[0],side[1]-midx,side[2]+midy,side[3]], 1, 1)

		assign(core,G,1,1)
				## right below
		print('#####right_ver_below')
		containNrun(core,[side[0]+midx,side[1],side[2],side[3]-midy], 1, 0)
				## right above
		print('#####right_ver_above')		
		containNrun(core,[side[0]+midx,side[1],side[2]+midy,side[3]], 1, 1)

		## smaller squares
		place(core,[side[0],side[1]-midx,side[2],side[3]-midy],n)	# left below
		place(core,[side[0],side[1]-midx,side[2]+midy,side[3]],n)	# left above
		place(core,[side[0]+midx,side[1],side[2],side[3]-midy],n)	# right below
		place(core,[side[0]+midx,side[1],side[2]+midy,side[3]],n)	# right above

def main():
	filename = sys.argv[1]
	QP1 = create(filename)
	n = 1
	side = [0,100,0,100]
	print()
	print('n = ', n)
	if (solveforx(QP1)):
		place(QP1,side,n)
		writeback(QP1, 'final.txt')

if __name__ == "__main__": main() 
