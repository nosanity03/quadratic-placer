#! /usr/bin/python3

## @package gqp_placer
# General Quadratic Placer

import sys

from copy import deepcopy
import numpy as np
from numpy import array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import math
from itertools import repeat

class MyError(Exception):
	## The constructor.
	def __init__(self, value):
		self.value = value

class mothercore:
	## The constructor.
	def __init__(self, side):
		## Dictionary of gates
		self.gate = dict([])	
		## Dictionary of pads
		self.pad = dict([])		
		## Number of gates
		self.numG = 0		
		## Number of pads	
		self.numP = 0			
		## Number of nets
		self.numN = 0			
		## List of 4 elements storing grid limits [xmin, xmax, ymin, ymax]
		self.side = side 		
		## Dictionary of X coordinates
		self.gateX = dict([])	
		## Dictionary of Y coordinates
		self.gateY = dict([])	
		## Dictionary of nets
		self.nets = dict([])	
		## List of sorted gates
		self.sortedgate = []

	## Helper function to return the number of gates.
	def get_numG(self):
		return deepcopy(self.numG)
	
	## Helper function to return the number of pads
	def get_numP(self):
		return deepcopy(self.numP)

	## Helper function to return the number of nets
	def get_numN(self):
		return deepcopy(self.numN)

	## Helper function to set the side variable
	def set_side(self, side):
		self.side = deepcopy(side)
	
	## Helper function to return a dictionary of gates
	def get_gate(self):
		return deepcopy(self.gate)

	## Helper function which returns the connections of given gate
	def get_gateconnections(self, gatenum):
		return deepcopy(self.gate[gatenum])

	## Helper function which returns the location of a given gate of number 'gatenum'
	def get_gatelocation(self, gatenum):
		return deepcopy([self.gateX[gatenum], self.gateY[gatenum]])

	## Helper function which returns the dictionary of pads
	def get_pad(self):
		return deepcopy(self.pad)

	## Helper function which returns the location of a given pad of number 'padnum'
	def get_padlocation(self, padnum):
		return deepcopy(self.pad[padnum])

	## Helper function which makes a new gate and adds list of connections
	def add_gate(self, gatenum, listofconnections):
		self.gate[gatenum] = deepcopy(listofconnections)
		self.numG += 1
		return 1

	## Helper function which makes a new pad and adds its connections and location
	def add_pad(self, padnum, netandlocation):
		self.pad[padnum] = deepcopy(netandlocation)
		self.numP += 1
		return 1

	## Helper function which makes a new net, if needed, and appends connections to the net 'netnum'
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

	## Helper function which adds a list of connections to a net of number 'netnum'
	def add_netconns(self, netnum, listofconnections):
		# 0 for gate, 1 for pad

		## if netnum doesn't already exist in dictionary
		if netnum not in self.nets:
			self.nets[netnum] = [[],[]] # first list for gates, second for pads
			self.numN += 1

		self.nets[netnum] = deepcopy(listofconnections)
		return 1

	## Helper function which returns the connections of a net 'netnum' other than the given gate 
	def get_otherconns(self, netnum, gatenum):
		# returns the other entities connected to Net "netnum", other than "gatenum".
		# returns a list of 2 lists: first list for gates, second for pads
		netconns = deepcopy(self.nets[netnum])

		# gatenum is not in netconnections
		if (gatenum not in netconns[0]): return [[],[]]

		# gatenum is in netconnections
		netconns[0].remove(gatenum)
		return netconns

	## Helper function which returns the dictionary of nets
	def get_nets(self):
		return deepcopy(self.nets)

	## Helper function which adds location values for given gate keys
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

	## Helper function to add a sorted list of gates
	def add_sorted(self, sort):
		self.sortedgate = deepcopy(sort)
		return 1

	## Helper function which returns a sorted list of gates
	def get_sorted(self):
		return deepcopy(self.sortedgate)

	## Helper function which returns a list of locations
	def get_location(self):
		l = [None]*2
		l[0] = deepcopy(self.gateX)
		l[1] = deepcopy(self.gateY)
		return l

## Creates the "mothercore" class from file
def create(filename):
	print('Creating data structure from file ...')

	###########################################################################################
	## Open file
	cf = open(filename,'r')

	###########################################################################################
	## New object
	returncore = mothercore([0,100,0,100])

	###########################################################################################
	## Add gates, nets and pads to object from file
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()

		if i == 0:
			## Number of gates and nets
			numG = int(v[0])
			numN = int(v[1])
			numP = 1
		elif i == numG+1:
			## Number of pads
			numP = int(v[0])
		elif i in range(1,numG+1):
			## Add gates and nets
			nets = [int(v[1])] # number of nets
			for j in range(2,int(v[1])+2):
				nets.append(int(v[j]))
				returncore.add_net(int(v[j]), int(v[0]), 0)
			returncore.add_gate(int(v[0]), nets)
		elif i in range(numG+2,numG+numP+2):
			## Add pads and nets
			nets = []
			for j in range(1,4):
				nets.append(int(v[j]))
			returncore.add_net(int(v[1]), int(v[0]), 1)
			returncore.add_pad(int(v[0]), nets)
		else:
			continue

	###########################################################################################
	## Close the file
	cf.close()

	###########################################################################################
	## Return
	print('Done. Added '+ str(returncore.get_numG()) + ' gates, '+ str(returncore.get_numP()) + ' pads, '+ str(returncore.get_numN()) + ' nets.')
	return returncore
	###########################################################################################

## Solves the given core for locations of gates
def solveforx(core):
	print('Solving for locations ...')

	###########################################################################################
	## Initializations
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
	## Calculating weights using number of connections of nets
	print('Calculating weights ...')
	for val in nets:
		k = len(nets[val][0])+len(nets[val][1])
		#print(val, k)
		#print(len(nets[val][0]), len(nets[val][1]))
		weights[val] = 1/(k-1)
	
	#print('nets = ', nets)
	#print('weights = ', weights)
	
	###########################################################################################
	## Gate numbers may vary, this dictionary keeps them in order
	gateorder = dict([])
	for i in range(0,G):
		gateorder[key[i]] = i

	###########################################################################################
	## Calculating C and A matrices, bx and by vectors
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
	## Derive R, C, V matrices for sparse matrix generation
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
	## Solve for x and y vectors
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
	## Update new location values in the class
	update_location(core, x, y,key)
	print('Done.')	
	return 1
	###########################################################################################

## Read and write back to given filename
def writeback(core,filename):
	G = core.get_location()
	wf = open(filename,'w')
	for l in range(1,core.numG+1):
		print(l," ",G[0][l]," ", G[1][l], file=wf)
	wf.close()
	###########################################################################################

## Sorts the gates from the given core according to hORv and lORr
# 1 horizontal sort, then 1 vertical sort if needed
def assign(core, G, hORv, lORr):	
	print('Assignment ...')

	###########################################################################################
	## horizontal sort
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

	###########################################################################################
	## vertical sort (if necessary)
	if (hORv == 1):

		# vertical sort after horizontal depends on side -> left or right
		sortdata = deepcopy(var_sorted)
		gateslORr = sortdata[0+(midgate*lORr):midgate+(balance*lORr)]
		
		x = 1
		var = deepcopy(G[x])
		var2 = deepcopy(G[1-x])
		keys = deepcopy(var)

		# remove gates which belong to the other side -> right or left
		for i in keys:
			if i not in gateslORr:
				del var[i]
				del var2[i]

		# merge 2 sort keys into 1 using (100000*y+x). This will work as long
		# as there are less than 100k gates.
		for i in var.keys():
			var2[i] = var[i]*100000+var2[i]

		# sort var using the modified values in var2
		var_sorted = sorted(var, key=var2.__getitem__)
		if (len(var_sorted) == 0):
			return 1	

	###########################################################################################
	## Add the sorted gates back to class
	core.add_sorted(deepcopy(var_sorted))		
	return 1
	###########################################################################################

## Updates the coordinate according to bounding box
def update_coordinates(coordinate, side):
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
	###########################################################################################

## Updates the location values given in the class
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
	###########################################################################################

## Contains half the gates in core within the bounding box and solves for gate locations
def containNrun(core, side, hORv, lORr):
	print('Containment ...')

	###########################################################################################
	## create a new core
	new = mothercore(side)

	###########################################################################################
	## add gates to the new core, which are on appropriate side
	sortgate = core.get_sorted()
	#print(sortgate)
	midgate = math.floor(len(sortgate)/2)
	balance = len(sortgate) - midgate
	# print('sorteddata =', sortgate)

	# new has gates on the given side
	presentgates = []
	for gate in sortgate[0+(midgate*lORr):midgate+(balance*lORr)]:
		presentgates.append(gate)
		new.add_gate(gate, core.get_gateconnections(gate))

	###########################################################################################	
	## add pads and nets to the new core
	donepads = []
	donegates = []
	donenets = []
	# i recurses through the gates on given side
	for i in range(0, len(presentgates)):
		pgatenum = presentgates[i]
		pgateconns = new.get_gateconnections(pgatenum)[1:] # first one is discarded, it has num connections
		# k recurses through the nets in each presentgate
		for k in range(0, len(pgateconns)):
			pnet = pgateconns[k]
			if pnet in donenets: continue
			else: donenets.append(pnet)

			# for each net get all the other connections apart from the presentgate
			netconn = core.get_otherconns(pnet, pgatenum)
			new.add_net(pnet, pgatenum, 0)

			# process connected pads as new pads
			for l in range(0, len(netconn[1])):
				newconn = netconn[1][l]	# pad number
				new.add_net(pnet, newconn, 1)

				# skipping if newconn has been processed, else go ahead.
				if newconn in donepads: continue
				else: donepads.append(newconn)

				padtemp = [None]*3
				padtemp[0] = pnet # net

				padtemp = core.get_padlocation(newconn)
				if (padtemp[0] != pnet):
					raise MyError('returned pad is not connected to the net')
				else:
					# correct coordinates which are outside the bounding box
					padtemp[1:3] = update_coordinates(padtemp[1:3], side)
					# add new pad
					new.add_pad(newconn, padtemp)

			# process connected gates as new pads
			for l in range(0, len(netconn[0])):
				newconn = netconn[0][l] # gate number

				# skipping if newconn is a gate inside the bounding box
				if newconn in presentgates: 
					new.add_net(pnet, newconn, 0)
					continue
				new.add_net(pnet, newconn+100000, 1)

				# skipping if newconn has been processed, else go ahead.
				if newconn in donegates: continue
				else: donegates.append(newconn)

				padtemp = [None]*3
				padtemp[0] = pnet # net
				padtemp[1:3] = core.get_gatelocation(newconn)

				# correct coordinates of gates if they are outside the bounding box
				padtemp[1:3] = update_coordinates(padtemp[1:3], side)
				
				# correct coordinates of gates if they are inside the bounding box
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

				# add new pad
				new.add_pad(newconn+100000, padtemp) # Offset to remove clash between padnum and gatenum

	print('Done. Added '+ str(new.get_numG()) + ' gates, '+ str(new.get_numP()) + ' pads, '+ str(new.get_numN()) + ' nets.')
	#print(new.get_nets())
	#print('gates = ',new.get_gate())
	#print('pads = ',new.get_pad())
	#print('nets = ',new.get_nets())

	###########################################################################################	
	## solve "new" for locations of gates inside "new" and add them to the original core
	if (solveforx(new)):
		# if solveforx returns successfully, add all the new locations to core
		xcoord = deepcopy(list(new.gateX.values()))
		ycoord = deepcopy(list(new.gateY.values()))
		gatekeys = deepcopy(list(new.get_gate().keys()))
		update_location(core, xcoord, ycoord, gatekeys)
	else:
		print('Error')
	
	return 1
	###########################################################################################	

## Iterative Placer
def place(core,side,n):
	###########################################################################################	
	## stopping condition
	if n >= 8:
		return 1

	###########################################################################################
	## if n is within bounds
	else:
		side = deepcopy(side)
		n = n*2
		print()
		print('n = ', n)
		G = deepcopy(core.get_location())
		
		# remove the gate locations not within the bounding box
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
		
		# sort overall horizontally
		assign(core,G, 0, 0)	
		midx = (side[1]-side[0])/2
		
		# left containment
		print('#####hor_left')	
		containNrun(core,[side[0],side[1]-midx,side[2],side[3]], 0, 0)
		
		# right containment
		print('#####hor_right')	
		containNrun(core,[side[0]+midx,side[1],side[2],side[3]], 0, 1)

		# sort left half vertically
		assign(core,G,1,0)
		midy = (side[3]-side[2])/2
		
		# left below
		print('#####left_ver_below')		
		containNrun(core,[side[0],side[1]-midx,side[2],side[3]-midy], 1, 0)
		
		# left above
		print('#####left_ver_above')		
		containNrun(core,[side[0],side[1]-midx,side[2]+midy,side[3]], 1, 1)

		# sort right half vertically
		assign(core,G,1,1)
		
		# right below
		print('#####right_ver_below')
		containNrun(core,[side[0]+midx,side[1],side[2],side[3]-midy], 1, 0)
		
		# right above
		print('#####right_ver_above')		
		containNrun(core,[side[0]+midx,side[1],side[2]+midy,side[3]], 1, 1)

		# smaller squares
		place(core,[side[0],side[1]-midx,side[2],side[3]-midy],n)	# left below
		place(core,[side[0],side[1]-midx,side[2]+midy,side[3]],n)	# left above
		place(core,[side[0]+midx,side[1],side[2],side[3]-midy],n)	# right below
		place(core,[side[0]+midx,side[1],side[2]+midy,side[3]],n)	# right above

	###########################################################################################

## Main
def main():
	## create mothercore class
	filename = sys.argv[1]
	QP1 = create(filename)
	n = 1
	side = [0,100,0,100]
	print()
	print('n = ', n)
	if (solveforx(QP1)):
		## recursively place
		place(QP1,side,n)
		## write final placement to a file
		writeback(QP1, 'final.txt')

if __name__ == "__main__": main() 
