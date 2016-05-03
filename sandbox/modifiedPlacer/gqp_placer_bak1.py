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
	def __init__(self, numG, numP, numN, side):
		self.gate = dict([])
		self.pad = dict([])
		self.numG = numG
		self.numP = numP
		self.numN = numN
		self.side = side
		self.gateX = dict([])
		self.gateY = dict([])
		self.nets = dict([])
		for l in range(1,numG+1):
			self.gateX[l] = [None]
			self.gateY[l] = [None]
			self.gate[l] = [None]
		for l in range(1,numP+1):
			self.pad[l] = [None]
		for l in range(1,numN+1):
			# self.nets[0] -> gate1
			# self.nets[1] -> gate2, if exists or 0
			# self.nets[2] -> pad, if exists or 0
			self.nets[l] = [0,0,0]
		self.sortedgate = [None]*numG	
	
	def get_gate(self):
		return deepcopy(self.gate)

	def get_pad(self):
		return deepcopy(self.pad)

	def add_gate(self, gatenum, listofconnections):
		self.gate[gatenum] = listofconnections

	def add_pad(self, padnum, listofconnections):
		self.pad[padnum] = listofconnections

	def add_net(self, netnum, connection, gateorpad):
		# 0 for gate, 1 for pad
		if (gateorpad == 1):
			self.nets[netnum][2] = connection
		else:
			if (self.nets[netnum][0] == 0):
				self.nets[netnum][0] = connection
			else:
				self.nets[netnum][1] = connection

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

	def get_location(self):
		l = [None]*2
		l[0] = deepcopy(self.gateX)
		l[1] = deepcopy(self.gateY)
		return l

def create(filename):
	print('Creating data structure from file ...')
	cf = open(filename,'r')

	## Find the number of gates, nets and pads for creating object
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
			break
		else:
			continue
	cf.close()

	## New object
	returncore = mothercore(numG,numP,numN,[0,100,0,100])

	## Add gates, nets and pads to object from file
	cf = open(filename,'r')
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()

		## Add gates and nets
		if i in range(1,numG+1):
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
	print('Done. Added '+ str(numG) + ' gates, '+ str(numP) + ' pads, '+ str(numN) + ' nets.')
	return returncore

def solveforx(core):
	print('Solving for locations ...')
	gate = core.get_gate()
	G = core.numG
	key = list(gate.keys())
	pad = core.get_pad()
	P = core.numP
	if (G == 0):
		return 1

#	Conmat = [0]*G			## Connectivity matrix with cells containing the name of net connected between two gates
	CM = [0]*G
	for i in range(0,G):
		CM[i] = [0]*G

#	for i in range(0,G):
#		if (i%10 == 0): print('0 -> Gate ',i,' of ',G)
#		Conmat[i] = [[] for j in repeat(None, G)]

	net = dict([])			## Dictionary for connected nets for calculating weight of net
	#print(key)
	print('Searching for connected gates ...')
	for i in range(0,G):
		if (i%10 == 0): print('1 -> Gate ',i+1,' of ',G)
		for j in range(0,G):
			if (i == j):
				continue
			else:
				for k in range(1,gate[key[i]][0]+1):
					for l in range(1,gate[key[j]][0]+1):
						#print('i = ',key[i],'j = ',key[j], 'val = ',gate[key[j]][l])
						if (gate[key[i]][k] == gate[key[j]][l]):
							#print('i = ',key[i],'j = ',key[j], 'val = ',gate[key[i]][k])
							val = gate[key[i]][k]
							if val not in net: net[val] = [0,[]]	## 1st 0 to calculate weight, 2nd list for pads
							if key[i] not in net[val]: net[val].append(key[i])
							if key[j] not in net[val]: net[val].append(key[j])
#							if val not in Conmat[i][j]: Conmat[i][j].append(val)

	bx = []
	by = []
	print('Searching for connected pads ...')
	for i in range(0,G):
		if (i%10 == 0): print('2 -> Gate ',i+1,' of ',G)
		for k in range(1,gate[key[i]][0]+1):
			for m in range(1,P+1):
				#print('i = ',key[i],'m = ',m-1)
				#print(pad[m])
				if (gate[key[i]][k] == pad[m][0]):
					#print('i = ',key[i],'m = ',m, 'val = ',gate[key[i]][k])
					val = gate[key[i]][k]
					if val not in net: net[val] = [0,[]]
					if key[i] not in net[val]: net[val].append(key[i])
					if m not in net[val][1]: net[val][1].append(m)
#					if val not in Conmat[i][i]: Conmat[i][i].append(val)

	### Calculating weights
	print('Calculating weights ...')
	for val in net:
		k = len(net[val])-2+len(net[val][1])
		print(val, k)
		print(len(net[val]), len(net[val][1]))
		net[val][0] = 1/(k-1)
	
	#print('net = ', net)
	
	### Conmat
	print('Calculating valid contributions to the cost function ...')

	Conmat = [[] for j in repeat(None, G)]
	for i in range(0,G):
		if (i%10 == 0): print('3 -> Gate ',i+1,' of ',G)
		for j in range(0,G):
			for k in range(1,gate[key[i]][0]+1):
				val = gate[key[i]][k]
				if (i == j):
					for m in range(1,P+1):
						if (gate[key[i]][k] == pad[m][0]):
							if val not in Conmat[i]: Conmat[i].append(val)
				else:
					for l in range(1,gate[key[j]][0]+1):
						#print('i = ',key[i],'j = ',key[j], 'val = ',gate[key[j]][l])
						if (gate[key[i]][k] == gate[key[j]][l]):
							if val not in Conmat[j]: Conmat[j].append(val)
		#print('Conmat = ', Conmat)					
		CMrowsum = 0
		for j in range(0,G):
			CM_temp = 0
			if (i != j):
				while (len(Conmat[j]) > 0):
					val = Conmat[j].pop()
					CM_temp += net[val][0]*1
				CM[i][j] = -CM_temp
			CMrowsum += CM_temp
		CM_temp = 0
		bxtemp = 0
		bytemp = 0
		while (len(Conmat[i]) > 0):
			val = Conmat[i].pop()
			while (len(net[val][1]) > 0):
				valpad = net[val][1].pop()
				#print('i ', i, ' val ',val,'valpad',valpad)
				CM_temp += net[val][0]*1
				bxtemp += net[val][0]*pad[valpad][1]
				bytemp += net[val][0]*pad[valpad][2]
		CM[i][i] = CMrowsum+CM_temp
		if CM[i][i] == 0: 
			print('ZERO ',i)
			CM[i][i] = 0.0001
		bx.append(bxtemp)
		by.append(bytemp)

	del Conmat	
	del net
	#print('CM = ',CM)	
	print('Forming R, C, V matrices ...')
	R = []
	C = []
	V = []
	for i in range(0,G):
		if (i%10 == 0): print('4 -> Gate ',i+1,' of ',G)
		for j in range(0,G):
			if (CM[i][j] != 0) & (math.isnan(CM[i][j]) is False) & (math.isinf(CM[i][j]) is False):
				#print(CM[i][j])
				R.append(i)
				C.append(j)
				V.append(CM[i][j])
	del CM			
	#print('R =', R)
	#print('C =', C)
	#print('V =', V)
	#print('bx =', bx)
	#print('by =', by)
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
	try:
		if (core.add_location(deepcopy(x),deepcopy(y),key) is False):
			raise MyError('failed to add location')
	except MyError as e:
		print('My exception occurred, value:', e.value)
	print('Done.')	
	return 1

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

def containNrun(core, side, hORv, lORr):
	print('Containment ...')
	new = deepcopy(core)
	new.side = deepcopy(side)

	## Update the pad coordinates of existing pads according to bounding box
	for l in range(1,new.numP+1):
		# X less than Xmin
		if new.pad[l][1] < side[0]:
			new.pad[l][1] = side[0]
		# X more than Xmax
		if new.pad[l][1] > side[1]:
			new.pad[l][1] = side[1]
		# Y less than Ymin
		if new.pad[l][2] < side[2]:
			new.pad[l][2] = side[2]
		# Y more than Ymax
		if new.pad[l][2] > side[3]:
			new.pad[l][2] = side[3]
	
	sortgate = new.sortedgate
	print(sortgate)
	midgate = math.floor(len(sortgate)/2)
	balance = len(sortgate) - midgate
	# print('sorteddata =', sortgate)

	# i in leftside or rightside depending on lORr
	for i in range(0+(midgate*lORr),midgate+(balance*lORr)):
		# j in otherside from i depending on lORr
		for j in range(midgate-(midgate*lORr),len(sortgate)-(balance*lORr)):
			# k recurses through all connected nets in gate[sortgate[i]]
			for k in range(1,new.gate[sortgate[i]][0]+1):
				# l recurses through all connected nets in gate[sortgate[j]]
				for l in range(1,new.gate[sortgate[j]][0]+1):
					#print('i = ', sortgate[i],'k = ', k , 'j = ',sortgate[j])

					## If nets match, propagate connected gate to cutline
					if new.gate[sortgate[i]][k] == new.gate[sortgate[j]][l]:
						newtemp = [None]*3
						newtemp[0] = new.gate[sortgate[j]][l]
						newtemp[1] = new.gateX[sortgate[j]]
						newtemp[2] = new.gateY[sortgate[j]]
						
						## Correct coordinates which are outside the bounding box
						# X less than Xmin
						if newtemp[1] < side[0]:
							newtemp[1] = side[0]
						# X more than Xmax	
						if newtemp[1] > side[1]:
							newtemp[1] = side[1]
						# Y less than Ymin        
						if newtemp[2] < side[2]:
							newtemp[2] = side[2]
						# Y more than Ymax        
						if newtemp[2] > side[3]:
							newtemp[2] = side[3]

						## Correct coordinates which may be inside the bounding box
						if (hORv == 0):
							if (lORr == 0):
								newtemp[1] = side[1]
							else:
								newtemp[1] = side[0]
						else:
							if (lORr == 0):
								newtemp[2] = side[3]
							else:
								newtemp[2] = side[2]

						## When gates outside bounding box becomes pads
						if newtemp not in new.pad.values(): 		
							new.numP += 1
							new.pad[new.numP] = newtemp
					else:
						continue

	## remove gates which are on the other side, dictated by lORr
	sortdata = sortgate[0+(midgate*lORr):midgate+(balance*lORr)]
	oldgate = deepcopy(new.gate)
	for gate in oldgate:
		if gate not in sortdata:
			del new.gate[gate]
			new.numG -= 1

	print('gates = ',len(list(new.gate.keys())))
	print('Done. Added '+ str(new.numG) + ' gates, '+ str(new.numP) + ' pads, '+ str(new.numN) + ' nets.')
	#print('pads = ',new.pad)
	
	## if solveforx returns successfully, add all the new locations to core
	if (solveforx(new)):
		try:
			print(deepcopy(list(new.gate.keys())))
			print(deepcopy(list(core.gate.keys())))
			print(deepcopy(list(new.gateX.values())))
			print(deepcopy(list(new.gateY.values())))
			if (core.add_location(deepcopy(list(new.gateX.values())),deepcopy(list(new.gateY.values())),deepcopy(list(core.gate.keys()))) is False):
				raise MyError('failed to add location')
			else:
				del new
		except MyError as e:
			print('My exception occurred, value:', e.value)
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
