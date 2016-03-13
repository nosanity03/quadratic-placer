#! /usr/bin/python3
from copy import deepcopy
import numpy as np
from numpy import array
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import math

class MyError(Exception):
	def __init__(self, value):
	        self.value = value

class mothercore:
	def __init__(self, gate, pad, numG, numP, numN, side):
		self.gate = gate
		self.pad = pad
		self.numG = numG
		self.numP = numP
		self.numN = numN
		self.side = side
		self.gateX = dict([])
		self.gateY = dict([])
		for l in range(1,numG+1):
			self.gateX[l] = [None]
			self.gateY[l] = [None]
		self.sortedgate = [None]*numG	
	
	def get_gate(self):
		return deepcopy(self.gate)

	def get_pad(self):
		return deepcopy(self.pad)

	def add_location(self, x, y, data):	# data is required to know the keys of x,y values
		if (len(data) == len(x)):
			for l in range(0,len(data)):
				self.gateX[data[l]] = x[l]
				self.gateY[data[l]] = y[l]
			return 1
		else:
			return 0

	def add_sorted(self, sort):
		self.sortedgate = sort
		return 1

	def get_location(self):
		l = [None]*2
		l[0] = deepcopy(self.gateX)
		l[1] = deepcopy(self.gateY)
		return l

def create(filename):
	cf = open(filename,'r')
	gate = dict([])
	pad = dict([])
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
		else:
			if i in range(1,numG+1):
				gate[int(v[0])] = [int(v[j]) for j in range(1,int(v[1])+2)]
			elif i in range(numG+2,numG+numP+2):
				pad[int(v[0])] = [int(v[j]) for j in range(1,4)]

	cf.close()
	return mothercore(gate,pad,numG,numP,numN,[0,100,0,100])

def solveforx(core):
	gate = core.get_gate()
	G = core.numG
	key = list(gate.keys())
	pad = core.get_pad()
	P = core.numP
	Conmat = [0]*G			## Connectivity matrix with cells containing the name of net connected between two gates
	CM = [0]*G
	net = dict([])			## Dictionary for connected nets for calculating weight of net
	#print(key)
	for i in range(0,G):
		Conmat[i] = [0]*G
		CM[i] = [0]*G
		for j in range(0,G):
			Conmat[i][j] = []
			if (i == j):
				continue
			else:
				for k in range(1,gate[key[i]][0]+1):
					for l in range(1,gate[key[j]][0]+1):
						#print('i = ',key[i],'j = ',key[j], 'val = ',gate[key[j]][l])
						if (gate[key[i]][k] == gate[key[j]][l]):
							val = gate[key[i]][k]
							if val not in net: net[val] = [0,[]]	## 1st 0 to calculate weight, 2nd list for pads
							if key[i] not in net[val]: net[val].append(key[i])
							if key[j] not in net[val]: net[val].append(key[j])
							if val not in Conmat[i][j]: Conmat[i][j].append(val)
	bx = []
	by = []
	for i in range(0,G):
		for k in range(1,gate[key[i]][0]+1):
			for m in range(1,P+1):
				#print('i = ',key[i],'m = ',m)
				if (gate[key[i]][k] == pad[m][0]):
					val = gate[key[i]][k]
					if val not in net: net[val] = [0,[]]
					if key[i] not in net[val]: net[val].append(key[i])
					if m not in net[val][1]: net[val][1].append(m)
					if val not in Conmat[i][i]: Conmat[i][i].append(val)

	### Calculating weights
	for val in net:
		k = len(net[val])-2+len(net[val][1])
		net[val][0] = 1/(k-1)
	
	#print('net = ', net)
	#print('Conmat = ', Conmat)
	
	### Conmat
	for i in range(0,G):
		CMrowsum = 0
		for j in range(0,G):
			CM_temp = 0
			if (i != j):
				for val in Conmat[i][j]:
					CM_temp += net[val][0]*1
				CM[i][j] = -CM_temp
			CMrowsum += CM_temp
		CM_temp = 0
		bxtemp = 0
		bytemp = 0
		for val in Conmat[i][i]:
			for valpad in net[val][1]:
				#print('i ', i, ' val ',val,'valpad',valpad)
				CM_temp += net[val][0]*1
				bxtemp += net[val][0]*pad[valpad][1]
				bytemp += net[val][0]*pad[valpad][2]
		CM[i][i] = CMrowsum+CM_temp
		bx.append(bxtemp)
		by.append(bytemp)
	#print(CM)	
	R = []
	C = []
	V = []
	for i in range(0,G):
		for j in range(0,G):
			if CM[i][j] != 0:
				R.append(i)
				C.append(j)
				V.append(CM[i][j])
	#print('R =', R)
	#print('C =', C)
	#print('V =', V)
	#print('bx =', bx)
	#print('by =', by)
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
	print('x = ', x)
	print('y = ', y)
	try:
		if (core.add_location(deepcopy(x),deepcopy(y),key) is False):
			raise MyError('failed to add location')
	except MyError as e:
		print('My exception occurred, value:', e.value)
	return 1

def writeback(core,filename):
	G = core.get_location()
	wf = open(filename,'w')
	for l in range(1,core.numG+1):
		print(l," ",G[0][l]," ", G[1][l], file=wf)
	wf.close()

def assign(core, hORv):
	G = core.get_location()
	var = G[hORv]		# 0 for horizontal (x) and 1 for vertical (y)
	var_temp = sorted(var.values())
	midgate = math.floor(core.numG/2)
	try:
		#print(midgate)
		got = 0
		var_sorted = sorted(var, key=var.__getitem__)
		if (var_temp[midgate-1] == var_temp[midgate]):
			raise MyError('equal var')
			var2 = G[1-hORv]
			var2_temp = sorted(var2, key=var2.__getitem__)
			lessx = var_sorted[midgate-1]
			morex = var_sorted[midgate]
			for i in range(0,core.numG):
				if var2_temp[i] == lessx:
					break
				elif var2_temp[i] == morex:
					got = 1
					break
		if got == 1:
			var_sorted[midgate-1] = morex
			var_sorted[midgate] = lessx
		core.add_sorted(deepcopy(var_sorted))
	except MyError as e:
		print('My exception occurred, value:', e.value)
	return 1

def containNrun(core, side, hORv, lORr):
	new = deepcopy(core)
	new.side = side
	for l in range(1,new.numP+1):
		if new.pad[l][1] < side[0]:
			new.pad[l][1] = side[0]
		elif new.pad[l][1] > side[1]:
			new.pad[l][1] = side[1]
		elif new.pad[l][2] < side[2]:
			new.pad[l][2] = side[2]
		elif new.pad[l][2] > side[3]:
			new.pad[l][2] = side[3]
	
	sortgate = new.sortedgate
	midgate = math.floor(new.numG/2)
	balance = new.numG - midgate
	#print('sorteddata =', sortgate)
	for i in range(0+(midgate*lORr),midgate+(balance*lORr)):
		for j in range(midgate-(midgate*lORr),new.numG-(balance*lORr)):
			for k in range(1,new.gate[sortgate[i]][0]+1):
				for l in range(1,new.gate[sortgate[j]][0]+1):
					#print('i = ', sortgate[i],'k = ', k , 'j = ',sortgate[j])
					if new.gate[sortgate[i]][k] == new.gate[sortgate[j]][l]:
						newtemp = [None]*3
						newtemp[0] = new.gate[sortgate[j]][l]
						newtemp[1] = new.gateX[sortgate[j]]
						newtemp[2] = new.gateY[sortgate[j]]
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
						if newtemp not in new.pad.values(): 		
							new.numP += 1
							new.pad[new.numP] = newtemp
					else:
						continue
	g = deepcopy(new.numG)
	for j in range(midgate-(midgate*lORr),g-(balance*lORr)):
		del new.gate[sortgate[j]]	
		new.numG -= 1
	#print()
	#print('gates = ',list(new.gate.keys()))
	#print('pads = ',new.pad)
	if (solveforx(new)):
		try:
			if (core.add_location(deepcopy(list(new.gateX.values())),deepcopy(list(new.gateY.values())),deepcopy(list(core.gate.keys()))) is False):
				raise MyError('failed to add location')
		except MyError as e:
			print('My exception occurred, value:', e.value)
	else:
		print('Error')
	
	return 1

def main():
	QP1 = create('benchmarks/3QP/toy2')
	if (solveforx(QP1)):
		writeback(QP1, 'QP1.txt')
		print()
		assign(QP1,0) 				# 0 for horizontal sort, 1 for vertical sort
		containNrun(QP1,[0,50,0,100], 0, 0)	# right 0 for left containment (horizontal or vertical) 
		writeback(QP1, 'QP2.txt')
		print()
		containNrun(QP1,[50,100,0,100], 0, 1)   # right 1 for right containment (horizontal or vertical)
		writeback(QP1, 'QP3.txt')

if __name__ == "__main__": main() 
