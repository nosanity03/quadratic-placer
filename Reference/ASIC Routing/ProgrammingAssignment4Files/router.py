#! /usr/bin/python3
from copy import deepcopy

class route:
	def __init__(self, grid, net, numX, numY, bend_pen, via_pen):
		self.grid = grid
		self.net = net
		self.numX = numX
		self.numY = numY
		self.bend = bend_pen
		self.via = via_pen
	
	def get_grid(self):
		return deepcopy(self.grid)

	def get_net(self):
		return deepcopy(self.net)

def create(filename):
	print('Creating data structure from file ...')
	cf = open(filename+'.grid','r')
	grid = dict([])
	k = 1
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()
		if i == 0:
			numX = int(v[0])
			numY = int(v[1])
			bend_pen = int(v[2])
			via_pen = int(v[3])
			grid[k] = [0]*numY
		else:
			if i == k*numY+1:	
				k += 1
				grid[k] = [0]*numY
			grid[k][i-1-((k-1)*numY)] = [int(v[j]) for j in range(0,numX)]

	cf.close()
	#print('grid = ',grid)

	cf = open(filename+'.nl','r')
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()
		if i == 0:
			net = [0]*int(v[0])
		else:
			net[i-1] = (int(v[1]),int(v[2]),int(v[3]),int(v[4]),int(v[5]),int(v[6]))
	cf.close()
	#print('net = ',net)
	return route(grid,net,numX,numY,bend_pen,via_pen)

def iterator(ng,grid,Src,T,bend,via):
	tags = ('N','S','E','W','U','D')
	Src_ret = []
	if len(Src) < 1:
		return 0

	for S in Src:
		if (S == T):
			continue
		else:
			l = S[0]
			x = S[1]
			y = S[2]
			val = ng[l][y][x][0]
			S_pred = ng[l][y][x][1]

			## Moving towards south
			ch = 0
			tx = x
			ty = y-1
			if (ty >= 0):
				tval = grid[l][ty][tx]
				fval = val+tval
				if (S_pred == tags[1]) | (S_pred == tags[2]) | (S_pred == tags[3]): fval += bend
				if (tval != -1):
					if (type(ng[l][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[l][ty][tx][0]): ch = 1
						if (fval == ng[l][ty][tx][0]):
							#delx = T[1] - tx
							dely = T[2] - ty
							if (dely < 0): ch = 1
			if (ch == 1):
				ng[l][ty][tx] = (fval,tags[0])
			        #print(fval)
				if (l,tx,ty) not in Src_ret: Src_ret.append((l,tx,ty))

			## Moving towards north
			ch = 0
			tx = x
			ty = y+1
			if (ty < len(grid[l])):
				tval = grid[l][ty][tx]
				fval = val+tval
				if (S_pred == tags[0]) | (S_pred == tags[2]) | (S_pred == tags[3]): fval += bend
				if (tval != -1):
					if (type(ng[l][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[l][ty][tx][0]): ch = 1
						if (fval == ng[l][ty][tx][0]):
							#delx = T[1] - tx
							dely = T[2] - ty
							if (dely > 0): ch = 1
			if (ch == 1):
				ng[l][ty][tx] = (fval,tags[1])
				#print(fval)
				if (l,tx,ty) not in Src_ret: Src_ret.append((l,tx,ty))

			## Moving towards west
			ch = 0
			tx = x-1
			ty = y
			if (tx >= 0):
				tval = grid[l][ty][tx]
				fval = val+tval
				if (S_pred == tags[0]) | (S_pred == tags[1]) | (S_pred == tags[3]): fval += bend
				if (tval != -1):
					if (type(ng[l][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[l][ty][tx][0]): ch = 1
			if (ch == 1):
				ng[l][ty][tx] = (fval,tags[2])
				#print(fval)
				if (l,tx,ty) not in Src_ret: Src_ret.append((l,tx,ty))
	
			## Moving towards east
			ch = 0
			tx = x+1
			ty = y
			if (tx < len(grid[l][ty])):
				tval = grid[l][ty][tx]
				fval = val+tval
				if (S_pred == tags[0]) | (S_pred == tags[1]) | (S_pred == tags[2]): fval += bend
				if (tval != -1): 
					if (type(ng[l][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[l][ty][tx][0]): ch = 1
			if (ch == 1):				
				ng[l][ty][tx] = (fval,tags[3])
				#print(fval)
				if (l,tx,ty) not in Src_ret: Src_ret.append((l,tx,ty))

			## Moving up
			ch = 0
			tx = x
			ty = y
			tl = l+1
			if (tl <= 2):
				tval = grid[tl][ty][tx]
				fval = val+tval+via
				if (tval != -1):
					if (type(ng[tl][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[tl][ty][tx][0]): ch = 1
			if (ch == 1):
				ng[tl][ty][tx] = (fval,tags[5])
				#print(fval)
				if (tl,tx,ty) not in Src_ret: Src_ret.append((tl,tx,ty))

			## Moving down
			ch = 0
			tx = x
			ty = y
			tl = l-1
			if (tl > 0):
				tval = grid[tl][ty][tx]
				fval = val+tval+via
				if (tval != -1):
					if (type(ng[tl][ty][tx]) == int):
						ch = 1
					else:
						if (fval < ng[tl][ty][tx][0]): ch = 1
			if (ch == 1):
				ng[tl][ty][tx] = (fval,tags[4])
				#print(fval)
				if (tl,tx,ty) not in Src_ret: Src_ret.append((tl,tx,ty))

	#print('Src : ' , Src_ret)
	iterator(ng,grid,Src_ret,T,bend,via)

def traceback(pred,loc):
	l = loc[0]
	x = loc[1]
	y = loc[2]
	if pred == 'N':
		y += 1
	elif pred == 'S':
		y -= 1
	elif pred == 'E':
		x += 1
	elif pred == 'W':
		x -= 1
	elif pred == 'U':
		l += 1
	elif pred == 'D':
		l -= 1
	return (l,x,y)

def wavefront(grid,net,bend,via):
	S = (net[0],net[1],net[2])
	T = (net[3],net[4],net[5])
	print('S = ',S,'T = ',T)
	
	ng = deepcopy(grid)
	ng[S[0]][S[2]][S[1]] = (1,'SRC') ## (Cost, Predecessor)
	Src = [S]
	iterator(ng,grid,Src,T,bend,via)
	#print(ng)

	## Backtrace
	path = [T]

	l = T[0]
	x = T[1]
	y = T[2]
	try:
		present = ng[l][y][x][1]
		print('Path Cost = ', ng[l][y][x][0])
		cell = (l,x,y)
		while (present != 'SRC'):
			cell = traceback(present,cell)
			path.append(cell)
			present = ng[cell[0]][cell[2]][cell[1]][1]
	except:
		print('Unroutable!!!')
		path.append(S)
	return path

def router(problem,filename):
	A = deepcopy(problem)
	tags = ('N','S','E','W','U','D')
	grid_temp = A.get_grid()
	wf = open(filename+'.txt','w')

	nets = A.net
	print(len(nets),file=wf)
	for i, net in enumerate(nets):
		print()
		print('Net # ',i+1)
		grid_temp[net[0]][net[2]][net[1]] = 1
		grid_temp[net[3]][net[5]][net[4]] = 1
		path = wavefront(grid_temp,net,A.bend,A.via)
		#print('Path = ',path)
		print(i+1,file=wf)
		while len(path) > 0:
			cell = path.pop()
			grid_temp[cell[0]][cell[2]][cell[1]] = -1
			print(cell[0],cell[1],cell[2], file=wf)
		print(0,file=wf)	

	wf.close()
	return 1

def main():
	filename = 'bench4'
	QP1 = create(filename)
	router(QP1,filename)

if __name__ == "__main__": main() 
