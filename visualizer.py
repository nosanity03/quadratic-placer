import sys
import matplotlib.pyplot as plt
import itertools

filename = sys.argv[1] 

if (len(sys.argv) > 2):
	filename2 = sys.argv[2]

cf = open(filename,'r')
xc = list()
yc = list()
for i, line in enumerate(cf):
	line = line.strip()
	if line == '': continue
	v = line.split()
	xc.append(float(v[1]))
	yc.append(float(v[2]))

cf.close()

if (len(sys.argv) > 2):
	cf = open(filename2,'r')
	xp = list()
	yp = list()
	for i, line in enumerate(cf):
		line = line.strip()
		if line == '': continue
		v = line.split()
		xp.append(float(v[1]))
		yp.append(float(v[2]))
	cf.close()

# s is width of dot
plt.autoscale(enable=False, axis="xc")
plt.autoscale(enable=False, axis="yc")
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.scatter(xc,yc, marker='o', s=20, color='red') 

if (len(sys.argv) > 2):
	plt.scatter(xp,yp, marker='s', s=30, color='blue') 

plt.grid()
plt.show()
