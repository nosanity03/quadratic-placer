import sys
import matplotlib.pyplot as plt
import itertools

filename = sys.argv[1] 

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

# s is width of dot
plt.autoscale(enable=False, axis="xc")
plt.autoscale(enable=False, axis="yc")
plt.xlim(-10,110)
plt.ylim(-10,110)
plt.scatter(xc,yc, marker='o', s=30, color='red') 

plt.grid()
plt.show()