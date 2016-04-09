
import matplotlib.pyplot as plt
import itertools

cf = open('final.txt','r')
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
plt.scatter(xc,yc, marker='o', s=30, color='red') 

plt.grid()
plt.show()