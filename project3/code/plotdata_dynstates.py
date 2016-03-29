#test.py
#!/usr/bin/python


import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d

data = np.genfromtxt(
	'dynamicstates.dat',
	usecols = (0,1,2,3,4),
	names=['time','x','y','z','r']
)


fig1 = pl.figure()

fig1.suptitle('Mercury')

ax = fig1.add_subplot(211)
plot1 = ax.plot(data['time'],data['r'], 'k-')
ax.grid(True)
ax.set_xlabel('time')
ax.set_ylabel('r')

ax = fig1.add_subplot(2,1,2,projection='3d')
#ax = fig1.gca(projection='3d')
#ax.set_title('Mercury')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.view_init(elev=45,azim=340)
ax.dist=10

ax.scatter(
	data['x'], data['y'], data['z'],
	color='purple',
	marker='o',
	s=30
)

pl.savefig('mercury.pdf', bbox_inches='tight')
pl.show()

