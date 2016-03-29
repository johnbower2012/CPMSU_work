import time
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import axes3d

pl.axis([0,1000,0,1])
pl.ion()
pl.show()

for i in range(1000):
	y= np.random.random()
	pl.scatter(i,y)
	pl.draw()
	time.sleep(0.05)

data = np.genfromtxt(
	'test.dat',
	usecols = (0,9,10,12),
	names=['time','x','y','r']
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
ax.set_xlabel('time')
ax.set_ylabel('x')
ax.set_zlabel('y')

ax.view_init(elev=20,azim=80)
ax.dist=10

ax.scatter(
	data['time'], data['x'], data['y'],
	color='purple',
	marker='o',
	s=30
)

pl.savefig('mercury.pdf', bbox_inches='tight')
pl.show()

