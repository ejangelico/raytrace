import sys
import os
import numpy as np 
import matplotlib.pyplot as plt


zs = np.linspace(0, 100, 100)
phi = 50*np.pi/180.0
th = 10*np.pi/180.0
x = []
y = []
for z in zs:
	x.append(z*np.tan(th)*np.cos(phi))
	y.append(z*np.tan(th*np.sin(phi)))

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(x, y)
ax.set_xlim([-4, 4])
ax.set_ylim([-4, 4])
plt.show()