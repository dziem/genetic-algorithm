from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
import math

'''
cos note
The value passed in this function should be in radians.
https://www.geeksforgeeks.org/python-math-cos-function/
'''
'''
function 1
def z_function(x, y):
	a1 = 0
	a2 = 0
	for i in range(1,6):
		a1 += i * math.cos(math.radians((i + 1) * x + 1))
		a2 += i * math.cos(math.radians((i + 1) * y + 1))
	a1 = a1 * -1
	return a1 * a2
'''
def z_function(x, y):
	a = -1 * math.cos(math.radians(x))
	b = math.cos(math.radians(y))
	c = math.exp((-1 * ((x - math.pi) ** 2)) - ((y - math.pi) ** 2))
	return a * b * c
f2 = np.vectorize(z_function)
print(f2(20,20))
'''
x = np.linspace(-100, 100, 200)
y = np.linspace(-100, 100, 200)

X, Y = np.meshgrid(x, y)
Z = f2(X, Y)
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_wireframe(X, Y, Z, color='green')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
'''