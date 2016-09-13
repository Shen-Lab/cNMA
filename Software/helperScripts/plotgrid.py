import matplotlib.pyplot as plt
import numpy as np

# make these smaller to increase the resolution
dx, dy = 2.0, 2.0

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(-3, 3 + dy, dy),
                slice(-3, 3 + dx, dx)]
z = (1 - x / 2. + x ** 5 + y ** 3) * np.exp(-x ** 2 - y ** 2)
# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
z_min, z_max = -np.abs(z).max(), np.abs(z).max()


print x
print y
print z

plt.subplot(2, 2, 1)
plt.pcolormesh(x, y, z)
plt.title('pcolor')
# set the limits of the plot to the limits of the data
#plt.axis([x.min(), x.max(), y.min(), y.max()])
#plt.colorbar()

plt.show()