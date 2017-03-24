#!/usr/bin/python
# TODO: add docstring

# test usage
import sys
if len(sys.argv) != 3:
  print("Usage: {} <input file> <limit>".format(sys.argv[0]))
  sys.exit()



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# read in the lines of the log file -- TODO harden into a script parameter
A = np.genfromtxt(sys.argv[1], delimiter=" ")

# compute the relative error for each
realdiff = A[:, 1] - A[:, 3]
imagdiff = A[:, 2] - A[:, 4]
numerator = np.sqrt(realdiff * realdiff + imagdiff * imagdiff)
denom = np.sqrt(A[:, 3] * A[:, 3] + A[:, 4] * A[:, 4])
relerr = numerator / denom


# do a histogram
#bins = 50
#fig, ax = plt.subplots()
# so most of the poins are fine, but there are a handful (relatively speaking
# that have much higher incidental errors)
#n, b, patches = ax.hist(relerr, bins, log=1)
#fig.tight_layout()
#plt.show()

# do the scatter plot
#select = relerr > 0.1
#select2 = relerr < 1.0
#subset = A[np.logical_and(select, select2), :]
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(subset[:,5], subset[:,6], subset[:,7])
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#plt.show()


# select points with relative error too high
select = relerr > float(sys.argv[2])
errs = relerr[select]
subset = A[select, :]
print("The following indices have a relative error greater than {}".format(
      sys.argv[2]))
for i in np.arange(subset.shape[0]):
  print("{} {}".format(int(subset[i,0]), errs[i]))

