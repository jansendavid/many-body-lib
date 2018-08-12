import numpy as np
import sys
import matplotlib.pyplot as plt
filename=str(sys.argv[1])
filename2=str(sys.argv[2])
x, y =np.loadtxt(filename, unpack=True)
z, w =np.loadtxt(filename2, unpack=True)
plt.plot(x, y, 'r.')
plt.plot(z, w, 'b.')
plt.show()
plt.savefig("timeplot.png")
