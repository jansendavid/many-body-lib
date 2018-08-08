import numpy as np
import sys
import matplotlib.pyplot as plt
filename=str(sys.argv[1])
x, y =np.loadtxt(filename, unpack=True)
plt.plot(x, y, 'r.')
plt.show()
plt.savefig("timeplot.png")
