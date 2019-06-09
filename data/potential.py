from matplotlib import pyplot as plt
import numpy as np

data = np.fromfile('potential.bin')
R = 1000
domain = np.linspace(R/len(data), R + R/len(data), len(data))
v = lambda x : -(x + 1) * np.exp(-2*x) + 1

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
ax.plot(domain, data)
ax.plot(domain, v(domain))

plt.show()
