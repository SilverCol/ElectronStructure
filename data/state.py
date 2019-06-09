from matplotlib import pyplot as plt
import numpy as np

data = np.fromfile('state.bin')
R = 1000
domain = np.linspace(R/len(data), R + R/len(data), len(data))

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
ax.plot(domain, data)

plt.show()
