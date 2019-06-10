from matplotlib import pyplot as plt
import numpy as np

R = 15
data = np.fromfile('state' + str(R) + '.bin')
domain = np.linspace(R/len(data), R + R/len(data), len(data))

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
ax.plot(domain, data)

plt.show()
