from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate

R = 15
state = np.fromfile('state' + str(R) + '.bin')[100::100]
stateH = np.fromfile('stateH' + str(R) + '.bin')[100::100]
domain = np.linspace(R / len(state), R + R / len(state), len(state))

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
ax.set_yscale('log')
He, = ax.plot(domain, state**2)
H, = ax.plot(domain, stateH**2)
ax.legend((He, H), ('$|\\phi_{He}(r)|^2$', '$|\\phi_{H}(r)|^2$'))
ax.set_xlabel('$r$')

plt.show()
