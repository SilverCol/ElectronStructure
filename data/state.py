from matplotlib import pyplot as plt
import numpy as np
from scipy import integrate

R = 15
state = np.fromfile('state' + str(R) + '.bin')[100::100]
pot = np.fromfile('potential' + str(R) + '.bin')[100::100]
domain = np.linspace(R / len(state), R + R / len(state), len(state))
radial = domain * state
radial /= np.sqrt(integrate.simps(radial**2, domain))

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
# sline, = ax.plot(domain, state**2)
uline, = ax.plot(domain, radial)
pline, = ax.plot(domain, pot)
xline, = ax.plot(domain, (3 / (2 * np.pi**2))**(1/3) * np.power(np.absolute(state), 2/3))
ax.legend((uline, pline, xline), ('$u(r)$', '$U(r)$', '$V_{xc}(r)$'))
ax.set_xlabel('$r$')
ax.set_ylim(-.05, 1.05)
ax.set_title('H atom')

plt.show()
