from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import simps

R = 15
state = np.fromfile('state' + str(R) + '.bin')[100::100]
stateH = np.fromfile('stateH' + str(R) + '.bin')[100::100]
domain = np.linspace(R / len(state), R + R / len(state), len(state))


def f(r, a, b):
    return a * np.exp(-r/b)


# HeOpt, HeCov = curve_fit(f, domain, state)
# print("He | $r = %.2e \\pm %.2e$" % (HeOpt[1], np.sqrt(HeCov[1][1])))
# HOpt, HCov = curve_fit(f, domain, stateH)
# print("H | $r = %.2e \\pm %.2e$" % (HOpt[1], np.sqrt(HCov[1][1])))

rHe = simps(domain**3*state**2, domain)
print("He | r = %.2e" % rHe)
rH = simps(domain**3*stateH**2, domain)
print("H | r = %.2e" % rH)

plt.rcParams.update({'font.size': 15})
fig = plt.figure()
ax = fig.subplots()

ax.grid()
ax.set_yscale('log')
He, = ax.plot(domain, state**2)
# ax.plot(domain, f(domain, HeOpt[0], HeOpt[1])**2, '--', lw=1)
H, = ax.plot(domain, stateH**2)
# ax.plot(domain, f(domain, HOpt[0], HOpt[1])**2, '--', lw=1)
ax.legend((He, H), ('$|\\phi_{He}(r)|^2$', '$|\\phi_{H}(r)|^2$'))
ax.set_xlabel('$r$')

plt.show()
