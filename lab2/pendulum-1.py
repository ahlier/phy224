import numpy as np
from matplotlib import pyplot as plt
import scipy.optimize as op

uncertainty = np.loadtxt('uncertainty.txt', skiprows=1, unpack=True)
period_1, period_2, position = np.loadtxt('data.txt', skiprows=1, unpack=True)

length = 1.00131  # In meter
length_uncer = 0.0015   # Three measurements were taken, and uncertainty for
# each measurement was 0.0005m

uncer = np.ones(len(period_1))*np.std(uncertainty)
# uncertainty for period of 16 oscillation


def model_func(x, a, b):
    # the relation between fine adjustment mass position and period is
    # approximated linear when the mass displacement is small
    return a*x+b


popt1, pcov1 = op.curve_fit(model_func, position, period_1, absolute_sigma=True, sigma=uncer)
popt2, pcov2 = op.curve_fit(model_func, position, period_2, absolute_sigma=True, sigma=uncer)
graph = plt.plot(position, period_1, 'ro', label='Period 1')
graph2 = plt.plot(position, period_2, 'bo', label='Period 2')
graph3 = plt.plot(position, model_func(position, *popt1), 'r-')
graph4 = plt.plot(position, model_func(position, *popt2), 'b-')

plt.xlabel('Position of mass (cm)')
plt.ylabel('Period of 16 oscillations (s)')
plt.title('Change in Period of Oscillation With Respect to \n'
          'Change in Fine Adjustment Mass Position')
plt.legend()
plt.show()

root = (popt2[1]-popt1[1])/(popt1[0]-popt2[0])
# this calculates where these two trend lines intersect
period = model_func(root, *popt2)/16    # period of one oscillation when the
# length between pivots points are effective length
period_uncer = (abs(pcov1[0][0]/popt1[0])+abs(pcov1[1][1]/popt1[1])+
              abs(pcov2[0][0]/popt2[0])+abs(pcov2[1][1]/popt2[1]))*period
# uncertainty of period

g = (2*np.pi)**2*(length/period**2)
g_uncer = g*((period_uncer/period)*2+length_uncer/length)  # uncertainty of g
chi_square1 = (1/(len(position)-2))*sum(((period_1-model_func(position, *popt1))/uncer)**2)
chi_square2 = (1/(len(position)-2))*sum(((period_2-model_func(position, *popt2))/uncer)**2)
# chi-square values for both trend lines

print(period, period_uncer, g, g_uncer, chi_square1, chi_square2)

