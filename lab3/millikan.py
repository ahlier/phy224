import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

po = 875.3
pa = 0
g = 9.8
dv = 1.827*10**(-5)
d = 0.006
constant = 6 * np.pi * dv * d * np.sqrt(9 * dv / (2 * g * (po - pa)))
charge = [] # this stores the experimental charge for each trial
uncertainty = []  # uncertainty for charge
oil_radius = []
for i in range(1, 57):
    frame = pd.read_excel('data/{}.xlsx'.format(i))
    trial, voltage_up, stop_voltage = np.loadtxt('voltage.txt', unpack=True,
                                                 skiprows=1)
    frame_lst = []
    for index, row in frame.iterrows():
        # because some entries are '#NV', which needs to be filtered out
        if row[0] != '#NV':
            frame_lst.append(row[0])

    b = np.array(frame_lst[:21])    # for calculation of velocity going
    c = np.array(frame_lst[-20:])
    v_up = abs(((b[1:] - b[:-1])  / (520*1000)) / 0.1).mean()
    v_ter = abs(((c[1:] - c[:-1]) / (520*1000)) / 0.1).mean()
    Q = float(constant * (v_ter + v_up) * v_ter ** 0.5 / voltage_up[i-1])
    charge.append(Q)
    uncertainty.append(Q*(5/520+0.5/voltage_up[i-1]))
    # the uncertainty for frame is 520 +/- 1, therefore the uncertainty of
    # velocity obtained from the frame is is velocity multiply by (2/520), and
    # the smallest digit for voltage reading is 1V, therefore, the uncertainty
    # is 0.5 V.
    oil_radius.append(np.sqrt(9*dv*v_ter/(2*g*(po-pa)))*1000)

    v_up_std = np.std(((b[1:] - b[:-1])  / (520*1000)) / 0.1)
    v_ter_std = np.std(((c[1:] - c[:-1]) / (520 * 1000)) / 0.1)


count = []
total_count = []
smallest_division = 10**(-20)

for i in range(10, 100):
    # The charges are dividing from 1*10^(-19)C to 1*10^(-18)C, with smallest
    # division of 1*10^(-20)C

    total_count.append(0)
    for q in range(len(charge)):
        if charge[q]%(i*smallest_division) < uncertainty[q] or \
                charge[q]//(i*smallest_division) != \
                (charge[q]+uncertainty[q])//(i*smallest_division):
            # this if statement check if the division fall within the range of
            # experimental charge with uncertainty
            count.append(i)
            total_count[-1] += 1
_ = plt.plot(np.arange(10, 100), total_count)
distributin = plt.hist(count, bins=30, rwidth=0.9)
plt.xlabel('Charge (C)')
plt.ylabel('Frequency')
plt.title('Count of Divisible Charges in Range 1*10^(-19) and 1*10^(-18) \n'
          'With Air Buoyancy Being Zero')
plt.show()


# The following will produce a histogram of frequency of different radius of oil drop
average_radius = np.array(oil_radius).mean()
radius_distribution = plt.hist(oil_radius, bins = int(max(oil_radius)//0.0003)+1, rwidth=0.9)
# the histogram will have 0.0003mm as the width of each each bin
plt.xlabel('Radius of oil drop (mm)')
plt.ylabel('Frequency')
plt.title('Number of Occurrence of Different Radius of Oil Drop')
plt.show()
print(average_radius)



