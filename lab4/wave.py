import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

wavelength_2a = np.loadtxt('2a_wavelength.TXT', unpack=True, skiprows=1)
frequency_2a = np.loadtxt('2a_frequency.txt', skiprows=1, unpack=True)
depth_2b = np.loadtxt('2b_depth.txt',  skiprows=1, unpack=True)
wavelength_2b = np.loadtxt('2b_wavelength.txt', skiprows=1, unpack=True)
wavelength_2a, wavelength_2b = wavelength_2a/100, wavelength_2b/100
depth_2b = depth_2b/1000

def model_func(x, a):
    return a*x

# Exercise 2a
mean_lambda_2a = []     # wavelength
lambda_uncer_2a = []    # wavelength uncertainty
graph_domain = np.linspace(0.1, 0.23, 100)
for i in range(5):
    # 5 trails were taken
    mean_lambda_2a.append(wavelength_2a[i:5*(1+i)].mean())
    lambda_uncer_2a.append(np.std(wavelength_2a[i:5*(1+i)]))

plt.subplot(2, 1, 1)
# plotting the period vs frequency graph with its residual graph as a sub-graph
graph_2a = plt.errorbar([1/i for i in frequency_2a], mean_lambda_2a,
                        yerr=lambda_uncer_2a, fmt='o')
popt_2a, pcov_2a = curve_fit(model_func, [1/i for i in frequency_2a], mean_lambda_2a,
                             absolute_sigma=True, sigma=lambda_uncer_2a)
best_fit = plt.plot(graph_domain, model_func(graph_domain, *popt_2a))

plt.ylabel('Wavelength (m)')
plt.title('Change of Period of Ripple Generator With Respect \n'
          'to Change of Wavelength')

plt.subplot(2, 1, 2)
graph_residual = plt.errorbar([1/i for i in frequency_2a],
                              mean_lambda_2a-model_func(np.array([1/i for i in frequency_2a]), *popt_2a),
                              fmt='o', yerr=lambda_uncer_2a)
horizontal_line = plt.plot(graph_domain, np.zeros(len(graph_domain)))
plt.xlabel('Period (s)')
plt.ylabel('Residual')
plt.show()

chi_square = (1/(len(frequency_2a)-2))*sum((((mean_lambda_2a-
                                              model_func(np.array([1/i for i in frequency_2a]),
                                                         *popt_2a))/lambda_uncer_2a)**2))
print(chi_square, *popt_2a, *np.diag(pcov_2a))



# Exercise 2b
mean_lambda_2b = []     # wavelength for 2b
lambda_uncer_2b = []    # wavelength uncertainty
speed_2b = np.sqrt(9.8*depth_2b)
graph_domain_2b = np.linspace(0.2, 0.3, 100)
for i in range(4):
    mean_lambda_2b.append(wavelength_2b[i:4*(1+i)].mean())
    lambda_uncer_2b.append(np.std(wavelength_2b[i:4*(1+i)]))

plt.subplot(2, 1, 1)
graph_2b = plt.errorbar(speed_2b, mean_lambda_2b,
                        yerr=lambda_uncer_2b, fmt='o', label='Measured Data')
popt_2b, pcov_2b = curve_fit(model_func, speed_2b, mean_lambda_2b,
                             absolute_sigma=True, sigma=lambda_uncer_2b)
best_fit_2b = plt.plot(graph_domain_2b, model_func(graph_domain_2b, *popt_2b))

plt.ylabel('Wavelength (m)')
plt.title('Change in Wavelength With Respect to Change in \n'
          'Wave Speed by Varying Water Depth')

plt.subplot(2, 1, 2)
graph_residual = plt.errorbar(speed_2b, mean_lambda_2b-model_func(speed_2b, *popt_2b),
                              fmt='o', yerr=lambda_uncer_2b)
horizontal_line = plt.plot(graph_domain_2b, np.zeros(len(graph_domain_2b)))
plt.xlabel('Wave Speed (m/s)')
plt.ylabel('Residual')
plt.show()
chi_square = (1/(len(speed_2b)-2))*sum((((mean_lambda_2b-
                                          model_func(speed_2b, *popt_2b))/
                                         lambda_uncer_2b)**2))
print(popt_2b, pcov_2b, chi_square)



# Exercise 4:
plt.subplot(2, 1, 1)
slit_width, angle = np.loadtxt('Exercise_4.txt', skiprows=1, unpack=True)
slit_width = slit_width/100     # slit width was measured in cm
spread_width = np.sin(angle*np.pi/180)  # convert degree to radiant
spread_uncer = np.cos(angle*np.pi/180)*(5*np.pi/180)
frequency_4 = 20
depth_4 = 0.0092

x_range_4 = np.linspace(20, 70, 100)
popt_4, pcov_4 = curve_fit(model_func, [1/i for i in slit_width], spread_width,
                           absolute_sigma=True, sigma=spread_uncer)
best_fit_4 = plt.plot(x_range_4, model_func(x_range_4, *popt_4))
errorbar_4 = plt.errorbar([1/i for i in slit_width], spread_width, fmt='o',
                          yerr=spread_uncer)
plt.ylabel('Width of Spread (sin\u03F4)')
plt.title('Change in Inverse of Slit Width With Respect to \n'
          'Change in Width of Wave Spread')

plt.subplot(2, 1, 2)
residual_4 = plt.errorbar([1/i for i in slit_width], spread_width-
                          model_func(np.array([1/i for i in slit_width]), *popt_4),
                          fmt='o', yerr=spread_uncer)
horizontal_4 = plt.plot(x_range_4, np.zeros(len(x_range_4)))
plt.xlabel('Inverse of Slit Width (m^-1)')
plt.ylabel('Residual')
plt.show()
chi_square = (1/(len(spread_width)-2))*sum((((spread_width-
                                          model_func(np.array([1/i for i in slit_width]), *popt_4))/
                                         spread_uncer)**2))
print(chi_square, *popt_4, pcov_4)



# Exercise 5:
plt.subplot(2, 1, 1)
separation, m1, m2 = np.loadtxt('Exercise_5.txt', unpack=True, skiprows=1)
separation = separation/100
frequency_5 = 20
depth_5 = 0.0085
angle_5 =np.hstack([m1*np.pi/180, m2*np.pi/180])
angle_uncer = np.hstack([np.cos(m1)*(3*np.pi/180), np.cos(m2)*(3*np.pi/180)])
m_over_d = np.hstack([np.full(5, 0.5)/separation,np.full(5, 1.5)/separation])
x_range_5 = np.linspace(0, 35, 100)

popt_5, pcov_5 = curve_fit(model_func, m_over_d, np.sin(angle_5),
                           absolute_sigma=True, sigma=angle_uncer)
best_fit_5 = plt.plot(x_range_5, model_func(x_range_5, *popt_5))
errorbar_5 = plt.errorbar(m_over_d, np.sin(angle_5), fmt='o',
                          yerr=angle_uncer)
plt.ylabel('Angular Spread (sin\u03F4)')
plt.title('Order of Interference Minima Divides by Wave Sources \n'
          'Separation VS Angular Spread')

plt.subplot(2, 1, 2)
residual_5 = plt.errorbar(m_over_d, np.sin(angle_5)-
                          model_func(np.array(m_over_d), *popt_5),
                          fmt='o', yerr=angle_uncer)
horizontal_5 = plt.plot(x_range_5, np.zeros(len(x_range_5)))
plt.xlabel('Order of Minima / Wave Sources Separation (m^-1)')
plt.ylabel('Residual')
plt.show()

chi_square = (1/(len(m1)-2))*sum((((np.sin(angle_5)-
                                          model_func(np.array(m_over_d), *popt_5))/
                                   angle_uncer)**2))
print(chi_square, popt_5, pcov_5)


