import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

maxima, grating_angle = np.loadtxt('grating.txt', unpack=True, skiprows=1)
maxima_water, angle_water = np.loadtxt('water.txt', unpack=True, skiprows=3)
# 'water.txt' has data for 7 trials
freq = np.loadtxt('frequency.txt', skiprows=1)


def uncertainty_cal(uncer):
    # This function calculate uncertainty when doing addition/subtraction
    a = 0
    for i in uncer:
        a += i**2
    return np.sqrt(a)


grating_space = 0.00001016  # this is calculated from '2500 lines per inch'
grating_angle = grating_angle*np.pi/180  # converting from degree to radiant
angle_water = angle_water*np.pi/180
angle_uncer = 0.005*np.pi/180
# uncertainty for angle measurement was estimated to be 0.005 degree, this
# includes instrument incapability and human error
freq = freq*1000000  # the measured frequency was in MHz
wave_length = abs(grating_space*np.sin(grating_angle)/maxima)
wave_length = wave_length.mean()  # wavelength of sodium lamp
WL_uncer = abs(grating_space*np.sin(grating_angle)/maxima)*(angle_uncer/grating_angle)
WL_uncer = uncertainty_cal(WL_uncer)  # wavelength uncertainty for sodium lamp
inverse_freq = 1/freq

lambda_water = []  # sound wavelength underwater
lambda_uncer = []  # uncertainty of sound wavelength underwater


for i in range(len(freq)):
    # 7 trails of data are stored in 'water.txt', so each loop goes through one trial

    WL = abs(maxima_water[i*4:i*4+4]*wave_length/np.sin(angle_water[i*4:i*4+4]))
    # WL is wavelength of sound calculated from each measurement
    uncer_percent = np.ones(4)*WL_uncer/wave_length+\
                    abs(angle_uncer/angle_water[i*4:i*4+4])
    uncer = abs(maxima_water[i*4:i*4+4]*wave_length/
                np.sin(angle_water[i*4:i*4+4]))*(uncer_percent)
    # uncer is the uncertainty of sound wavelength
    WL = WL.mean()
    lambda_water.append(WL)
    lambda_uncer.append(uncertainty_cal(uncer))
lambda_uncer = np.array(lambda_uncer)
lambda_water = np.array(lambda_water)

def model_func(inverse_freq, speed):
    return speed*inverse_freq


popt, pcov = curve_fit(model_func, inverse_freq, lambda_water, sigma=lambda_uncer)
Best_fit = plt.plot(inverse_freq, model_func(inverse_freq, popt), label='Line of Best Fit')
Errorbar = plt.errorbar(inverse_freq, model_func(inverse_freq, popt),
                 yerr=lambda_uncer, linestyle='none', label='Measured Data')
chi_square = (1/(len(inverse_freq)-2))*sum(((lambda_water -
                                             model_func(inverse_freq, popt))
                                            / lambda_uncer)**2)
plt.xlabel('Inverse of Frequency (s)')
plt.ylabel('Wavelength (m)')
plt.title('Change in Sound Wavelength Underwater With Respect to \n'
          'Change in Inverse of Sound Frequency')
plt.legend()
plt.show()
print(wave_length, WL_uncer, popt[0], pcov[0][0], chi_square)
# this prints sodium lamp wavelength, sodium lamp wavelength uncertainty
# speed of sound underwater, the speed uncertainty, and chi-square
