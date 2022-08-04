# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 14:33:29 2022

@author: agerlt
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

# get a bunch of evenly spaced numbers between 0 and 8 pi
realspace_numbers = np.linspace(0, 8*np.pi, (2**12)+1)

# apply the cosine function
cos = np.cos(realspace_numbers)

# do a fft of the cos function
f_cos = np.fft.fft(cos)

# plot them to see what happened
plt.figure()
plt.plot(realspace_numbers, cos*500, c='r', label='realspace')
plt.plot(realspace_numbers[1:100]*40, f_cos[1:100], c='k', label='transform')
plt.title('real and ft of cosine function')
plt.legend()
plt.tight_layout()

# now do it with a more complex function
rsn = np.linspace(0, 20*np.pi, (2**12)+1)
complex_function = 4.1*np.cos(3*rsn) + 1.2*np.sin(-4*rsn)
f_complex = np.fft.fft(complex_function)

# plot them to see what happened
plt.figure()
plt.plot(rsn, complex_function*500, c='r', label='realspace')
plt.plot(rsn[1:200]*20, f_complex[1:200]*100, c='k', label='transform')
plt.title("real and ft of complex function")
plt.ylim(-3500, 3500)
plt.legend()
plt.tight_layout()


# alternately, try going from fourier space to real space.
f_data = np.zeros(realspace_numbers.shape)
# add some random spikes
f_data[0] = 10
f_data[10] = 10
f_data[11] = 10
f_data[20] = 10
f_data[21] = 10
f_data[30] = 10
f_data[45:60] = 5
# transform back
data = np.fft.ifft(f_data)

# plot them to see what happened
plt.figure()
plt.plot(realspace_numbers[1:100]*40, f_data[1:100], c='k', label='transform')
plt.plot(realspace_numbers, data*500, c='r', label='realspace')
plt.title("ft function and realspace equivilant")
plt.legend()
plt.tight_layout()
