# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:20:08 2022

@author: agerlt
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp2d

plt.close('all')
# load up data from spinodal structure
original_data = plt.imread("spinodal.png")[:, :, 0]

# display it (not important, just helpful to see)
plt.imshow(original_data)

# make x and y axis for regular data
x, y = np.arange(original_data.shape[0]), np.arange(original_data.shape[0])
xx, yy = np.meshgrid(x, y)
# make an interpolation function.
f = interp2d(x, y, original_data)


# create some non_uniform data to sample from your interpolation data
non_uniform_x = np.sort(np.random.random(20*(x.size))*(x.size))
non_uniform_y = np.sort(np.random.random(20*(y.size))*(y.size))
# get some data
nu_all_data = f(non_uniform_x, non_uniform_y)

# this is a LOT of data: cut it down to the same size as the original dataset
indices = np.random.choice(np.arange(nu_all_data.size), original_data.size)
nu_xx, nu_yy = np.meshgrid(non_uniform_x, non_uniform_y)
nu_x = (nu_xx.flatten())[indices]
nu_y = (nu_yy.flatten())[indices]
nu_data = (nu_all_data.flatten())[indices]

# plot it to make sure the function worked
plt.figure()
plt.contourf(nu_xx[::50, ::50], nu_yy[::50, ::50], nu_all_data[::50, ::50])
plt.xlim(0, 128)
plt.ylim(128, 0)
# and save out a uniform and non_uniform sampling of the data
np.savetxt("uniform_data.txt", np.vstack([
    xx.flatten(),
    yy.flatten(),
    original_data.flatten()
    ]).T)
np.savetxt("non_uniform_data.txt", np.vstack([
    nu_x,
    nu_y,
    nu_data
    ]).T)
