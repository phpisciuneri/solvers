# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:58:48 2014

@author: patrick
"""

import numpy as np
import matplotlib.pyplot as plt

# read points
profile = np.loadtxt( "profile.out" )
x  = profile[0,:]
yi = profile[1,:]
yf = profile[2,:]

plt.figure()

plt.plot( x, yf - yi, lw=3, label='Error' )

# label figure
plt.title( 'sin(x)', fontsize=18, fontweight='bold')
plt.xlabel( 'x', fontsize=16, fontweight='bold' )
plt.ylabel( 'y(x)', fontsize=16, fontweight='bold' )
plt.legend( loc='upper right' )

# plot
plt.box('on')
plt.grid('on')
plt.show()
