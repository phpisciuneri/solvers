# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 16:24:19 2014

@author: patrick
"""

import numpy as np
import matplotlib.pyplot as plt

# read points
profile = np.loadtxt( "profile.out" )
x  = profile[0,:]

plt.figure()

plt.plot( x, .5*np.tanh( 15*(x-1) ) + .5 , lw=3 )
plt.plot( x, -.5*np.tanh( 15*(x-2) ) + .5 , lw=3 )

# label figure
plt.title( 'sin(x)', fontsize=18, fontweight='bold')
plt.xlabel( 'x', fontsize=16, fontweight='bold' )
plt.ylabel( 'y(x)', fontsize=16, fontweight='bold' )
plt.legend( loc='upper right' )

# plot
plt.box('on')
plt.grid('on')
plt.show()
