# -*- coding: utf-8 -*-
"""
Created on Thu Feb 04 13:03:58 2016

@author: JWMize
"""

## Initialize
import random as random;
import time as time;
import matplotlib.pyplot as plt;
import numpy as np;
import math as m;

pi = np.pi; # Import the value of pi

x = []; # Array of horizontal data (m)
dx = 1; # Distance between x steps (m)
z = []; # Array of vertical data (m)
dirt = []; # Array of dirt elevations (m)
zideal = []; # Array of ideal vertical data (m)
net = []; # Array of net gain/loss of dirt at each point (m)
w = []; # Array of weathering rates at each spot
slope = []; # Array of slopes between x steps (unitless)
q = []; # Array of fluxes between x steps (m/yr)

noise = 0; # For adding topographic noise to the systems (m)

rho_soil = 2650; # Density of soil (kg/m^3)
rho_rock = 2750 # Density of granite (kg/m^3)
w_0 = 1; # Base weathering rate (m/yr)
k_eff = 0.01; # landscape movement efficency
kappa = k_eff / rho_soil; # topographic diffusivity
tmax = 1000; # Length of experiment (yr)
dt = 1; # Number of years between time steps (yr)
hstar = 0.3; # Characteristic length scale for dirt (m)

extent = 1000; # Length of arrays

for i in range(0,extent,dx): # Loop to initialize arrays
    x.append(i); # Populates x
    z.append(0); # Lengthens z array (because I don't know how else to do it)
    dirt.append(0); # See above
    w.append(0);
    net.append(0);
    q.append(0);
    slope.append(0);
    zideal.append(0);
    noise = noise + (random.random()-0.5)*10; # Generates noise, dependant on earlier noise
    z[i] = 1000 + noise; # Sets z plateau at 1000 elevation, adds noise
    dirt[i] = z[i]; # Initializes dirt at elevation of z
    
#    zideal[i] = w_0/(2*kappa)*(extent**2-x[i]**2) # Populates zideal

q.append(0); # Array q needs to be long

plt.plot(x,z, label = 'Model');

plt.axis([x[0], x[extent-1], 0, 1500]);
plt.title('Cross section of hillslope elevation');
plt.legend();
plt.show();
time.sleep(1); # Pauses 1 second

## Run

for j in range(0,tmax,dt): # Do actions once per time step
    plt.clf(); # Clear old stuff    
    z.append(0); # Need long z to calculate s
    
    for k in range(0,extent,dx): # Calculate full arrays once per time step
        
        w[k] = w_0*m.exp(-(dirt[i]-z[i])/hstar); # Sets value of weathering rate at each x
        slope[k] = (z[k+1]-z[k])/dx; # Sets value of slope between x points
        q[k+1] = slope[k] * -k_eff; # Sets value of flux between x points
        
        if (q[k] - q[k+1]) * dt > dirt[k] - z[k] + w[k] * dt: # Prevent non-existant dirt from being
            q[k+1] = dirt[k] - z[k] + w[k] / dt; # Removes all dirt from an area over a single time step.
        else:
            ();
        
    q[extent/dx] = 0; # Sets flux in/out from right side to 0
    q[0] = 0; # Sets flux in/out from left side to 0
    z.pop(extent); # Shorten z so it can be plotted against x

    for l in range(0,extent,dx): # These arrays can't be calculated until above arrays are calculated
        
        net[l] = dt * (q[l] - q[l+1]) + dt * w[l] * (rho_rock - rho_soil) / rho_rock; # Calculates gain/loss of soil at each x
        dirt[l] = dirt[l] + net[l];
        z[l] = z[l] - w[l]*dt
        
    dirt[0] = z[0]; # Stream sweeps away all the dirt
    dirt[extent/dx-1] = z[extent/dx-1]; # Stream sweeps away all the dirt
    
    plt.axis([x[0], x[extent-1], 0, 1500]);
    plt.title('Cross section of hillslope elevation');
    plt.plot(x,z); # Plot the new calculated z
    plt.plot(x,dirt); # Plot the dirt elevation
#    plt.plot(x,zideal); # Plot the ideal z curve
    plt.show(); # Show the plot
#    time.sleep(0.1) # Pauses 0.1 second