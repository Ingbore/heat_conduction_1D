# -*- coding: utf-8 -*-
## Introduction
'''
This code was written to model heat diffusion of a block of magma in two dimensions
By JWM and REH for 3/16/2016
'''

import landlab as ll
import numpy as np
import pylab

## Initialize

# Constants
Ts = 0. # Surface temperature, C
Tmagma = 1500. # Temperature of intruding magma, C
Qm = .04 # Mantle heat flow, W/m^2
k = 2. # Thermal conductivity, W/m*K
rho = 2700. # Density, kg/m^3
c = 2000. # Speific heat capacity, J/kg*K
xmax = 301 # Number of horizontal grid spaces
dx = 100. # Distance between grid spaces (both x and z), m
zmax = 301 # Number of vertical grid spaces
tmax = np.pi*10**7*200000. # End of simulation, s (200000 years)

# Variables
kappa = k/(rho*c) # Thermal diffusivity, m^2/s
dt = dx**2/(kappa*2)/4 # Maximum time step to ensure stability, s

# Arrays and Grids
mg = ll.RasterModelGrid(zmax,xmax,dx) # Temperature grid

z = mg.add_zeros('node','Depth')
x = mg.add_zeros('node','Extent')
T = mg.add_zeros('node','Temperature')
Q = mg.add_zeros('link','Flux')
#Tgrad = np.zeros(mg.number_of_links)
dQdA = np.zeros(mg.number_of_links)

# Initialize Arrays
z[:] = mg.node_y # Vertical array of all nodes
x[:] = mg.node_x # Horizontal array of all nodes
T[:] = z * Qm/k + Ts # Temperature array, C

for i in range (mg.number_of_nodes): # Loop to create a diapir in the bottom of the model
    if 10000 < int(x[i]):
        if int(x[i]) < 20000:
            if int(z[i]) > 20000:
                T[i] = Tmagma
            else:
                pass
        else:
            pass
    else:
        pass
    pass
# Boundary conditions
mg.set_closed_boundaries_at_grid_edges(False, True, True, True)

#Allgrads = mg.calculate_gradients_at_links(T)
#Tgrad[mg.active_links] = Allgrads[mg.active_links]



## Run
for j in range (int(tmax/dt)):
    Tgrad = mg.calculate_gradients_at_active_links(T) # Temperature gradient with respect to location
    Q = -k * Tgrad
    dQdd = mg.calculate_flux_divergence_at_nodes(Q)
    T += -1/(rho*c)*dQdd*dt
    """
    -What I want to do is first get a change of temperature with respect to 
    location (x and z). 
    -With the change in temperature based on location, a flux can be determined.
    -Use the flux to determine a change in flux with respect to location (x
    and z)
    -Finally, use the change in flux with respect to location, the density, the
    specific heat capacity, and the change in time to calculate Temperatures
    each time step.
    
    ---Plotting it every 10 or 50 time steps
    """
    
## Finalize
# Get a 2D array version of the elevations
Tr = mg.node_vector_to_raster(T)

# Create a shaded image
pylab.close()  # clear any pre-existing plot
im = pylab.imshow(Tr, cmap=pylab.cm.RdBu, extent=[0,xmax*dx,0,zmax*dx],
                                          origin='lower')
# add contour lines with labels
cset = pylab.contour(Tr, extent=[0,xmax*dx,zmax*dx,0], hold='on',origin='image')
pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
# add a color bar on the side
cb = pylab.colorbar(im)
cb.set_label('Temperature (C)')

# add a title and axis labels
pylab.title('Simulated cross section of heat diffusion')
pylab.xlabel('Distance (m)')
pylab.ylabel('Distance (m)')
# Display the plot
pylab.show()
print('Run time = '+str(dt*j)+' seconds')