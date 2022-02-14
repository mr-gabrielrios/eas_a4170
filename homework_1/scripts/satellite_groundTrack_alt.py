#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 10:30:26 2022

@author: gabriel
"""

# Purpose: Calculate and plot satellite orbits 
#          EAS417 (Satellite Meteorology): textbook Section 2.5
# Author: Johnny Luo
# Date: Feb 2022

# Setting up the python environments
#
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable 
from scipy import optimize
import cartopy.crs as ccrs

# Define constants
#
r_earth = 6.372e6
G = 6.673e-11
m_earth = 5.97e24

# Satellite parameters: 1 to turn on and 0 to turn off
#
noaa = 1
molniya = 0
satellite = 'NOAA' if noaa else 'Molniya'

if(noaa):
    # NOAA polar orbiter
    semi_major = 7.229606e6
    epsilon = 0.00119958
    i_angle = 98.97446 *np.pi/180
    omega_big_0 = 29.31059*np.pi/180
    omega_small_0 = 167.74754*np.pi/180

if(molniya):
    # Molniya
    semi_major = 2.6554e7
    epsilon = 0.72
    i_angle = 63.4 *np.pi/180
    omega_big_0 = 340*np.pi/180
    omega_small_0 = 270*np.pi/180
    
# Calculate period, set up time list and three anomalies
#
T = 2*np.pi*np.sqrt((semi_major**3/G/m_earth)) # Period (Equation 2.4)

print(36400/T)

n = 2*np.pi/T # Mean motion constant (Equation 2.9)

# Orbital perturbations (Section 2.5.1)
# (Equations 2.12, 2.13, and 2.14)
#
J2 = 1.08263e-3; # coefficient of the quadrupole term (Appendix E) 
r_ee = 6.378137e+6; # equatorial radius of the Earth

dMdt = n*(1+1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-1.5)*(1-1.5*(np.sin(i_angle))**2));
domega_big_dt = -dMdt*(1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-2)*np.cos(i_angle));
domega_small_dt = dMdt*(1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-2)*(2-2.5*(np.sin(i_angle))**2));

# Set up time as a list; all parameters will be evaluated at these time stamps
#
number_of_orbits = 2
time = np.linspace(0*T,number_of_orbits*T,1000) # Time increases evenly (chopped into 1000 pieces)

#
# All parameters are functions of time 
#
M = np.zeros(time.size) # Mean anomaly
e = np.zeros(time.size) # Eccentric anomaly
theta = np.zeros(time.size) # True anomaly
radius = np.zeros(time.size) # radius (from Earth to satellite)
omega_small = np.zeros(time.size)
omega_big = np.zeros(time.size)

# Variables in rotation of axis
#
x_0 = np.zeros(time.size); y_0 = np.zeros(time.size); z_0 = np.zeros(time.size)
x_1 = np.zeros(time.size); y_1 = np.zeros(time.size); z_1 = np.zeros(time.size)
x_2 = np.zeros(time.size); y_2 = np.zeros(time.size); z_2 = np.zeros(time.size)
x_3 = np.zeros(time.size); y_3 = np.zeros(time.size); z_3 = np.zeros(time.size)

r_s = np.zeros(time.size); delta_s = np.zeros(time.size); omega_s = np.zeros(time.size)
lat_s = np.zeros(time.size); lon_s = np.zeros(time.size)

# define function to calculate eccentric anamoly (Equation 2.8)
def eccentric_anomaly(e, epsn, M):
    return e - epsn*np.sin(e) - M

# Loop over time steps to calculate the three anomalies
#
for i in range(time.size):
    
    M[i] = n*time[i] # Mean anomaly (increases evenly with time)
    if M[i]>=2*np.pi:
        M[i] = np.mod(M[i],2*np.pi)
    
    # Eccentric anamoly (find the root of the Equation 2.8)
    e[i] = optimize.fsolve(eccentric_anomaly,0,args=(epsilon,M[i]))
    if e[i] >= 2*np.pi: 
        e[i] = np.mod(e[i],2*np.pi) # Make sure e[i] is within [0,2*np.pi]
    
    # True anamaly
    if e[i] <= np.pi:
        theta[i] = np.arccos((np.cos(e[i])- epsilon)/(1 - epsilon*np.cos(e[i])))
    else:
        theta[i] = 2*np.pi -np.arccos((np.cos(e[i])- epsilon)/(1 - epsilon*np.cos(e[i])))
        
    # Update omega_small and omega_big (Equation 2.21)
    #
    omega_small[i] = omega_small_0 + domega_small_dt*time[i]
    omega_big[i] = omega_big_0 + domega_big_dt*time[i]    
    #    
    # Cast into Cartesian coordinate (Equation 2.22): turn (r,theta) to (x, y, z)
    #
    radius[i] = semi_major*(1-epsilon**2)/(1+epsilon*np.cos(theta[i]));
    x_0[i] = radius[i]*np.cos(theta[i]);
    y_0[i] = radius[i]*np.sin(theta[i]);
    z_0[i] = 0;
    #
    # First rotation (Equation 2.23)
    #
    x_1[i] = x_0[i]*np.cos(omega_small[i])-y_0[i]*np.sin(omega_small[i]);
    y_1[i] = x_0[i]*np.sin(omega_small[i])+y_0[i]*np.cos(omega_small[i]);
    z_1[i] = z_0[i];
    #
    # Second rotation (Equation 2.24)
    #
    x_2[i] = x_1[i];
    y_2[i] = y_1[i]*np.cos(i_angle)-z_1[i]*np.sin(i_angle);
    z_2[i] = y_1[i]*np.sin(i_angle)+z_1[i]*np.cos(i_angle);    
    #
    # Third rotation (Equation 2.25)
    #
    x_3[i] = x_2[i]*np.cos(omega_big[i])-y_2[i]*np.sin(omega_big[i]);
    y_3[i] = x_2[i]*np.sin(omega_big[i])+y_2[i]*np.cos(omega_big[i]);
    z_3[i] = z_2[i];


    # Convert the Cartesian coordinate to radius-declination-right ascension
    # (Equation 2.26) 
    #
    r_s[i] = np.sqrt((x_3[i]**2+y_3[i]**2+z_3[i]**2));
    delta_s[i] = np.arcsin(z_3[i]/r_s[i]);
    omega_s[i] = np.arctan2(y_3[i],x_3[i]);
     
    lat_s[i] = delta_s[i];
    lon_s[i] = omega_s[i]-time[i]*7.2921e-5; #7.2921e-5 is the rotation rate of Earth
            
fig = plt.figure()
ax = plt.axes(projection=ccrs.PlateCarree())    
ax.coastlines(resolution='50m', color='black', linewidth=1)
ax.set_extent([-180, 180, -90, 90])
xtick = np.arange(-180, 180 + 60, 60)
ytick = np.arange(-90, 90 + 30,30)
xlabel = ['180W','120W','60W','0','60E','120E','180E']
ylabel = ['90S','60S','30S','0','30N','60N','90N']
ax.set_xticks(xtick)
ax.set_yticks(ytick)
ax.set_xticklabels(xlabel)
ax.set_yticklabels(ylabel)
ax.set_title('{0} Orbit - {1:.1f} complete orbits'.format(satellite, number_of_orbits))
points = ax.scatter(lon_s*180/np.pi,
                    lat_s*180/np.pi,
                    c = time/T,
                    cmap='plasma',
                    transform = ccrs.PlateCarree())

cax = make_axes_locatable(ax)
cax_ = cax.new_horizontal(size="3%", pad=0.15, axes_class=plt.Axes)
fig.add_axes(cax_)
colorbar = fig.colorbar(points, cax=cax_)
colorbar_label = colorbar.set_label('Orbit number', rotation=270, labelpad=20)

fig.tight_layout()