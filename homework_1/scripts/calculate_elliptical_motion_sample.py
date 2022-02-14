# Purpose: Calculate Keplerian Orbits (Elliptical Motions) 
#          EAS417 (Satellite Meteorology): textbook Section 2.2
# Author: Johnny Luo
# Date: Feb 2022

# Import libraries and modules
#
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# Define constants
#
r_earth = 6.372e6 # radius of Earth
G = 6.673e-11     # gravitational const.
m_earth = 5.97e24 # mass of Earth

# Satellite parameters; 0/1: turn on/off the satellite
#
noaa = 0
molniya = 1

if(noaa):
    # NOAA polar orbiter
    semi_major = 7.229606e6
    epsilon = 0.00119958
    inclination = 98.97446 *np.pi/180
    omega_big = 29.31059*np.pi/180
    omega_small = 167.74754*np.pi/180

if(molniya):
    # Molniya orbit
    semi_major = 2.6554e7
    epsilon = 0.72
    i_angle = 63.4 *np.pi/180
    omega_big = 0*np.pi/180
    omega_small = 270*np.pi/180

# Calculate period, set up time list and three anomalies
#
T = 2*np.pi*np.sqrt((semi_major**3/G/m_earth)) # Period (Equation 2.4)
n = 2*np.pi/T # Mean motion constant (Equation 2.9)

time = np.linspace(0,T,7)  # Track 7 time stamps from 0 to one period
var1 = np.zeros(time.size) # one of the three anomalies (you need to figure out which one it is)
var2 = np.zeros(time.size) # one of the three anomalies (you need to figure out which one it is)
var3 = np.zeros(time.size) # one of the three anomalies (you need to figure out which one it is)

# define a function to solve Equation 2.8
#
def equation_2pt8(eccn_anomaly, epsn, mean_anomaly):
    return eccn_anomaly - epsn*np.sin(eccn_anomaly) - mean_anomaly

# Loop over time steps to calculate the three anomalies
#
for i in range(time.size):
    var1[i] = n*time[i] # one of the three anomalies (you need to figure out which one it is)
    
    # Find the root of the Equation 2.8
    var2[i] = optimize.fsolve(equation_2pt8,0,args=(epsilon,var1[i]))
    
    # one of the three anomalies (you need to figure out which one it is)
    if var2[i] <= np.pi:
        var3[i] = np.arccos((np.cos(var2[i])- epsilon)/(1 - epsilon*np.cos(var2[i])))
    else:
        var3[i] = 2*np.pi -np.arccos((np.cos(var2[i])- epsilon)/(1 - epsilon*np.cos(var2[i])))


# Define a background ellipse to be drawn (100 points)
#
theta_1  = np.linspace(0,2*np.pi,100)
ellipse_r = semi_major*(1-epsilon**2)/(1+epsilon*np.cos(theta_1))

# The points to be drawn on the ellipse (7 points)
ellipse_r_2 = semi_major*(1-epsilon**2)/(1+epsilon*np.cos(var3))
      
# Making figures
#
plt.figure(figsize=(12,12))
# Plot the ellipse in gray
plt.plot(ellipse_r*np.cos(theta_1)/1000,ellipse_r*np.sin(theta_1)/1000,linewidth=6,color='gray')
# Plot the circumscribed circle in read
plt.plot((semi_major*np.cos(theta_1)-semi_major*epsilon)/1000, semi_major*np.sin(theta_1)/1000,linewidth=2,color='red')
# Plot the three angles
plt.plot((semi_major*np.cos(var1)-semi_major*epsilon)/1000, semi_major*np.sin(var1)/1000,'.k',markersize=40,label="Variable 1")
plt.plot((semi_major*np.cos(var2)-semi_major*epsilon)/1000, semi_major*np.sin(var2)/1000,'.g',markersize=40,label="Variable 2")
plt.plot(ellipse_r_2*np.cos(var3)/1000,ellipse_r_2*np.sin(var3)/1000,'.r',markersize=40,label="Variable 3")
# Plot the Earth at the origin
plt.legend(fontsize=16)
plt.plot(0,0,'bo',markersize=40)
plt.text(-3000,-5000,'Earth',fontsize = 30)
plt.tick_params(axis="x",labelsize=16)
plt.tick_params(axis="y",labelsize=16)
plt.xlabel('Km',fontsize=18)
plt.ylabel('Km',fontsize=18)
plt.title('Calculating Elliptical Motion',fontsize = 20)
plt.grid()

plt.show()
