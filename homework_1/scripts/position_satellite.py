# Purpose: Calculate and plot satellite orbits 
#          EAS417 (Satellite Meteorology): textbook Section 2.5.1
# Author: Johnny Luo
# Date: Feb 2022

# Load python libraries and packages 
#
import numpy as np
from scipy import optimize

# Define constants
#
r_earth = 6.372e6
G = 6.673e-11
m_earth = 5.97e24

# Satellite parameters
#
noaa = 0
molniya = 1

if(noaa):
    # NOAA polar orbiter
    semi_major = 7.229606e6
    epsilon = 0.00119958
    i_angle = 98.97446 *np.pi/180
    omega_big_0 = 29.31059*np.pi/180
    omega_small_90 = 167.74754*np.pi/180

if(molniya):
    # Molniya
    semi_major = 2.6554e7
    epsilon = 0.72
    i_angle = 63.4 *np.pi/180
    omega_big_0 = 340*np.pi/180
    omega_small_0 = 270*np.pi/180


# Calculate period and three anomalies (angles)
#
T = 2*np.pi*np.sqrt((semi_major**3/G/m_earth)) # Period (Equation 2.4)
n = 2*np.pi/T # Mean motion constant (Equation 2.9)
    
# Orbital perturbations (Section 2.5.1)
#
J2 = 1.08263e-3; # coefficient of the quadrupole term (Appendix E) 
r_ee = 6.378137e+6; # equatorial radius of the Earth
#
# (Equations 2.12, 2.13, and 2.14)
#
dMdt = n*(1+1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-1.5)*(1-1.5*(np.sin(i_angle))**2));
domega_big_dt = -dMdt*(1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-2)*np.cos(i_angle));
domega_small_dt = dMdt*(1.5*J2*(r_ee/semi_major)**2*(1-epsilon**2)**(-2)*(2-2.5*(np.sin(i_angle))**2));

# Time as a single variable
time = (1/8) * T

# define function to calculate eccentric anamoly (Equation 2.8)
def eccentric_anomaly(e, epsn, M):
    return e - epsn*np.sin(e) - M

# Compute satellite position (following section 2.5.1)
#
    
M = n*time # Mean anomaly (increases evenly with time)
if M>=2*np.pi: M = np.mod(M,2*np.pi)

# Eccentric anamoly (find the root of the Equation 2.8)
e = optimize.fsolve(eccentric_anomaly,0,args=(epsilon,M))

# True anamaly
if e <= np.pi:
    theta = np.arccos((np.cos(e)- epsilon)/(1 - epsilon*np.cos(e)))
else:
    theta = 2*np.pi -np.arccos((np.cos(e)- epsilon)/(1 - epsilon*np.cos(e)))


# Update omega_small and omega_big (Equation 2.21)
#
omega_small = omega_small_0 + domega_small_dt*time
omega_big = omega_big_0 + domega_big_dt*time    
#    
# Cast into Cartesian coordinate (Equation 2.22)
#
radius = semi_major*(1-epsilon**2)/(1+epsilon*np.cos(theta));
x_0 = radius*np.cos(theta); 
y_0 = radius*np.sin(theta); 
z_0 = 0;
#
# First rotation (Equation 2.23)
#
x_1 = x_0*np.cos(omega_small)-y_0*np.sin(omega_small); 
y_1 = x_0*np.sin(omega_small)+y_0*np.cos(omega_small);
z_1 = z_0;
#
# Second rotation (Equation 2.24)
#
x_2 = x_1;
y_2 = y_1*np.cos(i_angle)-z_1*np.sin(i_angle);
z_2 = y_1*np.sin(i_angle)+z_1*np.cos(i_angle);    
#
# Third rotation (Equation 2.25)
#
x_3 = x_2*np.cos(omega_big)-y_2*np.sin(omega_big);
y_3 = x_2*np.sin(omega_big)+y_2*np.cos(omega_big);
z_3 = z_2;

# Convert the Cartesian coordinate to radius-declination-right ascension
# (Equation 2.26) 
#
r_s = np.sqrt((x_3**2+y_3**2+z_3**2));
delta_s = np.arcsin(z_3/r_s);
omega_s = np.arctan2(y_3,x_3);

lat_s = delta_s;
lon_s = omega_s-time*7.2921e-5; #7.2921e-5 is the rotation rate of Earth

print('The satellite position is:\n')
print('X = ',x_3,'\n')
print('Y = ',y_3,'\n')
print('Z = ',z_3,'\n')
print('lat = ',lat_s*180/np.pi,'\n')
print('lon = ',lon_s*180/np.pi,'\n')





