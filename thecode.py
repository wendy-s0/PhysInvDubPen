##### Importing the necessary modules.
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

##### Defining the lengths and masses of the pendulums.
L1, L2 = 10, 10
m1, m2 = 10, 10
g = 9.81

##### Here, we are defining the coupled ODEs that we shall use. If you follow it carefully, all the maths presented here has been explained in the report. 
##### The variable names are quite self-explanatory - th1 = theta_1, w_1 = angular velocity 1 (or derivative of theta_1 here) etc.
def eqnprime(y, t, L1, L2, m1, m2):
    th1, w1, th2, w2 = y

    c, s = np.cos(th1-th2), np.sin(th1-th2)
    c2 = np.cos(2*(th1-th2))
    bot = ((2*m1)+m2-(m2*c2))
    
    th1p = w1
    th2p = w2
    w1p = (-g*(2*m1+m2)*np.sin(th1)-m2*g*s-2*s*m2*((th2p**2)*L2+(th1p**2)*L1*c))/(L1*bot)
    w2p = 2*s*(th1p**2*L1*(m1+m2)+g*(m1+m2)*np.cos(th1)+(th2p**2)*L2*m2*c)/(L2*bot)
    return th1p,w1p,th2p,w2p


##### The important line in which we change our initial conditions:
y0 = np.array([2.391, 0, 0.353, 0])

##### The lines below change how long (and at what time intervals) the 'solving'/integration of the ODEs is done at/for. For example, below, it is being done for
##### what would be 100 seconds in real life, at intervals of 0.01 seconds.
tmax, dt = 100, 0.01
t = np.arange(0, tmax+dt, dt)

##### The line which actually does the integration of the differential equations.
y = odeint(eqnprime, y0, t, args=(L1,L2,m1,m2))

##### Taking the outputs of the angles.
th1, th2 = y[:,0], y[:,2]

##### Defining the coordinates of the pendulums.
x1 = L1*np.sin(th1)
y1 = -L1*np.cos(th1)
x2 = x1+L2*np.sin(th2)
y2 = y1-L2*np.cos(th2)


##### The blocks below are to plot the graphs presented in the report.
##### First block: x-coordinates against time for both pendulums
plt.subplot(2,2,1)
plt.xticks([])
plt.yticks([])
timex = plt.plot(t,x1, color = 'black', linewidth=0.5), plt.plot(t,x2, color = 'orange', linewidth=0.5)

##### Second block: y-coordinates against time for both pendulums
plt.subplot(2,2,2)
plt.xticks([])
plt.yticks([])
timey = plt.plot(t,y1, color = 'black', linewidth=0.5), plt.plot(t,y2, color = 'orange', linewidth=0.5)

##### Third block: the path of the first pendulum
plt.subplot(2,2,3)
plt.xticks([])
plt.yticks([])
posa = plt.plot(x1,y1,color='red',linewidth=0.5)

##### Fourth block: the path of the second pendulum
plt.subplot(2,2,4)
plt.xticks([])
plt.yticks([])
posb = plt.plot(x2,y2,color='blue',linewidth=0.5)

##### This line saves the plot as a png (at 300dpi) to a folder. The file is saved as 'whateverfilename.png'.
plt.savefig('/Users/ratchatakornsotthivej/Desktop/untitled folder/whateverfilename.png', dpi=300)
