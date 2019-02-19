# Uses the Euler integration method to solve the Klein-Gordon equation d^2\phi / dt^2 - \grad^2\phi +
# dV/d\phi = 0, where \phi is a scalar field dependent on time t and one spatial dimension x, and V is
# dependent upon \phi and is the potential.

# In the case of the massive potential we have that V= 1/2 * m^2 * \phi^2.

# We consider a spatial box between 0 < x < L.

# Import packages
import numpy as np
import math

# Provide calculation parameters
x_step = int(input('How many divisions should we have between x=0 and x=L?'))

# Optimal to have t_step = 0.2* x_step
t_step = int(input('How many timesteps should we calculate?'))

# Specify the mass of our potential
m = float(input('What is the mass of our potential?'))

# Because our solution will be sinusoidal/cosine, we need to specify which mode we want
n = int(input('How many sine curves should we plot?'))

# Calculate spacings between successive space and time points
delta = 1 / float(x_step)
h = 0.2 * delta

# define arrays in which we will store our results.

# phi is an array with x running down the rows and t running from left to right along the columns (we've done it this
# way around to satisfy gnuplot)
phi = np.zeros((x_step + 1, t_step + 1))
phi_dot = np.ones((x_step + 1, t_step + 1))
phi_dx_dx = np.zeros((x_step + 1, t_step + 1))

# now define the x vector in increments from 0 to 1 (L).
x = [0]
# x is a column vector that has all the x-values from 0 at the top to 1 (L) at the bottom, we loop through stacking the
# elements on top of each other.
for i in range(1, x_step + 1):
    x = np.vstack((x, i * delta))
# we can then populate our initial condition on phi dot in order to get the solution we want: NOTE THAT THIS INITIAL
    # CONDITION MAY CHANGE.
    phi_dot[i, 0] = math.sin(math.pi * n * x[i])

# create the array where we will store our data which will then be written out into a file.
DataOut = np.zeros((x_step + 1, t_step + 2))
DataOut[:, 0] = x[:, 0]


# defne our potential function
def v_potential(psi):
    return m * m * psi


# We begin the RK4 method
for j in range(0, t_step):
    for i in range(1, x_step):
        # first we work out the second spatial derivative at this point
        # NOTE we have to change this depending on where in the mesh the point lies: need some if statements
        if i == 1:
            phi_dx_dx[i, j] = (phi[3, j] - 3 * phi[1, j] + 2 * phi[0, j]) / (4 * delta * delta)
        elif i == x_step:
            phi_dx_dx[i, j] = (2 * phi[x_step - 2, j] - phi[x_step - 1, j] - phi[x_step, j]) / (2 * delta * delta)
        else:
            phi_dx_dx[i, j] = (phi[i + 1, j] - 2 * phi[i, j] + phi[i - 1, j]) / (delta * delta)

        # then we begin the calculating the Runge-Kutta constants for this point

        k_0 = float(phi_dot[i, j])
        l_0 = float(phi_dx_dx[i, j] + v_potential(phi[i, j]))
        k_1 = float(phi_dot[i, j] + h * l_0 / 2)
        l_1 = float(phi_dx_dx[i, j] + v_potential(phi[i, j]) + h * k_0 / 2)
        k_2 = float(phi_dot[i, j] + h * l_1 / 2)
        l_2 = float(phi_dx_dx[i, j] + v_potential(phi[i, j]) + h * k_1 / 2)
        k_3 = float(phi_dot[i, j] + h * l_2)
        l_3 = float(phi_dx_dx[i, j] + v_potential(phi[i, j]) + h * k_2)

        # we finally calculate the new phi and phi_dot
        phi[i, j + 1] = phi[i, j] + h * (k_0 + 2 * k_1 + 2 * k_2 + k_3) / 6
        phi_dot[i, j + 1] = phi_dot[i, j] + h * (l_0 + 2 * l_1 + 2 * l_2 + l_3) / 6

    # Define the data to be written into the file.
    DataOut[:, j + 1] = phi[:, j]

# print(phi_dot)
# Then information is written to the file.
np.savetxt('RK4_massive.dat', DataOut)
print('max = ', np.max(phi), ', min = ', np.min(phi))
