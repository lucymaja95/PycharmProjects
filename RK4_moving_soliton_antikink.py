import numpy as np
import math
x_step = int(input('How many divisions should we have between x=0 and x=L?'))
t_step = int(input('How many timesteps should we calculate?'))
v = int(input('What is v?'))
lam = int(input('What is lambda?'))
x_0 = float(input('What is x_0?'))
delta = 1 / float(x_step)
h = 0.2 * delta

# define arrays in which we will store our results.

# phi is an array with x running down the rows and t running from left to right along the columns (we've done it this
# way around to satisfy gnuplot)
phi = np.zeros((x_step + 1, t_step + 1))
# phi dot is a column vector of length the same as the x vector.
phi_dot = np.zeros((x_step + 1, 1))

# now define the x vector in increments from 0 to 1 (L).
x = [0]
# x is a column vector that has all the x-values from 0 at the top to 1 (L) at the bottom, we loop through stacking the
# elements on top of each other.
for i in range(1, x_step + 1):
    x = np.vstack((x, i * delta))

# create the array where we will store our data which will then be written out into a file.
DataOut = np.zeros((x_step + 1, t_step + 2))
DataOut[:, 0] = x[:, 0]

# introducing our initial conditions
for i in range(0, x_step + 1):
    phi[i, 0] = v * math.tanh(math.sqrt(lam / 2) * v * (-x[i] + x_0))

phi[0, :] = v * math.tanh(math.sqrt(lam / 2) * v * x_0)
phi[x_step, :] = v * math.tanh(math.sqrt(lam / 2) * v * (-1 + x_0))

for j in range(0, t_step):
    for i in range(1, x_step):
        # first we work out the second spatial derivative at this point
        phi_dx_dx = (phi[i+1, j] - 2 * phi[i, j] + phi[i - 1, j]) / (delta * delta)

        # then we begin the calculating the Runge-Kutta constants for this point
        c_1 = float(phi_dot[i, j])
        d_1 = float(phi_dx_dx - lam * phi[i, j] * (phi[i, j] * phi[i, j] - v * v))
        c_2 = float(phi_dot[i, j] + h * d_1 / 2)
        d_2 = float(phi_dx_dx - lam * (phi[i, j] + h * c_1 / 2) * ((phi[i, j] + h * c_1 / 2) * (phi[i, j] + h * c_1 / 2) - v * v))
        c_3 = float(phi_dot[i, j] + h * d_2 / 2)
        d_3 = float(phi_dx_dx - lam * (phi[i, j] + h * c_2 / 2) * ((phi[i, j] + h * c_2 / 2) * (phi[i, j] + h * c_2 / 2) - v * v))
        c_4 = float(phi_dot[i, j] + h * d_3)
        d_4 = float(phi_dx_dx - lam * (phi[i, j] + h * c_3) * ((phi[i, j] + h * c_3) * (phi[i, j] + h * c_3) - v * v))

        # we finally calculate the new phi and phi_dot
        phi[i, j + 1] = phi[i, j] + h * (c_1 + 2 * c_2 + 2 * c_3 + c_4) / 6
        phi_dot[i, j + 1] = phi_dot[i, j] + h * (d_1 + 2 * d_2 + 2 * d_3 + d_4) / 6

    # Define the data to be written into the file.
    DataOut[:, j + 1] = phi[:, j]

# Then information is written to the file.
np.savetxt('AntikinkAnimation.dat', DataOut)
print('max = ', np.max(phi), ', min = ', np.min(phi))
# print(DataOut)