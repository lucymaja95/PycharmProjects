import numpy as np
import math
x_step = int(input('How many divisions should we have between x=0 and x=L?'))
t_step = int(input('How many timesteps should we calculate?'))
n = int(input('How many sine curves should we plot?'))
m = int(input('What is the mass?'))

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
# elements on top of each other.
# x is a column vector that has all the x-values from 0 at the top to 1 (L) at the bottom, we loop through stacking the
for i in range(1, x_step + 1):
    x = np.vstack((x, i * delta))
# we can then populate our initial condition on phi dot in order to get the solution we want.
    phi_dot[i, 0] = math.sin(math.pi * n * x[i])

# create the array where we will store our data which will then be written out into a file.
DataOut = np.zeros((x_step + 1, t_step + 2))
DataOut[:, 0] = x[:, 0]

# first line of code:
for i in range(1, x_step):
    phi[i, 1] = phi[i, 0] + h * phi_dot[i, 0]

for j in range(1, t_step):
    for i in range(1, x_step):
        phi[i, j + 1] = 2 * phi[i, j] - phi[i, j - 1] - (m * m) * (h * h) * phi[i, j - 1] + (h * h) * (phi[i + 1, j] - 2 * phi[i, j] + phi[i - 1, j]) / (delta * delta)

    # Define the data to be written into the file.
    DataOut[:, j + 1] = phi[:, j]

# Then information is written to the file.
np.savetxt('2_wave_massive_gif_data.dat', DataOut)

print(DataOut)
