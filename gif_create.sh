#!/usr/bin/env bash
#this file is not meant to be run within the terminal itself, but is merely to keep a record of the commands
#required in order to create a .gif from data.

#change directory to be within the GR_finite_differencing project
cd PycharmProjects/GR_finite_differencing/KG_zero_potential
#gnuplot -persist
#
##insert the commands below:
#-e "set term gif animate delay 2"
##set output 'TestAnimation.gif' #change this title for a different number of waves
##set xrange [0:1]
##set yrange [-0.15:0.15] #change this yrange depending on how tall the peaks are in the data
#-e "do for[i=2:999]{plot 'AntikinkAnimation.dat' u 1:i w lines}"
#quit
#
##change the number after the : depending on how many timesteps we took
##it should always be 2:t_step-1


echo "Boo"
echo > gnuplot.in
for FILE in *; do
    echo "set term gif animate delay 2" >> gnuplot.in
    echo "set output 'RK4_debug.gif'" >> gnuplot.in
    echo "set xrange [0:1]" >> gnuplot.in
    echo "set yrange [-0.1:0.1]" >> gnuplot.in
    echo "do for[i=2:201]{plot 'RK4_debug.dat' u 1:i w lines}" >> gnuplot.in
done
gnuplot -persist gnuplot.in





