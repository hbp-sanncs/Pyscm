#!/usr/bin/env python
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
This programm sweeps over a weight given in the command line. The output is
saved in a .txt file as CSV in the analsis direction. There is also a gnuplot
script which plots the data
"""
import json
import sys

import numpy as np
import matplotlib.pyplot as plt
import pyscm
import pynam.data as data
import pynam.network as netw
import pynnless as pynl
import pyscm.parameters as par

n_min = 40
n_max = 50

# check simulator
if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <SIMULATOR>")
    sys.exit(1)

# Read in neuron data
with open("data/neuron_data.json", 'r') as outfile:
    temp_dict = json.load(outfile)

data_params = temp_dict["data_params"]
params = temp_dict["neuron_params"]
delay = temp_dict["delay"]
terminating_neurons = temp_dict["terminating_neurons"]
flag = temp_dict["Simple_Network"]
optimise_params = temp_dict["optimise_params"]
I, I_norm, fp, fn, I_start, I_norm_start, fp_start, fn_start = np.zeros(
    (8, (n_max - n_min)))
for i in xrange(n_max - n_min):
    n_samples = n_min + i
    data_params["n_samples"] = n_samples
    # Optimisation process
    opti = par.WeightOptimisation(data_params, params, delay,
                                  optimise_params,
                                  sys.argv[1], terminating_neurons, flag)
    opti.do_standard_stuff()
    I[i], I_norm[i], fp[i], fn[i], I_start[i], I_norm_start[i], fp_start[i], \
    fn_start[i] = opti.full_run()

# Write the data to the respective file
with open("analysis/" + "sample_sweep.txt", 'w') as outfile:
    outfile.write(
        "#weight,I, I_norm, fp ,fn ,I_start , I_norm_start, fp_start, fn_start\n")
    for i in xrange(n_min, n_max):
        outfile.write("%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n" % (
            i, I[i], I_norm[i], fp[i], fn[i], I_start[i], I_norm_start[i],
            fp_start[i], fn_start[i]))
