#!/usr/bin/env python
#  -*- coding: utf-8 -*-

#   PySCM -- Python Spike Counter Model
#   Copyright (C) 2016 Christoph Jenzen, Andreas Stöckel
#
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

import json
import sys

import matplotlib.pyplot as plt
import pyscm
import pynam.data as data
import pynam.network as netw
import pynnless as pynl

# check simulator
if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <SIMULATOR>")
    sys.exit(1)

# Read in neuron data
with open("data/neuron_data.json", 'r') as outfile:
    dict = json.load(outfile)

data_params = dict["data_params"]
params = dict["neuron_params"]
delay = dict["delay"]
terminating_neurons = dict["terminating_neurons"]
flag = dict["Simple_Network"]

# Read in neuron weights, some still have to be set manually
with open("data/optimised_weights.json", 'r') as outfile:
    weights = json.load(outfile)

# Generate BiNAM
mat_in = data.generate(data_params["n_bits_in"], data_params["n_ones_in"],
                       data_params["n_samples"])
mat_out = data.generate(data_params["n_bits_out"], data_params["n_ones_out"],
                        data_params["n_samples"])
print "Data generated!"

# set up simulator
scm = pyscm.SpikeCounterModel(mat_in, mat_out)

sim = pynl.PyNNLess(sys.argv[1])
net, input_indices, _, input_times = scm.build(params=params, weights=weights,
                                               delay=delay,
                                               terminating_neurons=terminating_neurons,
                                               flag=flag)
print "Preparations done"

# Simulation
res = sim.run(net)
print "Simulation done"

# Plot
for pIdx, pop in enumerate(res):
    if (not "spikes" in pop):
        continue
    output_times = pop["spikes"]
    fig = plt.figure()
    ax = fig.gca()
    for i, spikes in enumerate(output_times):
        ax.plot(spikes, [i + 1] * len(spikes), '.', color=[0, 0, 0])
    ax.set_xlabel("Spike time [ms]")
    ax.set_ylabel("Neuron index")
    ax.set_title("Population " + str(pIdx))
# plt.show()

# Analyse
output_times, output_indices = netw.NetworkInstance.match_static(input_times,
                                                                 input_indices,
                                                                 res[0][
                                                                     "spikes"])

analysis = netw.NetworkAnalysis(input_times=input_times,
                                input_indices=input_indices,
                                output_times=output_times,
                                output_indices=output_indices,
                                data_params=data_params,
                                mat_in=mat_in, mat_out=mat_out)

I, I_norm, fp, fn, I_start, I_norm_start, fp_start, fn_start = pyscm.scm_analysis(
    analysis, res[2]["spikes"], delay, flag)
