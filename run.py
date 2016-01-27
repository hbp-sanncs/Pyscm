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
import pyscm
import pynam.data as data
import pynam.entropy as entropy
import pynam.network as netw
import pynnless as pynl
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <SIMULATOR>")
    sys.exit(1)

params = {
    "v_rest": -50,
    "v_thresh": -1,
    "e_rev_E": 0,
    "e_rev_I": -100,
    "tau_syn_E": 5,
    "tau_syn_I": 5,
    "tau_m": 5
}

# weights = {
#     "wCH": 6.27,
#     "wCA": 18.0,
#     "wCSigma": -0.1,
#     "wCTExt": 1.046,
#     "wCTInh": -0.01,
#     "wAbort": -2000.0
# }
# n_bits_in = 4
# n_bits_out = 20
# n_ones_in = 2
# n_ones_out = 4
# n_samples = 6

weights = {
    "wCH": 3.25,
    "wCA": 12.6,
    "wCSigma": -0.1,
    "wCTExt": 0.5,
    "wCTInh": -0.002,
    "wAbort": -2000.0
}
n_bits_in = 100
n_bits_out = 100
n_ones_in = 4
n_ones_out = 4
n_samples = 200

mat_in = data.generate(n_bits_in, n_ones_in, n_samples)
mat_out = data.generate(n_bits_out, n_ones_out, n_samples)

scm = pyscm.SpikeCounterModel(mat_in, mat_out)

sim = pynl.PyNNLess(sys.argv[1])
net, input_indices, _, input_times = scm.build(params=params, weights=weights)
res = sim.run(net)

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

output_times, output_indices = netw.NetworkInstance.match_static(input_times,
                                                                 input_indices,
                                                                 res[0]["spikes"])

# Missing parameters?
analysis = netw.NetworkAnalysis(input_times=input_times, input_indices=input_indices,
                                output_times=output_times,
                                output_indices=output_indices,
                                mat_in=mat_in, mat_out=mat_out)

I, mat_out_res, errs = analysis.calculate_storage_capactiy(netw.OutputParameters(burst_size=10))

print "Information:", I
print mat_out_res
print errs

plt.show()
