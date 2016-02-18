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
import pynam.network as netw
import pyscm.parameters as par
import sys

if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <SIMULATOR>")
    sys.exit(1)

# params = {
#     "v_rest": -70,
#     "v_thresh": -1,
#     "e_rev_E": 0,
#     "e_rev_I": -100,
#     "tau_syn_E": 5,
#     "tau_syn_I": 5,
#     "tau_m": 5
# }
# delay = 0.1
#
# n_bits_in = 100
# n_bits_out = 100
# n_ones_in = 4
# n_ones_out = 4
# n_samples = 200
# data_params = netw.DataParameters(n_bits_in=n_bits_in, n_bits_out=n_bits_out,
#                                   n_ones_in=n_ones_in, n_ones_out=n_ones_out,
#                                   n_samples=n_samples)
# dict={"neuron_params":params,"data_params" : data_params}
# with open("data/neuron_data.json", 'w') as outfile:
#     json.dump(dict, outfile, indent=4)

with open("data/neuron_data.json", 'r') as outfile:
    dict = json.load(outfile)

data_params = dict["data_params"]
params = dict["neuron_params"]
delay = dict["delay"]

mat_in = data.generate(data_params["n_bits_in"], data_params["n_ones_in"],
                       data_params["n_samples"])
mat_out = data.generate(data_params["n_bits_out"], data_params["n_ones_out"],
                        data_params["n_samples"])
print "Data generated!"
scm = pyscm.SpikeCounterModel(mat_in, mat_out)

simulator = sys.argv[1]

weights = {
    "wCH": 0.00,
    "wCA": 0.00,
    "wCSigma": -0.00,
    "wCTExt": 0,  # 0.02,
    "wCTInh": -0.001,
    "wAbort": 0  # -1.0
}
wCH_min, wCH_max = 0.5, 10
wCA_min, wCA_max = 0.5, 10
wCSigma_min, wCSigma_max = -1.0, -0.001

weights["wCH"] = par.optimise_wCH(params, weights, delay, scm, simulator,
                                  wCH_min, wCH_max, data_params["n_ones_out"])
print weights["wCH"]

weights["wCA"] = par.Binam_wCA(params, weights, delay, scm, simulator, wCA_min,
                               wCA_max, data_params["n_ones_out"])
print weights["wCA"]

weights = par.optimise_wCA(params, weights, delay, scm, simulator,
                           data_params["n_ones_out"])

with open("data/optimised_weights.json", 'w') as outfile:
    json.dump(weights, outfile, indent=4)
