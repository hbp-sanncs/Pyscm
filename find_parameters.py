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
import pyscm
import pyscm.parameters as par
import pynam.data as data
import pynnless.pynnless_isolated as pynl

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

# Generate BiNAM
mat_in = data.generate(data_params["n_bits_in"], data_params["n_ones_in"],
                       data_params["n_samples"])
mat_out = data.generate(data_params["n_bits_out"], data_params["n_ones_out"],
                        data_params["n_samples"])
print "Data generated!"

# set up simulator
scm = pyscm.SpikeCounterModel(mat_in, mat_out)
sim = pynl.PyNNLessIsolated(sys.argv[1])

# Initial values
weights = {
    "wCH": 0.00,
    "wCA": 0.00,
    "wCSigma": -0.00,
    "wCTExt": 0,  # 0.02,
    "wCTInh": dict["optimise_params"]["wCTInh"],
    "wAbort": 0  # -1.0
}
wCH_min, wCH_max = dict["optimise_params"]["wCH_min"], dict["optimise_params"][
    "wCH_max"]


# Optimisation process
weights["wCH"] = par.optimise_wCH(params, weights, delay, scm, sim,
                                  wCH_min, wCH_max, data_params["n_ones_out"])
print "Calculated 1 of 3"


weights["wCA"] = 2.0 * weights["wCH"]
weights = par.optimise_wCSigma(params, weights, delay, scm, sim, -wCH_max,
                               data_params["n_ones_out"])

print "Calculated 2 of 3"
weights = par.optimise_wCT(params, weights, delay, scm, sim, 0.001, 10)
print "Calculated 3 of 3"


# Write to file
with open("data/optimised_weights.json", 'w') as outfile:
    json.dump(weights, outfile, indent=4)
