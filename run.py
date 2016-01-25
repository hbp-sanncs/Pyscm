# -*- coding: utf-8 -*-

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
import pynnless as pynl
import sys

if len(sys.argv) != 2:
    print("Usage: " + sys.argv[0] + " <SIMULATOR>")
    sys.exit(1)

n_bits_in = 16
n_bits_out = 16
n_ones_in = 3
n_ones_out = 3
n_samples = entropy.optimal_sample_count(n_bits_in, n_bits_out, n_ones_in, n_ones_out)

mat_in = data.generate(n_bits_in, n_ones_in, n_samples)
mat_out = data.generate(n_bits_in, n_ones_in, n_samples)

scm = pyscm.SpikeCounterModel(mat_in, mat_out)

sim = pynl.PyNNLess(sys.argv[1])
res = sim.run(scm.build())
print res
