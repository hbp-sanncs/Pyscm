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

import pynam
import pynnless as pynl


class SpikeCounterModel:
    def __init__(self, mat_in, mat_out, params={}, _type=pynl.TYPE_IF_COND_EXP):
        # Copy the given reference to the input and output data vectors
        self.mat_in = mat_in
        self.mat_out = mat_out
        self.params = params
        self._type = _type

        # Train the matrices CH and CA
        self.n_bits_in = mat_in.shape[1]
        self.n_bits_out = mat_out.shape[1]
        self.mat_CH = pynam.BiNAM().train_matrix(mat_in, mat_out)
        self.mat_CA = pynam.BiNAM().train_matrix(mat_out, mat_out)

    def build(self):
        net = pynl.Network()

        # Create the individual cell assemblies
        pop_C = pynl.Population(count=self.n_bits_out, _type=self._type, params=self.params,
                                     record=[pynl.SIG_SPIKES])
        pop_CS = pynl.Population(pop_C)
        pop_CT = pynl.Population(count=1, _type=self._type, params=self.params, record=[pynl.SIG_SPIKES])
        pop_source = pynl.Population(count=self.n_bits_in, _type=pynl.TYPE_SOURCE)

        # Create the connections
        def connections_from_matrix(mat, pSrc, pTar, w, delay=0.0):
            connections = []
            for i in xrange(mat.shape[0]):
                for j in xrange(mat.shape[1]):
                    if mat[i, j] != 0:
                        connections.append(((pSrc, i), (pTar, j), w, delay))
            return connections

        def connections_all_to_all(m, n, pSrc, pTar, w, delay=0.0):
            connections = map(lambda _: [], xrange(m * n))
            for i in xrange(m):
                for j in xrange(n):
                    connections[i * n + j] = (((pSrc, i), (pTar, j), w, delay))
            return connections

        wExc = 1.0
        wInh = -1.0
        iC = 0
        iCS = 1
        iCT = 2
        iSource = 3
        connections = (
            connections_from_matrix(self.mat_CH, iSource, iC, wExc) +
            connections_from_matrix(self.mat_CH, iSource, iCS, wExc) +
            connections_from_matrix(self.mat_CH, iC, iCS, wExc) +
            connections_from_matrix(self.mat_CA, iC, iC, wExc) +
            connections_all_to_all(self.n_bits_out, self.n_bits_out, iCS, iCS, wInh) +
            connections_all_to_all(self.n_bits_out, 1, iCS, iCT, wInh) +
            connections_all_to_all(1, self.n_bits_out, iCT, iCS, wInh) +
            connections_all_to_all(self.n_bits_out, 1, iC, iCT, wExc) +
            connections_all_to_all(1, self.n_bits_out, iCT, iC, wInh)
        )

        return pynl.Network(populations = [pop_C, pop_CS, pop_CT, pop_source], connections=connections)
