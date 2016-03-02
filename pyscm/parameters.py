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

import numpy as np
import pynnless.pynnless_isolated as pynl
import pyscm
import pynam.data as data
import json


class WeightOptimisation:
    '''
    This class contains all functions to optimise the weights in the SCM and
    depends on the parameters given in the initialiser.
    It simulates the network and evaluates with the first sample of a
    randomly set BiNAM. For details look into the PyNAM.
    The standard evaluation procedure is to set wCH first, double the value for
    wCA and then optimise for wCSigma. wAbort is set in dependence of wCA.

    TLDR: do_standard_stuff() should work if you gave all the params
    '''

    def __init__(self, data_params, params, delay, optimise_params, simulator,
                 terminating_neurons=1, flag=False):
        '''
        Calculation of the BiNAM and set up of the simulator
        :param data_params: BiNAM data parameters
        :param params: neuron params
        :param delay: delay in the synapses
        :param optimise_params: needs wCH_min, wCH_max, wCTInh = -0.001
        :param simulator: nest, nmmmc1, ...
        :param terminating_neurons: Number of terminating neurons
        :param flag: change between SCM (False) and a single neuron in the
                CS(igma) population
        '''
        # Generate BiNAM
        self.mat_in = data.generate(data_params["n_bits_in"],
                                    data_params["n_ones_in"],
                                    data_params["n_samples"])
        self.mat_out = data.generate(data_params["n_bits_out"],
                                     data_params["n_ones_out"],
                                     data_params["n_samples"])
        print "Data generated!"

        # set up simulator
        self.scm = pyscm.SpikeCounterModel(self.mat_in, self.mat_out)
        self.sim = pynl.PyNNLessIsolated(simulator)

        # Initial values
        self.weights = {
            "wCH": 0.00,
            "wCA": 0.00,
            "wCSigma": -0.00,
            "wCTExt": 0,  # 0.02,
            "wCTInh": optimise_params["wCTInh"],
            "wAbort": 0  # -1.0
        }
        self.wCH_min, self.wCH_max = optimise_params["wCH_min"], \
                                     optimise_params[
                                         "wCH_max"]
        self.params = params
        self.n_ones_out = data_params["n_ones_out"]
        self.delay = delay
        self.terminating_neurons = terminating_neurons
        self.flag = flag
        self.simulator = simulator
        print "Optimisation set up done!"

    # bisection for wCH
    def bisect_wCH(self, min, max):
        self.weights["wCH"] = (max + min) / 2.0
        net, _, _, _ = self.scm.build(params=self.params, weights=self.weights,
                                      delay=self.delay,
                                      terminating_neurons=self.terminating_neurons,
                                      flag=self.flag)
        res = self.sim.run(net, duration=105)
        # n counts the number of ones in output to terminate early if system is to active
        n = 0
        for i in res[0]["spikes"]:
            if i != []:
                n += 1
        if (n >= self.n_ones_out):
            return min, self.weights["wCH"]
        return self.weights["wCH"], max

    # uses the bisection to find minimal wCH
    def optimise_wCH(self):
        wCH_min, wCH_max = self.bisect_wCH(self.wCH_min, self.wCH_max)
        N = 0
        while True:
            wCH_min, wCH_max = self.bisect_wCH(wCH_min, wCH_max)
            N += 1
            if (wCH_max - wCH_min < 0.001):
                break
        self.weights["wCH"] = wCH_max
        print "wCH ready!"
        return self

    # Bisection to find a minimal wCA / maximal wCsigma which keeps the system active
    def find_min(self, min, max, accuracy, str):
        self.weights[str] = (max + min) / 2.0
        net, _, _, _ = self.scm.build(params=self.params, weights=self.weights,
                                      delay=self.delay,
                                      terminating_neurons=self.terminating_neurons,
                                      flag=self.flag)
        res = self.sim.run(net, duration=125)

        # Search for the last spike. If there was no spike np.max throws exception
        try:
            time = np.max(np.max(res[0]["spikes"]))
        except:
            return self.weights[str], max
        if time < 101.8:
            return self.weights[str], max

        # Search in the interval [time-1*delay, time]
        time -= 1 * self.delay
        # Count for number of ones in interval
        n = 0
        for i in res[0]["spikes"]:
            for j in i:
                if j >= time:
                    n += 1
                    # If to active...
                    if (n > accuracy * self.n_ones_out):
                        return min, self.weights[str]
        # Either not active enough or good
        if (n == self.n_ones_out):
            return self.weights[str], self.weights[str]
        return self.weights[str], max

    # Set wCSigma if you want to optimise wCA
    def set_wCSigma(self, weight=-0.01):
        self.weights["wCSigma"] = weight
        return self

    # Search for the optimal wCA, wCSigma should be set before
    def optimise_wCA(self, max=None):
        if (self.weights["wCSigma"] == 0) or (self.weights["wCH"] == 0):
            raise Exception("Set wCSigma and wCH before optimising wCA!")
        if (max == None):
            wCA_max = self.wCH_max
        wCA_min = self.weights["wCH"]
        while True:
            wCA_min, wCA_max = self.find_min(wCA_min, wCA_max, 2, "wCA")
            if (wCA_max - wCA_min < 0.001):
                break
        self.weights["wCA"] = wCA_max
        print "wCA ready!"
        return self

    # Set a standard wCA instead of searching for it, seems to work pretty well
    def set_wCA(self, weight=None):
        if (weight == None):
            self.weights["wCA"] = 2.0 * self.weights["wCH"]
        else:
            self.weights["wCA"] = weight
        return self

    # TODO: Optimisations for the simple model
    # Search for a optimal wCSigma, wCA should be set beforehand
    def optimise_wCSigma(self, wCSigma_min=None, wCSigma_max=-0.001):
        if (self.weights["wCA"] == 0):
            raise Exception("Set wCa before optimising wCSigma!")
        if (wCSigma_min == None):
            wCSigma_min = - self.weights["wCA"]
        while True:
            wCSigma_min, wCSigma_max = self.find_min(wCSigma_min, wCSigma_max,
                                                     2, "wCSigma")
            print wCSigma_min, wCSigma_max
            if (np.abs(wCSigma_min - wCSigma_max) < 0.001):
                break
        self.weights["wCSigma"] = wCSigma_min
        print "wCSigma ready!"
        return self

    # Bisection for wCT, differs a bit from the others.
    # It assures that the system spikes at least the min_spike_count
    # + that the termination neuron spikes at least one time
    def bisect_wCT(self, min, max, min_spike_count=5):
        self.weights["wCTExt"] = (max + min) / 2.0
        net, _, _, _ = self.scm.build(params=self.params, weights=self.weights,
                                      delay=self.delay,
                                      terminating_neurons=self.terminating_neurons,
                                      flag=self.flag)
        res = self.sim.run(net, duration=125)

        # If termination neuron does not spike, increase weight
        if (len(res[2]["spikes"][0]) == 0):
            return self.weights["wCTExt"], max
        # Estimate the number of spikes of an active neuron
        try:
            max_spikes = np.max([len(i) for i in res[0]["spikes"]])
        except:
            return min, self.weights["wCT"]
        # If the number of counts is achieved, set to -1 to stop the process
        if (max_spikes > min_spike_count):
            return -1, 0
        return min, self.weights["wCTExt"]

    # Set wAbort to standard value
    def set_wAbort(self, weight=None):
        if (weight == None):
            self.weights["wAbort"] = -1.5 * self.weights["wCA"] / float(
                self.terminating_neurons)
        else:
            self.weights["wAbort"] = weight
        return self

    # Optimising wCT with the bisect_wCT
    def optimise_wCT(self, wCT_min=0.001, wCT_max=None, minimal_spike_count=5):
        if (self.weights["wAbort"] == 0):
            raise Exception(
                "It makes no sense to optimise wCT before setting wAbort!")
        if (wCT_max == None):
            wCT_max = self.weights["wCA"]
        while True:
            wCT_min, wCT_max = self.bisect_wCT(wCT_min, wCT_max,
                                               minimal_spike_count)
            print wCT_min, wCT_max
            if (np.abs(wCT_min - wCT_max) < 0.001):
                break
            if (wCT_min < 0):
                break
        print "wCT ready!"
        return self

    # Make the standard optimisation
    def do_standard_stuff(self):
        return self.optimise_wCH().set_wCA().optimise_wCSigma().set_wAbort().optimise_wCT()

    # Write all weights to a json file, which can be read in from the scm
    def write_to_file(self, file="data/optimised_weights.json"):
        with open(file, 'w') as outfile:
            json.dump(self.weights, outfile, indent=4)
