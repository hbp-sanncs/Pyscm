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


# bisection for wCH
def bisect_wCH(params, weights, delay, scm, sim, min, max, n_ones_out):
    weights["wCH"] = (max + min) / 2.0
    net, _, _, _ = scm.build(params=params, weights=weights,
                             delay=delay)
    res = sim.run(net, duration=105)
    # n counts the number of ones in output to terminate early if system is to active
    n = 0
    for i in res[0]["spikes"]:
        if i != []:
            n += 1
    if (n >= n_ones_out):
        return min, weights["wCH"]
    return weights["wCH"], max


# uses the bisection to find minimal wCH
def optimise_wCH(params, weights, delay, scm, sim, min, max, n_ones_out):
    wCH_min, wCH_max = bisect_wCH(params, weights, delay, scm, sim, min,
                                  max, n_ones_out)
    N = 0
    while True:
        wCH_min, wCH_max = bisect_wCH(params, weights, delay, scm, sim,
                                      wCH_min, wCH_max, n_ones_out)
        N += 1
        if (wCH_max - wCH_min < 0.001):
            break
    return wCH_max


# Bisection to find a minimal wCA which keeps the system active with given weights
def find_min_wCA(params, weights, delay, scm, sim, min, max, n_ones_out,
                 accuracy):
    weights["wCA"] = (max + min) / 2.0
    net, _, _, _ = scm.build(params=params, weights=weights,
                             delay=delay)
    res = sim.run(net, duration=125)

    # Search for the last spike. If there was no spike np.max throws exception
    try:
        time = np.max(np.max(res[0]["spikes"]))
    except:
        return weights["wCA"], max
    if time < 120:
        return weights["wCA"], max

    # Search in the interval [time-2*delay, time]
    time -= 2 * delay
    # Count for number of ones in interval
    n = 0
    for i in res[0]["spikes"]:
        for j in i:
            if j >= time:
                n += 1
                # If to active...
                if (n >= accuracy * n_ones_out):
                    return min, weights["wCA"]
    # Either not active enough or good
    return weights["wCA"], max


# Search for the optimal wCA for the non-inhibitory part of the SCM
def Binam_wCA(params, weights, delay, scm, sim, min, max, n_ones_out):
    N = 1
    wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, sim, min,
                                    max, n_ones_out, 3)
    while True:
        wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, sim,
                                        wCA_min, wCA_max, n_ones_out, 3)
        N += 1
        if (wCA_max - wCA_min < 0.001):
            break
    return wCA_max


# Search for the optimal wCA with inhibitory part
def optimise_wCA(params, weights, delay, scm, sim, wCA_max, n_ones_out):
    # For now the strategy is to set wCSigma to a value and make wCA bigger step by step
    wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, sim,
                                    weights["wCA"], wCA_max, n_ones_out, 2)
    N = 0
    while True:
        wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, sim,
                                        wCA_min, wCA_max, n_ones_out, 2)
        N += 1
        if (wCA_max - wCA_min < 0.001):
            break
    weights["wCA"] = wCA_max
    return weights
