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

import pynnless as pynl
import numpy as np


# bisection for wCH
def bisect_wCH(params, weights, delay, scm, simulator, min, max, n_ones_out):
    weights["wCH"] = (max + min) / 2.0
    net, _, _, _ = scm.build(params=params, weights=weights,
                             delay=delay)
    sim = pynl.PyNNLess(simulator)
    res = sim.run(net, duration=105)

    # if np.asarray(res[0]["spikes"]).size:
    #    return min, weights["wCH"]
    n = 0
    for i in res[0]["spikes"]:
        if i != []:
            n += 1
    if (n >= n_ones_out):
        return min, weights["wCH"]
    return weights["wCH"], max


def optimise_wCH(params, weights, delay, scm, simulator, min, max, n_ones_out):
    wCH_min, wCH_max = bisect_wCH(params, weights, delay, scm, simulator, min,
                                  max, n_ones_out)
    N = 0
    while True:
        wCH_min, wCH_max = bisect_wCH(params, weights, delay, scm, simulator,
                                      wCH_min, wCH_max, n_ones_out)
        N += 1
        if (wCH_max - wCH_min < 0.001):
            break
    print N
    return wCH_max


def find_min_wCA(params, weights, delay, scm, simulator, min, max, n_ones_out):
    weights["wCA"] = (max + min) / 2.0
    net, _, _, _ = scm.build(params=params, weights=weights,
                             delay=delay)
    sim = pynl.PyNNLess(simulator)
    res = sim.run(net, duration=125)

    # if np.asarray(res[0]["spikes"]).size:
    #    return min, weights["wCH"]
    n = 0
    try:
        time = np.max(np.max(res[0]["spikes"]))
    except:
        return weights["wCA"], max
    if time < 120:
        return weights["wCA"], max
    time -= 2 * delay
    for i in res[0]["spikes"]:
        for j in i:
            if j >= time:
                n += 1
                if (n >= 5 * n_ones_out):
                    return min, weights["wCA"]
    return weights["wCA"], max


def Binam_wCA(params, weights, delay, scm, simulator, min, max, n_ones_out):
    N = 1
    wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, simulator, min,
                                    max, n_ones_out)
    while True:
        wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, simulator,
                                        wCA_min, wCA_max, n_ones_out)
        N += 1
        if (wCA_max - wCA_min < 0.001):
            break
    print N
    print wCA_min, wCA_max
    return wCA_max


def optimise_wCA(params, weights, delay, scm, simulator, n_ones_out):
    # w_CA_fallback = weights["wCA"]
    weights["wCSigma"] = -0.1
    wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, simulator,
                                    weights["wCA"], 30, n_ones_out)
    N = 0
    while True:
        wCA_min, wCA_max = find_min_wCA(params, weights, delay, scm, simulator,
                                        wCA_min, wCA_max, n_ones_out)
        N += 1
        if (wCA_max - wCA_min < 0.01):
            break
    print N
    weights["wCA"] = wCA_max
    return weights

# def optimise_wCSigma(weights, min, max, n_ones_out):
#     weights["wCSigma"] = (max + min) / 2.0
#     net, _, _, _ = scm.build(params=params, weights=weights,
#                              delay=delay)
#     sim = pynl.PyNNLess(sys.argv[1])
#     res = sim.run(net, duration=125)
#
#     # if np.asarray(res[0]["spikes"]).size:
#     #    return min, weights["wCH"]
#     n = 0
#     for i in res[0]["spikes"]:
#         for j in i:
#             if np.amax(j) >= 120:
#                 n += 1
#                 if (n > n_ones_out):
#                     return min, weights["wCSigma"], False
#     if (n < n_ones_out):
#         return weights["wCSigma"], max, False
#     return weights["wCSigma"], max, True
#
#
# N = 0
# stop = 0
# while True:
#     wCSigma_min, wCSigma_max, test = optimise_wCSigma(weights, wCSigma_min,
#                                                       wCSigma_max, n_ones_out)
#     print wCSigma_min, wCSigma_max
#     if (test == True):
#         break
#     N += 1
#     if (wCSigma_max - wCSigma_min < 0.001):
#         if (stop > 20):
#             print "No optimal parameter could be found"
#             break
#         weights["wCA"] += 0.1
#         stop += 1
#         wCSigma_min, wCSigma_max = -1.0, -0.001
#         # TODO values from the beginning
#         continue
# print N, stop
# print wCSigma_min, wCSigma_max
# weights["wCSigma"] = wCSigma_min
