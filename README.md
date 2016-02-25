PySCM
======

No, this is not pySCM, the python simple climate model! This is PySCM!
Instead of weather we do some neural networks.

PySCM stands for Python Spike Counter Model.
It implements a model of associative momeory with spiking neural networks.
If you want see some details: It was developed by Andreas Knoblauch and 
Günther Palm: [1]

## Usage

If you simply want to run the SCM, there are example files ready. Simply type

    ./run.py <SIMULATOR> 

where `SIMULATOR` is the simulator that should be used for execution (see below).
The program reads in the json files in data, namely 'neuron_data.json' containing
all the information about the neurons, the networks and some options. The second
file is the 'optimised_weights.json' file, which contains the synaptic weights of
the different populations.

If you want to simulate with your own network data, then you can change it 
in the 'neuron_data.json'. Since the finding of the weight is a pain, the tool
'find_parameters.py' may help you. It gives you a first idea where your 
parameters should be. Again you can simply type in:
	
	./find_paramters.py <SIMULATOR> 

In Addition, the flag 'Simple_Network' in 'neuron_data.json' allows to switch
between the SCM and a simpler model, which contains only one controlling 
inhibitory neuron. 

At the moment, 'spikey' is not supported. Furthermore, other systems could 
contain a small amount of bugs which reduces the amount of usability of this
model. Nevertheless, it works on 'Nest'.

## Simulators

Possible simulators are:

* (`spikey`)
* `nest`
* `nmmc1`
* `nmpm1`
* `ess`

## Authors

This project was established by Christoph Jenzen and Andreas Stöckel
at Bielefeld University in the [Cognitronics and Sensor Systems Group]
(http://www.ks.cit-ec.uni-bielefeld.de/). This work is 
part of the [Human Brain Project, SP 9](https://www.humanbrainproject.eu/neuromorphic-computing-platform).

## Reference

[1] Andreas Knoblauch, Günther Palm, Pattern separation and synchronization in spiking associative memories and visual areas, Neural Networks, Volume 14, Issues 6–7, 9 July 2001, Pages 763-780, ISSN 0893-6080, http://dx.doi.org/10.1016/S0893-6080(01)00084-3.
(http://www.sciencedirect.com/science/article/pii/S0893608001000843)

## License

This project and all its files are licensed under the
[GPL version 3](http://www.gnu.org/licenses/gpl.txt) unless explicitly stated
differently.


