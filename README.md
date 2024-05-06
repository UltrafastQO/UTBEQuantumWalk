# UTBEQuantumWalk
Model of an ultrafast time-bin encoding (UTBE) quantum walk using Strawberryfields (https://strawberryfields.readthedocs.io/).

Requires packages numpy and strawberryfields.

Reference: K. Fenwick et al. https://arxiv.org/abs/2404.02238

## How it works

The UTBE quantum walk being simulated is a sequence of aBBO crystals whose angle alternates between 0 deg and 45. 

As an example, consider the case of a 2-step walk. Using Strawberryfields, the code implements the circuit shown below:

<p align="center">
    <img src="graph.svg"/>
</p>

The input state is assumed to be a two-mode squeezed vacuum state in modes 0 and 1, and a coherent state in the mode 2. Mode 0 is a the heralding mode and does not participate in the walk.
