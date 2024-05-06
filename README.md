# UTBEQuantumWalk
Model of an ultrafast time-bin encoding (UTBE) quantum walk using Strawberryfields (https://strawberryfields.readthedocs.io/).

For details, see: K. Fenwick et al. https://arxiv.org/abs/2404.02238

Requires packages numpy and strawberryfields.

## How it works

As an example, consider a 2-step quantum walk. The code implements the circuit shown below in Strawberryfields.

<img src="graph.svg"/>

The input state is a two-mode squeezed vacuum state in modes 0 and 1, and a coherent state in the mode 2. Mode 0 is a the heralding mode and does not participate in the walk. The transformation U_BS is a 50:50 beam splitter while U_SWAP is a SWAP operation (i.e. a beam splitter with unit reflectivity). To ensure light does not propagate backwards in time, the SWAP operations are applied on the later time bins first.
