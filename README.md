[![Build Status](https://travis-ci.com/danijoo/WHAM.svg?branch=master)](https://travis-ci.com/danijoo/WHAM)

Weighted Histogram Analysis Method (WHAM)
===
This is an fast implementation of the weighted histogram analysis method
written in Rust. It allows the calculation of multidimensional free energy profiles
from umbrella sampling simulations.

Features
---
- Fast, especially for small systems
- Multidimensional
- Unit tested

Usage
---
wham has a convenient command line interface. You can see all options with
```wham -h```:

To run the two dimensional example (simulation of dialanine phi and psi angle):
```bash
wham --max 3.14,3.14 --min -3.14,-3.14 -T 300 --bins 100,100 --cyclic -f example/2d/metadata.dat       
> Supplied WHAM options: Metadata=example/2d/metadata.dat, hist_min=[-3.14, -3.14], hist_max=[3.14, 3.14], bins=[100, 100] verbose=false, tolerance=0.000001, iterations=100000, temperature=300, cyclic=true
> Reading input files.
> 625 windows, 624262 datapoints
> Iteration 10: dF=0.389367172324539
> Iteration 20: dF=0.21450559607810152
(...)
> Iteration 620: dF=0.0000005800554892309461
> Iteration 630: dF=0.00000047424278621817084
> Finished. Dumping final PMF
(... pmf dump ...)

```
After convergence, final bias offsets (F) and the free energy will be dumped to stdout and the output file is written.


The output file contains the free energy and probability for each bin. Probabilities are normalized to sum to P=1.0 and
the smallest free energy is set to 0 (with other free energies based on that).
```
#coord1         coord2       Free Energy Probability
-3.109590    	-3.109590    10.330312   0.000095
-3.046770    	-3.109590    8.907360    0.000168
-2.983950    	-3.109590    7.431969    0.000303
-2.921130    	-3.109590    6.170882    0.000502
-2.858310    	-3.109590    4.982956    0.000809
-2.795490    	-3.109590    3.584741    0.001417
-2.732670    	-3.109590    3.025337    0.001773
(...)
```

Examples
---
The example folder contains input and output files for two simple test systems:

- 1d: Phi torsion angle of dialanine in vaccum
- 2d: Phi and psi torsion angles of the same system


TODO
---
- Multithreading (?)
- Error analysis / bootstrapping
- Better error messages during file I/O
- Autocorrelation
- Replica exchange
- Unit tests for 2d/Nd WHAM

License 
---
WHAM is licensed under the GPLv3 license. Please read the LICENSE file in this
repository for more information.

Parts of this work, especially some perfomance optimizations and the I/O format, are inspired by the
implementation of A. Grossfield (*Grossfield, Alan, "WHAM: the weighted histogram analysis method", http://membrane.urmc.rochester.edu/content/wham*).
