[![Build Status](https://travis-ci.com/danijoo/WHAM.svg?branch=master)](https://travis-ci.com/danijoo/WHAM)

Weighted Histogram Analysis Method (WHAM)
===
This is an fast implementation of the weighted histogram analysis method
written in Rust. It allows the calculation of multidimensional free energy profiles
from umbrella sampling simulations. For more details on the method, I suggest Roux, B.
(1995). The calculation of the potential of mena force using computer simulations, CPC, 91(1), 275-282.

Features
---
- Fast, especially for small systems
- Multidimensional
- Error analysis
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
#coord1    coord2    Free Energy    +/-    Probability    +/-
-3.108600    	-3.108600    10.331716    0.000000    0.000095    0.000000
-3.045800    	-3.108600    8.893231    0.000000    0.000170    0.000000
-2.983000    	-3.108600    7.372765    0.000000    0.000312    0.000000
-2.920200    	-3.108600    6.207354    0.000000    0.000498    0.000000
-2.857400    	-3.108600    4.915298    0.000000    0.000836    0.000000
-2.794600    	-3.108600    3.644738    0.000000    0.001392    0.000000
-2.731800    	-3.108600    3.021743    0.000000    0.001787    0.000000
-2.669000    	-3.108600    2.827463    0.000000    0.001932    0.000000
-2.606200    	-3.108600    2.647531    0.000000    0.002076    0.000000
(...)
```

Error analysis
---
WHAM can perform error analysis using the bayesian bootstrapping method. Every simulation window is assumed to be an
individual set of data point. By calculating probabilities N times with randomly assigned weights for each window,
one can estimate the error as standard deviation between the N bootstrapping runs. For more details see
van der Spoel, D. et al. (2010). g_wham—A Free Weighted Histogram Analysis Implementation Including Robust Error and
Autocorrelation Estimates, JCTC, 6(12), 3713-3720.

To perform bayesian bootstrapping in WHAM, use the ```-bt <RUNS>``` flag to perform <RUNS> individual bootstrapping
runs. The error estimates of bin probabilities and free energy will be given as separate column (+/-) in the output file.
If no error analysis is performed, these columns are set to 0.0.

Examples
---
The example folder contains input and output files for two simple test systems:

- 1d: Phi torsion angle of dialanine in vaccum
- 2d: Phi and psi torsion angles of the same system


TODO
---
- Multithreading (?)
- Autocorrelation
- Replica exchange

License 
---
WHAM is licensed under the GPL-3.0 license. Please read the LICENSE file in this
repository for more information.

Parts of this work, especially some perfomance optimizations and the I/O format, are inspired by the
implementation of A. Grossfield (*A. Grossfield, "WHAM: the weighted histogram analysis method", http://membrane.urmc.rochester.edu/content/wham*).
