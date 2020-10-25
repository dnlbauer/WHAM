[![Build Status](https://travis-ci.com/danijoo/WHAM.svg?branch=master)](https://travis-ci.com/danijoo/WHAM) [![crates.io](https://img.shields.io/badge/crates.io-orange.svg?longCache=true)](https://www.crates.io/crates/wham) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1488597.svg)](https://doi.org/10.5281/zenodo.1488597)



Weighted Histogram Analysis Method (WHAM)
===
This is an fast implementation of the weighted histogram analysis method
written in Rust. It allows the calculation of multidimensional free energy profiles
from umbrella sampling simulations. For more details on the method, I suggest *Roux, B.
(1995). The calculation of the potential of mean force using computer simulations, CPC, 91(1), 275-282.*

Features
---
- Fast, especially for small systems
- Multithreaded
- Multidimensional
- Autocorrelation to remove correlated samples
- Error analysis
- Unit tested

Installation
---
WHAM requires the GSL library to be installed: 
```bash
# on debian/ubuntu:
sudo apt-get install libgsl0-dev
```

Installation from source via cargo:
```bash
# cargo installation
curl -sSf https://static.rust-lang.org/rustup.sh | sh

cargo install wham
```

Usage
---
wham has a convenient command line interface. You can see all options with
```wham -h```:

```
wham 0.9.9
D. Bauer <bauer@bio.tu-darmstadt.de>
wham is a fast implementation of the weighted histogram analysis method (WHAM) written in Rust. It currently supports
potential of mean force (PMF) calculations in multiple dimensions at constant temperature.

Metadata file format:
    /path/to/timeseries_file1  x_1  x_2  x_N  fc_1  fc_2  fc_N
    /path/to/timeseries_file2  x_1  x_2  x_N  fc_1  fc_2  fc_N
    /path/to/timeseries_file3  x_1  x_2  x_N  fc_1  fc_2  fc_N
The first column is a path to a timeseries file _relative_ to the metadata file (see below). This is followed by the
position of the umbrella potential x in N dimensions and the force constant fc in each dimension. Lines starting with a
# are treated as comments and will not be parsed.

Timeseries file format:
    time  x_1  x_2  x_N
    time  x_1  x_2  x_N
    time  x_1  x_2  x_N
The first column will be ignored and is followed by N reaction coordinates x.

Shipped under the GPLv3 license.

USAGE:
    wham [FLAGS] [OPTIONS] --bins <BINS> --max <HIST_MAX> --file <METADATA> --min <HIST_MIN> --temperature <temperature>

FLAGS:
    -c, --cyclic     For periodic reaction coordinates. If this is set, the first and last coordinate bin in each
                     dimension are treated as neighbors for the bias calculation.
    -h, --help       Prints help information
    -g, --uncorr     Estimates statistical inefficiency of each timeseries via autocorrelation and removes correlated
                     samples (default is off).
    -V, --version    Prints version information
    -v, --verbose    Enables verbose output.

OPTIONS:
    -b, --bins <BINS>                  Number of histogram bins (comma separated).
        --bt <bootstrap>               Number of bayesian bootstrapping runs for error analysis by assigning random
                                       weights (defaults to 0).
        --seed <bootstrap_seed>        Random seed for bootstrapping runs.
        --end <end>                    Skip rows in timeseries with an index larger than this value (defaults to 1e+20)
    -i, --iterations <ITERATIONS>      Stop WHAM after this many iterations without convergence (defaults to 100,000).
        --max <HIST_MAX>               Histogram maxima (comma separated). Also accepts "pi".
    -f, --file <METADATA>              Path to the metadata file.
        --min <HIST_MIN>               Histogram minima (comma separated for multiple dimensions). Also accepts "pi".
    -o, --output <output>              Free energy output file (defaults to wham.out).
        --start <start>                Skip rows in timeseries with an index smaller than this value (defaults to 0)
    -T, --temperature <temperature>    WHAM temperature in Kelvin.
    -t, --tolerance <TOLERANCE>        Abortion criteria for WHAM calculation. WHAM stops if abs(F_new - F_old) <
                                       tolerance (defaults to 0.000001).
```

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
*Van der Spoel, D. et al. (2010). g_wham—A Free Weighted Histogram Analysis Implementation Including Robust Error and
Autocorrelation Estimates, JCTC, 6(12), 3713-3720*.

To perform bayesian bootstrapping in WHAM, use the ```-bt <RUNS>``` flag to perform <RUNS> individual bootstrapping
runs. The error estimates of bin probabilities and free energy will be given as standard error (SE) in a 
separate column (+/-) in the output file. If no error analysis is performed, these columns are set to 0.0.

Autocorrelation analysis
---
With the ```--uncorr``` flag, WHAM calculates the autocorrelation time ```tau``` for all timeseries and all collective
variables. Timeseries are then filtered based on their highest autocorrelation time to remove correlated samples from
the dataset. This reduces the number of data points but can improve the accuracy of the result.

For filtering, the statistical inefficiency `g` is calculated: ```g = 1 + 2*tau```, and only every `g`th element of the
timeseries is used for unbiasing. A more detailed description of the method can be found in
*Chodera, J.D. et al. (2007). Use of the weighted histogram analysis method for the analysis of simulated and parallel
tempering simulations, JCTC 3(1):26-41*


Examples
---
The example folder contains input and output files for two simple test systems:

- 1d_cyclic: Phi torsion angle of dialanine in vaccum
- 2d_cyclic: Phi and psi torsion angles of the same system


TODO
---
- Replica exchange

License & Citing
---
WHAM is licensed under the GPL-3.0 license. Please read the LICENSE file in this
repository for more information.

There's no publication for this WHAM implementation. However, there is a citeabe DOI. If you use this software for your work, please consider citing it: *Bauer, D, WHAM - An efficient weighted histogram analysis implementation written in Rust, Zenodo.  https://doi.org/10.5281/zenodo.1488597*

Parts of this work, especially some perfomance optimizations and the I/O format, are inspired by the
implementation of A. Grossfield (*Grossfield, A, WHAM: the weighted histogram analysis method, http://membrane.urmc.rochester.edu/content/wham*).
