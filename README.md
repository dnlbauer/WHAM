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

License 
---
WHAM is licensed under the GPLv3 license. Please read the LICENSE file in this
repository for more information.

Parts of this work, especially some perfomance optimizations and the I/O format, are inspired by the
implementation of A. Grossfield (*Grossfield, Alan, "WHAM: the weighted histogram analysis method", http://membrane.urmc.rochester.edu/content/wham*).
