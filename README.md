# Introduction
Algorithms related to the computation of twin smooth integers.

This code was produced as part of my Master Thesis "Finding consecutive smooth integers by solving the Pell equation" in 2022.

I am currently working on an improved version of this (both from a mathematical and a programming standpoint) written entirely in PARI/GP, which will be released soon. This will also include more mathematical background on the methods used.

# How to use
- The `*.sage` files require [sagemath](https://www.sagemath.org/) to be run.
  - `quadratic_orders.sage` contains mostly experiments and some functions used for testing and displaying data related to infrastructure in real quadratic number fields, it wasn't used in the computation of the data.
  - `twin_smooth.sage` contains some more testing functions and functions that were used to compute suitable discriminants, which are also supplied in the `data` directory.
- The file `regulators.gp`is a [GP](https://pari.math.u-bordeaux.fr/) script and was used for computing the regulators of the real quadratic orders for the given discriminants.
- The file `twin_smooth.cc`, written in C++, uses the computational number theory library [LiDIA](https://github.com/mkoeppe/LiDIA), which to my knowledge is the only available computer algebra system, which can utilize suitable compact representations. Unfortunately development has been stopped more than a decade ago (without porting key features to other CAS) and this is the main reason, why this code is suboptimal.
