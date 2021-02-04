## Newton solver

In the first class, you learned about Newton's method, which can be suitably
described in about five bullet points. This representation, of course,
abstracts away the problem of computing derivatives (Jacobians) and factorizing
a sparse matrix. A typical implementation of Newton's method doesn't write
these subroutines from scratch, but interfacing with high-quality existing codes
for these purposes can be a nontrivial exercise in itself.

Here I've implemented Newton's method using the Ampl Solver Library (ASL) for
derivatives and HSL MA28 for factorizing the sparse unsymmetric Jacobian.
The source code needs to be compiled against these libraries, then can be called
using the Ampl interface (e.g. from the command line, Ampl, or Pyomo).

The following instructions assume you are using a Linux machine, and have been
tested on `linux.andrew.cmu.edu`, which you should all have access to via SSH
and your Andrew accounts.


### Dependencies

This code relies on ASL and MA28.

ASL can be obtained from the following link:
http://www.netlib.org/ampl/solvers.tgz

while MA28 must be requested from the following link:
http://www.hsl.rl.ac.uk/download/MA28/1.0.0/a/

Once the above tarballs have been obtained, they should be unpacked into
this directory. You should now have subdirectories `solvers` and `ma28-1.0.0`.

To compile ASL, navigate into `solvers` and run
```sh
$ ./configurehere
$ make
```
You should now have `amplsolver.a` in this directory.

To compile MA28, navigate into `ma28-1.0.0` and run
```sh
$ ./configure
$ make
```
You should now have `libma28.a` in `ma28-1.0.0/src`.


### Compile the solver

With the above dependencies compiled in their default locations,
you should be able to compile the solver by navigating into this
directory's `src` subdirectory and runnning
```sh
$ make
```

### Examples

Right now I have included one example of calling this solver
from Python. You will need to have Pyomo installed.
The example script assumes the solver executable is located
in its parent directory (the directory of this README file),
and can be run as is:
```sh
$ python simple.py
```
