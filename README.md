# Godunov/GRP scheme for Lagrangian/Eulerian hydrodynamics

What is it?
-----------

This is a implementation of fully explict forward Euler scheme for single/multi-fluid Euler equations of motion on Lagrangian/Eulerian coordinate.

How to use it?
-----------

*Programer's Guide* can be found in folder `doc`.

It is made by Doxygen and LaTeX,

Run the following command on the terminal.

```shell
cd src/hydrocode_*
doxygen Doxyfile
cd doc/hydrocode_*/latex
make
cd doc/hydrocode_*/Specification
xelatex Specification.tex
```

Open `doc/hydrocode_*/html/index.html` in a browser to view the specific instructions of this program.

Processing tools
---------

Tecplot(.tec .plt .mcr), ParaView(.vtk), Gmesh(.msh).

Compilers/Interpreters
---------

gcc -std=c99, g++ -std=c++20, Visual Studio 2022(.sln .vcxproj), MATLAB/Octave(.m), Python3(.py), Maple(.mw).

Program library
---------

OpenMP, HDF5 *[NuGet: hdf5-v120-complete]*, GNU Scientific Library (GSL) *[NuGet: gsl-msvc-x86]*, LAPACKE/OpenBLAS.

Debugging tools
---------

Autoconf, Make(.mk), gdb, gprof & gprof2dot *[rely pkg: Graphviz]*, Valgrind *[opt pkg: KCacheGrind, Massif-Visualizer]*, Cppcheck, GCOV(.gcov .gcda .gcno) & LCOV, perf & FlameGraph.

Debug, Run and Release
-----------

```shell
cd src/hydrocode_*
./hydrocode.sh
```

Dependency test
-----------

```shell
cd src/MAKE/Autoconf
autoconf
./configure
```

Licensing
---------

*GNU Lesser General Public License v3.0* or later.

Please see the file called `LICENSE`.

Contacts
--------

If you want more available support for this program, please send an email to [xinlei@cugb.edu.cn](mailto:xinlei@cugb.edu.cn).

Copyright
--------

> Part of this code is modified from the provision of Zhifang Du, Rui Chen & Jian Cheng @ IAPCM.
> Some source codes in the book *C Interfaces and Implementations* designed by David Hanson and 《常用算法程序及》 designed by 徐士良 are used.

Copyright © 2022 Xin Lei.
