# Godunov/GRP scheme for Lagrangian/Eulerian hydrodynamics

What is it?
-----------

This is an implementation of a fully explicit forward Euler scheme for single/multi-fluid Euler equations of motion on the Lagrangian/Eulerian coordinate.

How to use it?
-----------

*Programer's Guide* can be found in the folder `doc`.

It is made of Doxygen and LaTeX,

Run the following command on the terminal.

```shell
cd src/hydrocode_*
make doxygen
cd doc/Doxygen/hydrocode_*/latex
make
cd doc/Specification/hydrocode_*
xelatex Specification.tex
```

Open `doc/Doxygen/hydrocode_*/html/index.html` in a browser to view the specific instructions of this program.

- **Online Version**
  > [1D-HydroCODE](https://ximlel.github.io/zh-CN/2022/10/09/Radial-Lag-HydroCODE/1D-HydroCODE/html/index.html) &ensp;
  [2D-HydroCODE](https://ximlel.github.io/zh-CN/2022/10/09/Radial-Lag-HydroCODE/2D-HydroCODE/html/index.html) &ensp;
  [Radial-Lag-HydroCODE](https://ximlel.github.io/zh-CN/2022/10/09/Radial-Lag-HydroCODE/html/index.html)

Processing tools
---------

Tecplot(.tec .plt .mcr), ParaView(.vtk), Gmesh(.msh).

Compilers/Interpreters
---------

gcc/clang/icc -std=c99, g++/clang++ -std=c++20/icpc -std=c++17, Visual Studio 2022(.sln .vcxproj), MATLAB/Octave(.m), Python3(.py), Maple(.mw).

Program library
---------

OpenMP, OpenACC *[rely pkg: NVIDIA HPC SDK (pgcc/nvcc, pgc++/nvc++)]*, HDF5 *[NuGet: hdf5-v120-complete]*, GNU Scientific Library (GSL) *[NuGet: gsl-msvc-x86]*, LAPACKE/OpenBLAS.

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

Clone/Download HTTPS
--------

> <https://gitee.com/ximlel/HydroCODE.git>
> 
> <https://github.com/ximlel/HydroCODE.git>
> 
> <https://osredm.com/p16943850/HydroCODE.git>

Contacts
--------

If you want more available support for this program, please send an email to [xinlei@cugb.edu.cn](mailto:xinlei@cugb.edu.cn).

Copyright
--------

> Part of this code is modified from the provision of Zhifang Du, Rui Chen & Jian Cheng @ IAPCM.
> 
> Some source codes in the book *C Interfaces and Implementations* designed by David Hanson and 《常用算法程序集》 designed by 徐士良 are used.

Copyright © 2022 Xin Lei.
