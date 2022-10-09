AC_INIT([hydrocode], [0.4], [xinlei@cugb.edu.cn])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_OPENMP
#AC_OPENACC
AC_CHECK_PROGS([MAKE],       [$MAKE make])
AC_CHECK_PROGS([DOXYGEN],    [doxygen doxywizard])
AC_CHECK_PROGS([LATEX],      [pdflatex])
AC_CHECK_PROGS([XELATEX],    [xelatex])
AC_CHECK_PROGS([OCTAVE],     [octave])
AC_CHECK_PROGS([PYTHON],     [python3])
AC_CHECK_PROGS([GDB],        [gdb])
AC_CHECK_PROGS([GPROF],      [gprof])
AC_CHECK_PROGS([GPROF2DOT],  [gprof2dot])
AC_CHECK_PROGS([GCOV],       [gcov])
AC_CHECK_PROGS([LCOV],       [lcov])
AC_CHECK_PROGS([PERF],       [perf])
#AC_CHECK_PROGS([FLAMEGRAPH,  [flamegraph.pl stackcollapse-perf.pl])
AC_CHECK_PROGS([VALGRIND],   [valgrind])
AC_CHECK_PROGS([KCacheGrind],[kcachegrind])
#AC_CHECK_PROGS([Massif-Visualizer], [massif-visualizer])
AC_CHECK_PROGS([CPPCHECK],   [cppcheck])
AC_CHECK_PROGS([PARAVIEW],   [paraview])
AC_CHECK_PROGS([GRAPHVIZ],   [dot])
AC_CHECK_PROGS([GMESH],      [gmsh])

# Checks for libraries.
AC_CHECK_LIB([hdf5],      [H5Fopen])
AC_CHECK_LIB([gcov],      [__gcov_dump])
AC_CHECK_LIB([pthread],   [sem_init])
AC_CHECK_LIB([m],         [exp])
AC_CHECK_LIB([lapack],    [dgelsy_])
AC_CHECK_LIB([gslcblas],  [main])
AC_CHECK_LIB([gsl],       [gsl_expm1])

AC_OUTPUT
