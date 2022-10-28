# bistable_helix

Flagella simulations with the generalized IB method and viscoelastic fluids.

There is a sample Makefile that needs to be updated with the IBAMR source and build directories. Alternatively, you can use cmake in an out of source build. For example, you can run cmake in a directory with

cmake /path/to/bistable_helix/source -DIBAMR_ROOT=/path/to/ibamr/install

Note cmake can not be build in source because of the autotools sample Makefile.
