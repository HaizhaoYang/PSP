Configuring pspBLAS
===================

This document describes how to use the configure script of pspBLAS in various
situations and tune the parameters of the build system.


Recommended procedure
---------------------

We highly recommend the use of build directories. It guarantees you that the
source code is left untouched no matter what and lets you perform concurrent
builds from the exact same source, e.g. with different compilers or levels of
optimisation.

Here is how to build pspBLAS with a build directory:

  1. Go to the top source directory of pspBLAS.
  2. Create a new directory.
  3. Go to the new directory.
  4. Run `../configure` from there to set up the build.
  5. Run `make` to build the libraries.
  6. Run `make check` to run the test suite.
  7. Run `make install` to install the libraries.

Example:

    git clone http://www.e-cam2020.eu:10080/ESL/omm.git
    cd omm/pspBLAS
    ./autogen.sh
    mkdir my-build-with-gcc
    cd my-build-with-gcc
    ../configure --prefix="$HOME/my-libs/pspBLAS" CC="gcc" \
      FC="gfortran" CFLAGS="-O3 -march=native" FCFLAGS="-g -O2"
    make -j4   # Run make using 4 processors
    make check
    make install

If you are contributing to pspBLAS using Git, please note that build
directories starting with "tmp-" will be ignored by Git.

When running configure for the first time, you should type the following:

    ../configure --help

It will provide you with the full list of available options and some useful
information.


Optional features
-----------------

Optional features of pspBLAS are triggered through the `--enable` options of
configure. They can only be set to yes (`--enable-option`) or no
(`--disable-option`). To guarantee the reliability of the build system
behaviour, custom `--enable-option=not_yes_or_no_value` are explicitly
forbidden. All optional features are disabled by default.

The `--enable-debug` option is mainly used by developers and triggers verbose
output from the library at run-time.

The `--enable-timing` option is used for benchmarks and makes pspBLAS provide
timing information at run-time.


MPI parameters
--------------

### Default behaviour ###

If no MPI option is specified on the command line, the build system will check
whether the detected compilers support MPI. If the result is positive, MPI
support will be enabled, or disabled otherwise. If MPI options are specified,
the build system will enable MPI support if all corresponding tests pass, or
stop with an error otherwise, except for `--without-mpi` which disables MPI
support.

Please note that MPI is considered as working only if both the C and Fortran
compilers are MPI-aware.


### Using the --with-mpi option ###

Through the `--with-mpi` option, you can specify the installation prefix
of a working MPI implementation. The build system will automatically look for
MPI-aware compilers in the *bin/* subdirectory of the specified prefix. For
example, if you type:

    ../configure --with-mpi=/usr/local/my_mpi

the build system will look for */usr/local/my_mpi/bin/mpicc* and
*/usr/local/my_mpi/bin/mpif90*, check that they are executable, and set CC and
FC accordingly if they were previously unset.

If the `--with-mpi` option is specified without argument, the build system
will check that the compilers set in *$CC* and *$FC* are actually MPI-aware
and stop with an error if not. On the other hand, if the `--without-mpi`
option is used, MPI support will be completely disabled.

Please note that the `--with-mpi` option conflicts with the use of MPI-related
environment variables.


### Using environment variables ###

Alternatively to the `--with-mpi` option, you can specify MPI-aware compilers
through environment variables, either via CC and FC, or MPICC and MPIFC.
Please note that when using MPICC and MPIFC, both variables must be specified
or the build system will stop with an error, and that using environment
variables conflicts with the use of the `--with-mpi` option.


### Special case: adjustable MPI implementations ###

MPI-aware compilers, such as *mpicc*, *mpicxx*, and *mpif90*, are actually
wrappers which call the serial compilers with all necessary parameters to have
MPI working. A few MPI vendors even allow users to specify the compilers to
use at run-time, instead of hard-coding them in these wrappers when they are
installed. For instance, if you set `CC=gcc` before running *mpicc*, the
wrapper will use GCC to compile your programs, while if you set `CC=icc`, it
will use the Intel compiler.

This is very convenient from the perspective of a system administrator,
because the MPI user interface completely encapsulates the implementation and
hides its peculiarities. A supercomputer can be rearranged and upgraded while
still providing the same interface. This is however a problem for the
Autotools, because they require the C, C++, and Fortran compilers to be
accessed through CC, CXX, and FC, respectively. To use MPI, they require e.g.
*CC* to be set to *mpicc*.

A possible solution to circumvent this issue is to add another layer by
wrapping the wrappers. This is what the build system of pspBLAS does. If a
MPI implementation needs e.g. *CC* to be set to a particular serial compiler,
this definition will be put into a shell script that will call the MPI-aware
wrapper afterwards. This how the shell script schematically looks like:

    CC="gcc"
    export CC
    mpicc $*

The Autotools are then pointed to this script as if it were the real compiler.

You can activate this feature by setting *CC* and *FC* to actual serial
compilers and either use the `--with-mpi` option or the *MPICC* and *MPIFC*
environment variables. The build system will detect this double setting and
activate internal MPI wrapping.


Linear algebra parameters
-------------------------

### Default behaviour ###

Linear algebra is a requirement of pspBLAS. If you do not specify any linear
algebra-related option, or specify `--with-linalg` without arguments, the build
system will look for the include files and libraries using the information made
available by the current development environment. If it fails to detect them,
it will stop with an error.

Please note, that, for obvious reasons, specifying `--without-linalg` is
forbidden and will give an error.


### Using the --with-linalg option ###

Through the `--with-linalg` option, you can specify the installation prefix of
a working linear algebra installation. The build system will automatically look
for include files in the *include/* subdirectory of the specified prefix and
for libraries in the *lib/* subdirectory. For example, if you type:

    ../configure --with-linalg=/usr/local/my_linalg

the build system will look for includes in */usr/local/my_linalg/include/* and
for libraries in */usr/local/my_linalg/lib/*.

With this option, the names of the libraries are set by the build system and
take the MPI configuration into account: if MPI is enabled, the build system
will look for MPI-aware linear algebra libraries, otherwise it will only
look for the serial ones. You can however specify non-standard library names
and compiler-specific parameters through environment variables.

Please note that the `--with-linalg` option conflicts with the use of
linear algebra-related environment variables.


### Using environment variables ###

Alternatively to the `--with-linalg` option, you can specify include flags and
library flags through the *LINALG_INCLUDES* and *LINALG_LIBS* environment
variables. For example:

    ../configure \
      LINALG_INCLUDES="-I/usr/local/my_linalg/include" \
      LINALG_LIBS="-L/usr/local/my_linalg/lib -llapack -lf77blas -lcblas -latlas"

Please note that both variables must be specified or the build system will
stop with an error, and that using environment variables conflicts with the
use of the `--with-linalg` option.


Troubleshooting
---------------

When the configure script stops with an error, it means that the build
environment does not provide sufficient and/or consistent components for the
build of pspBLAS to succeed as you have requested. Here is some advice to
diagnose and possibly fix the issue.

### Diagnosing the issue ###

The first step is to determine what precisely caused the configure script to
stop. To achieve it, the first source of available information is the very
output of configure and the location of the error message within it. The error
message itself should already point quite explicitly to the cause of the error
in most cases. If it occurs at the beginning of the configuration process,
chances are that the options you specified are inconsistent and/or
conflicting. When happening later, it is probable that a required component is
missing or malfunctioning in your build environment.

If the origin of the issue is still not clear, the other main source of
information is the **config.log** file, located in the top build directory,
which tracks down and reports the whole configuration process. Open the file
and look at its contents above the last occurence of *failed*. This will tell
you why the last configuration test performed before the error failed.


### Fixing the issue ###

In many cases, basic instructions on how to fix the issue are provided with
the error message itself. They often take the form of proposed alternatives.
In such a case, select the one which is the most relevant to you. If it
involves the installation of third-party software, please consult your
operating system vendor and/or the component vendor should you need help, in
particular if you do not have the sufficient privileges to install software on
your system.

The most frequent origin of configuration failures is an incomplete and/or
inconsistent build environment. Here is a list of environment variables that
we invite you to check and adjust:

  * LD_LIBRARY_PATH consists in a column-separated list of directories
    containing shared libraries (*.so* files), ordered by search priority; it
    has been used for decades by compilers to dynamically link programs with
    their external dependencies;
  * LIBRARY_PATH is a refinement with respect to LD_LIBRARY_PATH and has the
    same syntax; it is used by modern compilers at link-time, while
    LD_LIBRARY_PATH is only used at run-time; whatever you add to the latter
    should be copied as well into LIBRARY_PATH;
  * PATH also consists in a column-separated list of directories, and
    indicates where to search for executable programs, including compilers.

For example, if linear algebra is installed in */usr/local/my_linalg*, you
should set LD_LIBRARY_PATH and LIBRARY_PATH as follows:

    LD_LIBRARY_PATH="/usr/local/my_linalg/lib:$LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH
    LIBRARY_PATH="/usr/local/my_linalg/lib:$LIBRARY_PATH"
    export LIBRARY_PATH

if you use a Bourne-like shell, or:

    setenv LD_LIBRARY_PATH "/usr/local/my_linalg/lib:$LD_LIBRARY_PATH"
    setenv LIBRARY_PATH "/usr/local/my_linalg/lib:$LIBRARY_PATH"

if you use a C-like shell.

Following these instructions, you will be able to fix nearly 95% of the issues
you will encounter with the configuration process.


### Making your fixes permanent ###

The instructions above will work only within your current session. You will
have to perform these actions again each time you want to build pspBLAS. One
possibility is to save all the shell commands needed for the configuration,
including the parameters of *configure*, into a script. If you already know
for sure that you will use only one version of each external dependency at a
time, a possible alternative is to save the commands above in the
configuration files of your preferred shell (hint: for Bash, the file is named
*.bash_profile*; for Tcsh, it is *.profile*).


### Making your fixes permanent and flexible ###

A popular and widely-used system to flexibly manage build-time and run-time
parameters for all kinds of programs and libraries is the **Environment
Modules** package. If you have used a HPC facility, you are likely familiar
with the `module load` and `module avail` commands. A modern implementation is
available through the Lmod package at: https://github.com/TACC/Lmod

You can also find the original implementation, as well as documentation and
tutorials, at the following address: http://modules.sourceforge.net/

The *Environment Modules* package is available as a standard package for most
Linux distributions and is usually called *lmod* (modern) or
*environment-modules* (legacy). To use it, you just need to create a directory
tree with the names of the components you want to use and store the
corresponding parameters into files named after the component version numbers.

Following is an example for pspBLAS. We will suppose that you have created a
directory named *$HOME/modulefiles/pspBLAS* and have created a file named
*git-version* there.

Here is an example of what you can put in the *git-version* file:

    #%Module1.0
    ## pspBLAS module
    ##
    
    # Package parameters
    set name    "pspBLAS"
    set version "git"
    set build   "gnu_4.9"
    set desc    "$name ($version, $build)"
    set url     "http://www.e-cam2020.eu:10080/ESL/omm.git"
    set root    /home/user/my_libs/$name-$version
    
    proc ModulesHelp { } {
      global name
      global desc
      global url
      puts stderr "This modulefile provides $desc.\n"
      puts stderr "More information about $name can be found at:"
      puts stderr "    $url\n"
    }
    
    module-whatis  "Sets the environment for $desc"
    
    setenv PSP_INCLUDES "-I$root/include"
    setenv PSP_LIBS "-L$root/lib -lpspBLAS"
    
    prepend-path LD_LIBRARY_PATH $root/lib
    prepend-path LIBRARY_PATH $root/lib
    prepend-path PATH $root/bin

Once done, you can issue the following commands:

    module use $HOME/modulefiles
    module avail

You should see *pspBLAS/git-version* appear among the listed packages. In this
case, if you have built and installed pspBLAS in *$HOME/my_libs/pspBLAS-git*,
you will be able to use it by typing `module load pspBLAS` and make it
unavailable by typing `module unload pspBLAS`.

If you want your custom module files to be permanently available, just add the
following to the *$HOME/.modulerc* file:

    module use /home/user/modulefiles

where you replace *user* by your login.


FAQ
---

### How to retrieve configure parameters from the last build? ###

You can find the configure parameters used for the last build of pspBLAS at
the top of the *config.log* file. To use them again, you just need to add
quotes at the appropriate locations.


### Why do I have to bother with C compilers? ###

Fortran is far from being the most standard and widespread language in the
world. On the other hand, there is a wealth of tools available for C, which is
used - among other things - to develop the operating system you are currently
using. Your Fortran compiler is mostly written in C and/or C++. The MPI
implementation you use is written in C and provides Fortran bindings. Several
high-performance linear algebra libraries have a C code base. This is why the
build system needs a C compiler to properly detect many of the components that
will let you successfully build your source code.

