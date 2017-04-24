Contributing to pspBLAS
=======================

pspBLAS follows the GNU Coding Standards as much as permitted by the
experience of its developers. If you are planning to contribute, we thus
recommend you to consult the official [GNU Coding Standards
homepage](https://www.gnu.org/prep/standards/) before doing so.


Programming languages
---------------------

The core of pspBLAS is written in Fortran. All other languages must be supported
through bindings.

Since we use the Autotools to build, install, and distribute the code, all C,
C++, or Fortran files, must be stored in language-specific subdirectories.
Mixes should not happen in the same directory, and even less in the same
library.


Optional features
-----------------

Optional features are triggered either through *HAVE_FEATURE* or *PSP_FEATURE*
preprocessing options, where *FEATURE* is the name of the feature. All these
preprocessing options are handled through *config.h*, with the exception of
the Autotools-reserved *HAVE_CONFIG_H* option.

Within the source code, access to the optional features is granted in a
boolean way, i.e.:

    #if defined HAVE_FEATURE
    ... do something optional ...
    #endif

This allows keeping the build system as simple and maintainable as possible.


Keeping the build system updated
--------------------------------

When changing the contents of a file, the build system does not need to be
updated. However, when adding, removing, or renaming files, you have to inform
the build system of what happened. The corresponding procedures are described
in [Hacking the build system](hacking-the-build-system.html).


Testing
-------

When adding new functions, it is extremely important to provide simple test
programs, aka "unit tests", to check whether these functions are performing as
they should. If you do not feel comfortable with integrating these tests with
the build system, please notify the other developers.

Systematically writing unit tests is not only essential to maintain the
overall quality of the code. It will greatly help efficiently design and
structure your contributions, since you will always have a concrete use
example at hand to feed your thoughts.

