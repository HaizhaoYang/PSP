Hacking the build system
========================

Using the Autotools
-------------------

The build system of pspBLAS is based on the Autotools, which automate a lot
of common operations related to the build, distribution, and installation, of
software packages, as soon as they comply with a reasonable set of software
programming good practices. By default, the structure of the source code has
to follow the [GNU Coding Standards](https://www.gnu.org/prep/standards/).
Having some familiarity with them will help understanding why some files must
be present within the source tree.

While end-users can safely ignore the Autotools, developers have to know a
minimum about them. For any question, you can always refer to their official
documentation:

  * [Autoconf](https://www.gnu.org/software/autoconf/manual/);
  * [Automake](https://www.gnu.org/software/automake/manual/);
  * [Libtool](https://www.gnu.org/software/libtool/manual/);
  * [GNU M4](https://www.gnu.org/software/m4/manual/);

or to these books and tutorials:

  * [Autotools Mythbuster](https://autotools.io/);
  * [Autotools Guide](http://www.freesoftwaremagazine.com/books/autotools_a_guide_to_autoconf_automake_libtool);
  * [Autotools Tutorial](https://www.lrde.epita.fr/~adl/autotools.html).

However, most of the common development-related tasks are sufficiently
self-explanatory to be performed without consulting these documents.

In the following sections, we describe the tasks a developer has to perform in
different situations, by increasing complexity. We will suppose that you are
already familiar with *make* and the use of *makefiles*. If not, please
consult the [GNU Make Manual](https://www.gnu.org/software/make/manual/) and
practice a little bit before continuing.


Synchronising the build system with the source tree
---------------------------------------------------

The simplest task a developer has to perform is to synchronise the build
system with the source tree every time the structure of the latter has
changed. It requires the execution of a series of predefined operations in a
particular order. In order to facilitate this task, a script called
*autogen.sh* is present at the top of the Git working tree of the package. It
is particularly useful to bootstrap a freshly cloned Git repository. To run
it, just type:

    ./autogen.sh

from the top source directory. Once the *configure* script and the
*Makefile.in* files have been created, later synchronisations can be performed
by running the *autoreconf* convenience script:

    autoreconf -i

The *-i* option tells the script to install or re-install possibly missing
files used by the Autotools. Most of these files are stored in the
*config/gnu/* subdirectory and must always be left untouched.

In very rare circumstances, the synchronisation may fail repeatedly, due to
some inconsistencies between the Autotools files. In such a case, you may run
the *wipeout.sh* script. Be warned though that this script will remove many
files and basically bring back to Git working tree to its pristine state. In
particular, it will remove any file or directory named *tmp*, starting with
*tmp-*, or ending in *.tmp*.

Please note that synchronisation has to occur only when there is a change in
the source tree structure: adding, removing, or renaming, a file or a
directory, or when changing something in the build system itself. It is
usually the last task to perform when one of these changes occurs.


Adding a source file
--------------------

When adding a source file, the build system has to be made aware of the new
file by adding it to the corresponding Makefile.am, usually in the same
directory. Where to add it mainly depends on the nature of the file.

Before going into the details, we would like to attract your attention on the
most common mistakes people do when modifying *Makefile.am* files:

  * using wrong names: in a *Makefile.am*, all names have a meaning for
    Automake, which means that you cannot change them because they look
    "inaesthetical" to you; what feels inadequate to you is actually a
    programming language;
  * forgetting the final '\\' continuation characters in multiline
    descriptions;
  * leaving a final '\\' continuation character on the last line of a
    multiline description.

All these mistakes will likely cause *Automake* to fail. We thus highly
recommend that you use an editor with syntax highlighting when modifying a
*Makefile.am*.

If the file is a ".F90" source, it will go to a *"_SOURCES"* variable. Finding
which one is just a question of common sense, since you already know what this
file is for.

If the file is a ".h" source, it will go to a *"_HEADERS"* variable. If it is
supposed to be installed, it will most likely fit into *"include_HEADERS"*. If
it is purely internal, then it should be registered within *"noinst_HEADERS"*.

Fortran modules - i.e. the *.mod* files - you want to install have to be
manually declared in *"include_DATA"*, since they are not standard files.

If the file can be considered as inert data from the perspective of the build
system, it should go to the *"EXTRA_DIST"* variable. ".c" files included in
other ".c" files enter in this category. Another typical example is a local
*README* file.

Once done with *Makefile.am*, synchronise the build system with the source
tree.

Automatic tracking of Fortran dependencies is a tricky process and is still
under development. Meanwhile, whenever you add a new file, you should update
*src/dependencies.mf* accordingly. When the dependency-generating feature
becomes mature enough, this file will be automatically generated and updated by
the build system.


Adding a library
----------------

Every library defined within pspBLAS is managed by Libtool. Therefore, it
must be declared as *"libname.la"*, where "name" is its distinctive name, in
one of the *"_LTLIBRARIES"* variables. If it is supposed to be installed, it
will have to be decalred in *"lib_LTLIBRARIES"*. If not, it should go into
*"noinst_LTLIBRARIES"*.

Any new library must have a corresponding *"libname_la_SOURCES"* variable,
containing the list of its source files. Please be careful not ot mix
different languages in the same directory, as it creates a lot of portability
issues. In particular, C and Fortran files must go into separate libraries
located in different directories.

If the new library depends on another internal pspBLAS library, it will
require a *"libname_la_LIBADD = libother.la"* declaration as well.

Once done with *Makefile.am*, synchronise the build system with the source
tree.


Adding a program
----------------

Every program must be declared in a *"_PROGRAMS"* variable, either:

  * in *"bin_PROGRAMS"* if it is installed;
  * in *"noinst_PROGRAMS"* if it is not installed;
  * in *"check_PROGRAMS"* if it is a test program executed by `make check`.

Just as libraries, a program named *myprog* requires a *"myprog_SOURCES"*
variable.

If the program depends on internal libraries, it will require a
*"myprog_LDADD = libname.a libother.a"* declaration as well.

Once done with *Makefile.am*, synchronise the build system with the source
tree.


Adding a subdirectory
---------------------

If a new subdirectory only contains data, it can then be managed from the
*Makefile.am* of its parent directory, through the *EXTRA_DIST* variable.

If it contains source files that will be compiled, then it must contain a
*Makefile.am*. To ease the process, you can copy a similar *Makefile.am* to
the new directory and modify it as explained in the previous sections.

Once done with the new *Makefile.am*, you have to add the new directory to the
*"SUBDIRS"* variable of the parent *Makefile.am*.

The following step is to add the *Makefile* (**NOT** the *Makefile.am*) to the
*AC_CONFIG_FILES* list in *configure.ac* (located at the end of the file).

Once done, synchronise the build system with the source tree.


Adding an optional feature trigger
----------------------------------

Optional features are controlled by the `--enable-feature` options of the
*configure* script and linked to the source code through the definition of
*HAVE_FEATURE* preprocessing options. All the relevant code is stored in
*configure.ac*, the source of the *configure* script.

Here is the procedure we recommend to add a new option:

  1. Inform the other developers of what you are up to.
  2. Look in *configure.ac* at the *"enable_debug"* and *"enable_timing"*
     options. They are represented respectively by the *"psp_debug_enable"*
     and *"psp_timing_enabl"* internal build-system variables.
  3. Once you understand enough how things work, duplicate the code associated
     to *"psp_debug_enable"*, renaming it to match your new option.
  4. Synchronise the build system with the source code and run *configure*
     twice, once without your new option and once with it. In both cases,
     look at the contents of *config.h* to confirm that your option is
     correctly defined.
  5. Build pspBLAS with `make` after configuring with your option enabled.
     Check that the corresponding code is actually built.


Adding an external dependency
-----------------------------

If you need to add a new external dependency to pspBLAS, it is necessary to
openly discuss it with the other developers. Unless you have a solid
experience in writing Autotools-based build systems, you will have to do the
corresponding work in tight collaboration with an expert and delegate a
significant part of it. We thus recommend you to define clear specifications
by answering to the following questions beforehand:

  1. Will the new dependency be mandatory or optional?
     Examples: FFTW3 is mandatory, while PFFT is optional.
  2. Will the new dependency be auto-detected or ignored by default?
     Examples: MPI is auto-detected by default, PFFT only if MPI is present.
  3. Is the new dependency available through libraries, or does it use another
     mechanism?
     Examples: FFTW3 and PFFT provide libraries, while MPI provides compilers.
  4. Does the new dependency require already known dependencies?
     Examples: MPI does not, PFFT depends on MPI.
  5. Is the configuration of the new dependency affected by other
     dependencies?
     Examples: PFFT has only one configuration, while FFTW3 provides different
     libraries for the serial and MPI cases.

Adding an external dependency requires the writing of new M4 macros, stored in
the *config/m4/* subdirectory of the source tree. Although M4 is a relatively
simple language, its mastery in the context of the Autotools requires a lot of
practice. This is why at least pair programming will be necessary to achieve
the creation of the new options and related detection mechanisms. Most
probably, the specifications agreed upon will be transmitted to experts and
the details of the implementation left up to them.

