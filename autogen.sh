#!/bin/sh
#
# Copyright (C) 2016 Yann Pouillon <notifications@materialsevolution.es>
#
# This file is part of pspBLAS.
#
# pspBLAS is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, version 3 of the License, or (at your option) any later
# version.
#
# pspBLAS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pspBLAS.  If not, see <http://www.gnu.org/licenses/> or write to
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301  USA.

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "src/pspBLAS.F90"; then
  echo "This is not a pspBLAS source tree - aborting now"
  exit 1
fi

# Create possibly missing directories
mkdir -p config/gnu config/m4

# Generate M4 macros
#echo "Generating M4 macros..."
#echo "done."

# Generate makefiles
#echo "Generating makefiles..."
#echo "done."

# Generate M4 includes
echo "Generating aclocal.m4..."
aclocal -I config/m4
echo "done."

# Generate configure auxiliary files
echo "Generating config.h.in..."
autoheader
echo "done."

# Generate configure
echo "Generating configure script..."
autoconf
echo "done."

# Generate libtool scripts
echo "Generating libtool scripts..."
my_libtoolize="libtoolize"
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  my_libtoolize="glibtoolize"
fi
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  echo "Error: could not find a working version of libtoolize" >&2
  exit 1
fi
${my_libtoolize} --automake --copy --force
echo "done."

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "Generating Makefile.in for each directory..."
automake --add-missing --copy
echo "done."
