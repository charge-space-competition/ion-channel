# ----------------------------------------------------------------------
# This source file is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

#
# Intel Maths Kernel: BLAS/LAPACK library 
#
# Set MKLVERSION to version part of MKL directory (if present)
MKLVERSION?=
# Set MKLROOTDIR to possible locations for MKL
MKLROOTDIRS+=/opt/intel/mkl/$(MKLVERSION) /usr/local/intel/Frameworks/mkl
MKLROOTDIRS+=$(subst :, ,$(subst /lib/intel64,,$(LIBRARY_PATH)))
MKLROOTDIRS+=$(subst :, ,$(subst /lib/em64t,,$(LIBRARY_PATH)))

LIBSUBDIRS=intel64 universal em64t x86_64

#
#  Should not be necessary to alter below here
#

ifndef MKLROOT
ifneq ($(MKL_PATH),)
MKLROOT:=$(MKL_PATH)
else
ifneq ($(MKL_INCLUDE),)
MKLROOT:=$(strip $(firstword $(foreach dir_,$(MKL_INCLUDE),$(shell test -e $(dir_)/mkl.h && dirname $(dir_)))))
else
ifneq ($(MKL_LIB),)
MKLROOT:=$(shell test -e $(MKL_LIB)/../.. && dirname $(MKL_LIB)/..)
else
MKLROOT:=$(strip $(firstword $(foreach dir_,$(MKLROOTDIRS),$(shell test -e $(dir_)/include/mkl.h && echo $(dir_)))))
endif
endif
endif
endif

ifeq ($(MKLROOT),)
MKL_INC_priv=$(error "Intel Maths library was requested but not found")
MKL_LIB_priv=$(error "Intel Maths library was requested but not found")
LAPACKLIB=$(error "Intel Maths library was requested but not found")
LAPACKINC=$(error "Intel Maths library was requested but not found")
MATHVER=$(error "Intel Maths library was requested but not found")
else

ifeq ($(MKL_INCLUDE),)
MKL_INC_priv:=-I$(MKLROOT)/include
else
MKL_INC_priv=$(foreach dir_,$(MKL_INCLUDE), -I$(dir_))
endif


ifeq ($(MKL_LIB),)
MKL_LIB_priv=$(strip $(foreach dir_,$(LIBSUBDIRS),$(shell test -d $(MKLROOT)/lib/$(dir_) && echo -L$(MKLROOT)/lib/$(dir_))))
else
MKL_LIB_priv=$(foreach dir_,$(MKL_LIB), -L$(dir_))
endif

LAPACKINC:=-DUSE_MKL $(MKL_INC_priv)
LAPACKLIB:=$(MKL_LIB_priv) -lmkl_intel_lp64 -lmkl_core

MKLROOT_STATIC:=$(strip $(firstword $(foreach dir_,$(LIBSUBDIRS),$(shell test -e $(MKLROOT)/lib/$(dir_)/libmkl_core.a && echo $(dir_)))))
LAPACKLIB_STATIC:= $(MKLROOT)/lib/$(MKLROOT_STATIC)/libmkl_intel_lp64.a $(MKLROOT)/lib/$(MKLROOT_STATIC)/libmkl_core.a  $(MKLROOT)/lib/$(MKLROOT_STATIC)/libmkl_sequential.a

ifdef serial_build
LAPACKLIB+= -lmkl_sequential
else 
ifdef GCCROOT
LAPACKLIB+= -lmkl_gnu_thread -lpthread
else
ifdef INTELROOT
LAPACKLIB+= -lmkl_intel_thread
else
LAPACKLIB+= -lmkl_sequential
endif
endif
endif


MATHVER:="call mkl_get_version_string(libver)"
endif


