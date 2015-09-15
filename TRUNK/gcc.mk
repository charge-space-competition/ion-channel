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
# Configure where your GCC compilers are
#
# set GCCVER to version part of gcc name (if present)
# (Here is minimum required version)
GCCVER?=-4.4
# Set GCCROOTDIR to possible locations for GCC
GCCROOTDIR+=/opt/bin /opt/local/bin /usr/bin /usr/local/bin
# you can also set GCCROOT to empty if in path
# GCCROOT:=

#
#  Should not be necessary to alter below here
#
ifndef GCCROOT
GCCROOT:=$(strip $(firstword $(foreach dir_,$(GCCROOTDIR),$(shell test -f $(dir_)/gcc$(GCCVER) && echo $(dir_)))))
GFCROOT:=$(strip $(firstword $(foreach dir_,$(GCCROOTDIR),$(shell test -f $(dir_)/gfortran$(GCCVER) && echo $(dir_)))))
endif
GCCROOTDIR:=

#
# Compiler characteristics
# use --std=c++0x if you want to use C++ random numbers (gcc version > 4.3)
# use -DHAVE_UNIQPTR=1  if std::unique_ptr is available
ifeq ($(strip $(GCCVER)),-4.2)
CXXFLAGS+= -fstrict-aliasing
FFLAGS+= -Wunderflow -std=f2003 -DHAVE_NO_IMPLICIT_ISNAN
DBGFFLAGS+= -fno-automatic -ffpe-trap=invalid,zero -fbounds-check -Wconversion -mieee-fp -fbounds-check
else
CXXFLAGS+= -std=c++0x -DHAVE_UNIQPTR=1
FFLAGS+= -mssse3 -msse4.1 -msse4.2 -cpp -Warray-temporaries -std=f2003 -fall-intrinsics -Wunderflow -fexternal-blas
DBGFFLAGS+= -fno-automatic -ffpe-trap=invalid,zero -fbounds-check  -Wconversion -mieee-fp -fbounds-check
endif
#
# Should not be necessary to alter below here.
#
CC=$(GCCROOT:=/)gcc$(GCCVER) -m64
CXX=$(GCCROOT:=/)g++$(GCCVER) -m64
DEATH= -ffast-math -floop-interchange -floop-block -floop-strip-mine -ftree-loop-distribution -fgcse-sm -fgcse-las -funsafe-loop-optimizations -Wunsafe-loop-optimizations

#
FC=$(GFCCROOT:=/)gfortran$(GCCVER) -m64
# compiler version information
CCVERSION:=$(CXX) $(shell $(CXX) -dumpversion)
CCTARGET:=$(shell $(CXX) -dumpmachine)
FCVERSION:=$(FC) $(shell $(FC) -dumpversion)
FCTARGET:=$(shell $(FC) -dumpmachine)


OPTFLAGS:=-O2 -fno-common -funroll-loops -finline -finline-limit=600 -fno-math-errno -fprefetch-loop-arrays -finline-functions 
#-flto

DBGFLAGS:=-O0 -DDEBUG=1 -ggdb -gdwarf-2
FLAGS:=-pipe -fmessage-length=0 -Wall -pedantic

# x86 arch SSE flags
TEST:=$(shell gcc -dM -E -xc /dev/null)
SSEFLAGS:=$(if $(findstring __MMX__,$(TEST)),-mmmx)
SSEFLAGS:=$(if $(findstring __SSE__,$(TEST)),-msse)
SSEFLAGS:=$(if $(findstring __SSE2__,$(TEST)),-msse2)
SSEFLAGS:=$(if $(findstring __SSE3__,$(TEST)),-msse3)
ifndef bad_mmx_register
MFPMATH=-mfpmath=sse,387
else
MFPMATH=-mfpmath=sse
endif
# OPTFLAGS+=$(SSEFLAGS) $(if $(findstring __SSE_MATH__,$(TEST)),$(MFPMATH))
OPTFLAGS+= -mmmx -msse -msse2 -msse3 -mtune=generic -march=nocona -ftree-vectorize
# OPTFLAGS+=-flto -mtune=generic -ftree-vectorize
TEST:=
ifndef serial_build
OPENMP:=-fopenmp
else
OPENMP:=
endif

OPTCFLAGS=$(OPTFLAGS) -Wno-format -DDEBUG=0 -Wno-trigraphs 
OPTCXXFLAGS=$(OPTFLAGS) -Wno-format -fvisibility-inlines-hidden -DDEBUG=0 -Wno-trigraphs 
OPTFFLAGS=$(OPTFLAGS) -fmodulo-sched

ifdef static_build
LDFLAGS+=-static-libgcc -static-libgfortran ./libstdc++.a
STATIC_DEFSH?=[ -f libstdc++.a ] || ln -s `$(CXX) -print-file-name=libstdc++.a` .
else
ifdef mac_os
SUFLIB=dylib
else
SUFLIB=so
endif
stdc++_lib_dir=$(dir $(shell $(CC) -print-file-name=libstdc++.$(SUFLIB)))
LDFLAGS+= -L$(stdc++_lib_dir) -lstdc++
endif

ifdef USE_LD_GROUP
STATIC_BEGIN?=-Wl,"-("
STATIC_END?=-Wl,"-)"
else
STATIC_BEGIN:=
STATIC_END:=
endif
