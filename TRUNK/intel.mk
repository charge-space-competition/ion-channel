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
# Compiler definitions
#
# set INTELVER to version part of icc name, usually blank
INTELVER?=
# Set INTELROOTDIR to possible locations for INTEL
INTELROOTDIR+=/opt/bin /opt/local/bin /usr/bin /usr/local/bin
# you can also set INTELROOT to empty if in path
# INTELROOT:=

#
#  Should not be necessary to alter below here
#
ifndef INTELROOT
INTELROOT:=$(strip $(foreach dir_,$(INTELROOTDIR),$(shell test -f $(dir_)/icc$(INTELVER) && echo $(dir_))))
endif
INTELROOTDIR:=

# set INTELARCH to the optimisation flag for your architecture, note
#    that Intel's -fast option cannot be used with scientific code.
INTELARCH=
# -m64
#
# Compiler characteristics
# use -DHAVE_LLRINT=1 if you get "error: function llrint(double) already declared"
# use -DHAVE_UNIQPTR=1 if std::unique_ptr is available
CXXFLAGS+=-DHAVE_LLRINT=1

#
# Should not need to change below here
#
CC=$(INTELROOT:=/)icc$(INTELVER)
CXX=$(INTELROOT:=/)icpc$(INTELVER)
CCVERSION=$(CXX)$(shell $(CXX) --version)
CCTARGET=$(shell $(CXX) -dumpmachine)
FC=$(INTELROOT:=/)ifort$(INTELVER)
FCVERSION=$(FC)$(shell $(FC) --version)
FCTARGET=$(shell $(FC) -dumpmachine)

OPTFLAGS:=-gdwarf-2 -O2 -fpp $(INTELARCH) -ansi-alias -unroll -mp1 
DBGFLAGS:=-gdwarf-2 -O1
ifeq ($(VARIANT),DEBUG)
FFLAGS+=-fpp -check bounds -check uninit -ftrapuv
endif
FLAGS:=-fp-model strict -fp-model except 
# -shared-intel
OBJ+=*__genmod.f90 *__genmod.mod

# x86 arch flags
TEST:=$(shell gcc -dM -E -xc /dev/null)
SSEFLAGS:=$(if $(findstring __SSE2__,$(TEST)),-axSSE2)
SSEFLAGS:=$(if $(findstring __SSE3__,$(TEST)),-axSSE3)
SSEFLAGS:=$(if $(findstring __SSSE3__,$(TEST)),-axSSSE3)
SSEFLAGS:=$(if $(findstring __SSE4.1__,$(TEST)),-axSSE4.1)
# -mssse3 -msse4.1 -msse4.2
TEST:=

OPTCFLAGS=$(OPTFLAGS) -Wall
OPTCXXFLAGS=$(OPTFLAGS) -Wall
OPTFFLAGS=$(OPTFLAGS) -warn all
# -pad

ifndef serial_build
OPENMP:= -openmp -parallel -mkl
else
OPENMP:= -mkl=sequential
endif

ifdef static_build
LDFLAGS+=-cxxlib
else
LDFLAGS+=-shared-intel -shared-libgcc -cxxlib
endif


STATIC_BEGIN=
STATIC_END=
