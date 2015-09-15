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
# GNU BLAS/LAPACK library 
#
# Set GSLVERSION to version part of GSL directory (if present)
GSLVERSION?=
# Set GSLROOTDIR to possible locations for GSL
GSLROOTDIR+=/opt /opt/local /usr /usr/local

#
#  Should not be necessary to alter below here
#
ifndef GSLROOT
GSLROOT=$(strip $(foreach dir_,$(GSLROOTDIR),$(shell test -d $(dir_)/include/gsl$(GSLVERSION) && echo $(dir_))))
endif

ifeq ($(GSLROOT),)
GSL_INC=$(error "GNU Scientific library was requested but not found")
GSL_LIB=$(error "GNU Scientific library was requested but not found") 
LAPACKINC=$(error "GNU Scientific library was requested but not found")
LAPACKLIB=$(error "GNU Scientific library was requested but not found")
MATHVER=$(error "GNU Scientific library was requested but not found")
else
GSL_INC:=$(foreach dir,$(GSLROOT),-I$(dir)/include/gsl$(GSLVERSION))
GSL_LIB:=$(foreach dir,$(GSLROOT),-L$(dir)/lib) 
LAPACKLIB:=$(GSL_LIB) -lgsl -llapack
LAPACKINC:=-DUSE_GSL $(GSL_INC)
MATHVER:="call math_lib_version(libver)"
endif

