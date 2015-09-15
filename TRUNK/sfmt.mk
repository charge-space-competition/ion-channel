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
# SSE Optimised Mersenne Twister 
#

# Set SFMTROOTDIR to possible locations for SFMT
SFMTROOTDIR+=/home/finnerty/code /Users/Scratch/finnerty_scr/code /opt/include /opt/local/include

# Set dSFMTVERSION to version part of dSFMT directory (if present)
dSFMTVER:=-src-2.1

# Set SFMTVERSION to version part of SFMT directory (if present)
SFMTVER:=-src-1.3.3

#
#  Should not be necessary to alter below here
#

ifndef SFMTROOT
# TRY dSFMT FIRST

# Set SFMT to version 
SFMT=dSFMT
# Set USESFMT to 1 for SFMT, 2 for dSFMT
USESFMT=2
# set version.
SFMTVERSION:=$(dSFMTVER)

SFMTROOT:=$(firstword $(strip $(foreach dir_,$(SFMTROOTDIR),$(shell test -d $(dir_)/$(SFMT)$(dSFMTVER) && echo $(dir_)))))

ifeq ($(SFMTROOT),)
# NO dSFMT try SFMT

# Set SFMT to version 
SFMT=SFMT
# Set USESFMT to 1 for SFMT, 2 for dSFMT
USESFMT=1
# set version.
SFMTVERSION:=$(SFMTVER)

SFMTROOT:=$(firstword $(strip $(foreach dir_,$(SFMTROOTDIR),$(shell test -d $(dir_)/$(SFMT)$(SFMTVER) && echo $(dir_)))))

endif
endif
SFMTROOTDIR:=

ifeq ($(SFMTROOT),)
SFMT_INC=$(error "Fatal error: Can not find SFMT or dSFMT directory!")
else
SFMT_INC:=-I$(SFMTROOT)/$(SFMT)$(SFMTVERSION)
endif

RANDFLAGS=-DUSE_SFMT=$(USESFMT) $(SFMT_INC)

