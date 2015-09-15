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
# DEFAULT CONFIGURATION SETTINGS.
#
# These values are used as default values if no other is given  
#
MAJORVERSION=17
MINORVERSION:=867

TOUCH=touch
VARIANT=RELEASE

# include machine specific settings
-include local.mk

ifdef static_build
MINORVERSION:=$(MINORVERSION).static
endif

# Add any extra files you want deleted with "make clean" here
OBJ:=
#
# Compiler data
#
ifeq ($(VARIANT),DEBUG)
CFLAGS+=$(FLAGS) $(DBGFLAGS)
CXXFLAGS+=$(FLAGS) $(DBGFLAGS)
FFLAGS+=$(FLAGS) $(DBGFLAGS)
endif

ifeq ($(VARIANT),RELEASE)
CFLAGS+=$(FLAGS) $(OPTCFLAGS) $(OPENMP)
CXXFLAGS+=$(FLAGS) $(OPTCXXFLAGS) $(OPENMP)
FFLAGS+=$(FLAGS) $(OPTFFLAGS) $(OPENMP)
endif

#
# Math lib data (defines LAPACKLIB and LAPACKINC)
#
CFLAGS+=$(LAPACKINC)
CXXFLAGS+=$(LAPACKINC)
ifdef static_build
LDFLAGS+=$(LAPACKLIB_STATIC)
else
LDFLAGS+=$(LAPACKLIB)
endif

#
# Random library (defines RANDINC), only necessary if
# you do not have C++ random numbers
CFLAGS+=$(RANDFLAGS)
CXXFLAGS+=$(RANDFLAGS)

#
#  Reconfiguring for your system should not require to alter below here
#

all:
	@echo "Ion Channel Version: $(MAJORVERSION).$(MINORVERSION)"
	@echo "Usage:"
	@echo " make channel (also generates channel$(MAJORVERSION).$(MINORVERSION))"
	@echo " make clean (removes temporary files)"
	@echo " make distclean (removes all generated files)"

#
# Create channelXXX
#
CLEANABLE += vers.f90
FORCE :
vers.o : FORCE
	echo "recreating vers.f90"
	@echo "module vers" > vers.f90
	@echo "character(*), parameter :: compnm=\"channel.$(MAJORVERSION).$(MINORVERSION)\"" >> vers.f90
	@echo "character(*), parameter :: compdt=\"$(shell "date")\"" >> vers.f90
	@echo "character(*), parameter :: compvr=\"exported\"" >> vers.f90
	@echo "character(*), parameter :: cmpclr=\"&" >> vers.f90
	@printf "& %s %s %s %s %s %s %s %s&\\n&" $(foreach opt,$(CCVERSION),"$(opt)") >> vers.f90
	@echo "\"" >> vers.f90
	@echo "character(*), parameter :: cmpctg=\"$(CCTARGET)\"" >> vers.f90
	@echo "character(*), parameter :: cmpclg=\"&" >> vers.f90
	@printf "& %s %s %s %s&\\n"  $(foreach opt,$(subst ",',$(CXXFLAGS)),"$(opt)") >> vers.f90
	@echo "&\"" >> vers.f90
	@echo "character(*), parameter :: cmpflr=\"&" >> vers.f90
	@printf "& %s %s %s %s %s %s %s %s&\\n" $(foreach opt,$(FCVERSION),"$(opt)") >> vers.f90
	@echo "&\"" >> vers.f90
	@echo "character(*), parameter :: cmpftg=\"$(FCTARGET)\"" >> vers.f90
	@echo "character(*), parameter :: cmpflg=\"&" >> vers.f90
	@printf "& %s %s %s %s&\\n"  $(foreach opt,$(subst ",',$(FFLAGS)),"$(opt)") >> vers.f90
	@echo "&\"" >> vers.f90
	@echo "contains" >> vers.f90
	@echo "function libver()" >> vers.f90
	@echo "implicit none" >> vers.f90
	@echo "character(256) :: libver" >> vers.f90
	@echo $(MATHVER)  >> vers.f90
	@echo "endfunction" >> vers.f90
	@echo "end module" >> vers.f90
	$(COMPILE.f) -o $@ vers.f90

CHNX=channel$(MAJORVERSION).$(MINORVERSION)
CHN90SRC=const.f90 grid.f90 strngs.f90 rnstat.f90 geom.f90 spec.f90 conf.f90 patch.f90 trial.f90 accum.f90 vers.f90 channel.f90 simstate.f90
const.o : const.f90
grid.o : grid.f90
simstate.o : simstate.f90
strngs.o : strngs.f90 | const.o
rnstat.o : rnstat.f90 | const.o
SPECSRC?=spec.f90
spec.o :: $(SPECSRC) | const.o strngs.o
geom.o :: geom.f90 | spec.o const.o strngs.o 
conf.o :: conf.f90 | spec.o geom.o const.o strngs.o
PATCHSRC?=patch.f90
patch.o :: $(PATCHSRC) | conf.o spec.o geom.o const.o strngs.o simstate.o $(LAPACKOBJ)
TRIALSRC?=trial.f90
trial.o :: $(TRIALSRC) | patch.o conf.o spec.o geom.o const.o strngs.o rnstat.o simstate.o
accum.o :: accum.f90 | trial.o patch.o conf.o spec.o geom.o const.o strngs.o simstate.o
channel.o :: channel.f90 | accum.o trial.o patch.o conf.o spec.o geom.o const.o strngs.o simstate.o
CHNXSRC=cutil.cpp 

$(CHNX) : DEPENDS+=$(CHNSRC) $(CHNXSRC) $(CHNCSRC) vers.f90
ifdef static_build
# To get static libstdc++ use of the following creates a symlink in the current directory
# to the compiler specific file.  Use -static-libstdc++ for gcc > 4.5(?)
# ln -s `g++ -print-file-name=libstdc++.a`
LDFLAGS:= $(STATIC_BEGIN) $(LDFLAGS)

$(CHNX) : $(CHN90SRC:.f90=.o) $(CHNSRC:.f=.o) $(CHNXSRC:.cpp=.o) $(CHNCSRC:.c=.o) $(LAPACKOBJ)
	$(STATIC_DEFSH)
	$(LINK.f) -o $@ $^ $(LDHOST) $(STATIC_LAPACK) $(STATIC_END)

else
LDFLAGS:= $(STATIC_BEGIN) $(LDFLAGS)
$(CHNX) : $(CHN90SRC:.f90=.o) $(CHNSRC:.f=.o) $(CHNXSRC:.cpp=.o) $(CHNCSRC:.c=.o) $(LAPACKOBJ)
	$(LINK.f) -o $@ $^ $(LDHOST) $(LAPACK) $(STATIC_END)
endif

channel : $(CHNX)
	ln -sf $(CHNX) channel

PROG+=$(CHNX)
OBJ+=$(CHN90SRC:.f90=.o) $(CHNSRC:.f=.o) $(CHNXSRC:.cpp=.o) $(CHNCSRC:.c=.o) 
CLEANABLE+=$(CHN90SRC:.f90=.mod)

#
#  Clean up directory
#
.PHONY : clean distclean
clean:
	$(RM) $(RMFLAG) $(OBJ) $(CLEANABLE)

distclean: clean
	$(RM) $(RMFLAG) $(PROG) $(OPTPROG) make.*.log
%.o %.mod : %.f90
	$(COMPILE.f) -o $@ $^ 
%.o %.mod : %.F90
	$(COMPILE.f) -o $@ $^ 

