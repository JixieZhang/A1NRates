########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################
# need to define MODELDIR, separated by ':' or space

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

########################################################################
EXECFILE    := a1nrates
LIBFILE     := libA1NRates.so
LIBNAME     := A1NRates
USERDICT    := $(LIBNAME)_Dict
VERSION     := 1.0.0

########################################################################
SRCDIR      := src
INCDIR      := include
OBJDIR      := obj.$(ARCH)

########################################################################
MODELDIR    := XSModel
#if this program need other libs, you need to define them outside, separated 
#by ':' or space. For example, 
#MODELDIR    := XSModel:HRSTransport
ifdef MODELDIR
  MODELLIST   := $(subst :, ,$(MODELDIR))
  MODELINCDIR := $(foreach n,$(MODELLIST),$(n) $(n)/include)
  INCDIR      += . $(MODELINCDIR)
  MODELLIBS   := $(foreach n,$(MODELLIST),-L$(n)/$(OBJDIR) -l$(n))
endif

########################################################################
# Compiler
AR          := ar
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCDIRS     := $(patsubst %,-I%,$(subst :, ,$(INCDIR)))
CFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS)
CXXFLAGS    := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS) 
FFLAGS      := -Wall -fPIC -O3 -g $(MODE) $(INCDIRS)
ifeq ($(MYOS),Darwin) 
#in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -fno-range-check -funroll-all-loops -fdollar-ok \
               -ffixed-line-length-none -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
GPPFLAGS    := -MM
LDFLAGS     := -O3 -g $(MODE)

########################################################################
# Generate obj file list
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)
# header files
HEADERS     := $(foreach n,$(subst :, ,$(INCDIR)),$(wildcard $(n)/*.hh))
# add .o to all the source files
OBJS        := $(addsuffix .o, $(basename $(SOURCES)))
OBJS        := $(patsubst  $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran
OTHERLIBS   := $(MODELLIBS)

########################################################################
# ROOT configure
ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLIBS    := $(shell root-config --libs)
ROOTGLIBS   := $(shell root-config --glibs) -lMinuit

CXXFLAGS    += $(ROOTCFLAGS) 
LIBS        := $(SYSLIBS) $(ROOTLIBS)
GLIBS       := $(SYSLIBS) $(ROOTGLIBS)

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean test vc
VPATH       := $(SRCDIR)

########################################################################
all: lib exe

########################################################################
# Make the $(TARGET).d file and include it.
$(OBJDIR)/%.d: %.c 
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.C 
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cc
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.CC
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cpp 
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

$(OBJDIR)/%.d: %.cxx
	@echo Making dependency for file $< ......
	@set -e; \
	$(CXX) $(GPPFLAGS) $(CXXFLAGS) $< | \
	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.f
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

#$(OBJDIR)/%.d: %.F
#	@echo Making dependency for file $< ......
#	@set -e; \
#	$(FF) -cpp $(GPPFLAGS) $(FFLAGS) $< | \
#	sed 's!$*\.o!$(OBJDIR)/& $@!' > $@; \
#	[ -s $@ ] || rm -f $@

ifneq ($(DEPS),)
  -include $(DEPS)
endif

########################################################################
exe: dir $(OBJS) $(OBJDIR)/Main.o $(OBJDIR)/$(USERDICT).o
	@$(LD) $(LDFLAGS) -o $(EXECFILE) $(OBJDIR)/Main.o $(OBJS)\
           $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@echo "Linking $(EXECFILE) ...... done!"

$(OBJDIR)/Main.o: Main.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
lib: dir $(OBJS) $(OBJDIR)/$(USERDICT).o
	@for model in $(MODELLIST); do \
		if [[ -d $$model ]]; then \
			make lib -s -C $$model; \
		fi; \
	done;
	@$(LD) -shared $(LDFLAGS) -o $(LIBFILE).$(VERSION) \
           $(OBJS) $(OBJDIR)/$(USERDICT).o $(LIBS) $(OTHERLIBS)
	@ln -sf $(LIBFILE).$(VERSION) $(LIBFILE)
	@echo "Linking $(LIBFILE) ...... done!"


$(USERDICT).cxx: $(HEADERS) $(LIBNAME)_LinkDef.h
		@echo "Generating dictionary $(USERDICT).cxx ......"
		@$(ROOTSYS)/bin/rootcint -f $@ -c $(CXXFLAGS) $^

$(OBJDIR)/$(USERDICT).o: $(USERDICT).cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

########################################################################
$(OBJDIR)/%.o: %.c
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CPP
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.CXX
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@  $(CXXFLAGS)

$(OBJDIR)/%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

$(OBJDIR)/%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@  $(FFLAGS)

########################################################################
dir:
	@if [ ! -d $(OBJDIR) ] ; then mkdir -p $(OBJDIR) ;fi

########################################################################
clean: dir
	@rm -f $(OBJDIR)/*
	@rm -f $(USERDICT).cxx $(USERDICT).h
	@rm -f $(EXECFILE) $(LIBFILE) $(LIBFILE).*
	@rm -f *~ *# */*~ */*#

distclean: clean
	@for model in $(MODELDIR); do \
		if [[ -d $$model ]]; then \
			make clean -s -C $$model; \
		fi; \
	done;


test:	
	@echo \\MYOS\:$(MYOS) \\ARCH\:$(ARCH)
	@echo ==================================
	@echo \\CFLAGS\:$(CFLAGS)
	@echo ==================================	
	@echo \\CXXFLAGS\:$(CXXFLAGS)   
	@echo ==================================     
	@echo \\FFLAGS\:$(FFLAGS)
	@echo ==================================
	@echo \\LDFLAGS\:$(LDFLAGS)
	@echo ==================================
	@echo \\SYSLIBS\:$(SYSLIBS)
	@echo ==================================
	@echo \\fsources\: $(FSOURCES)	
	@echo ==================================
	@echo \\sources\: $(SOURCES)
	@echo ==================================
	@echo \\headers\: $(HEADERS)
	@echo ==================================
	@echo \\objs\: $(OBJS)	
	@echo ==================================
	@echo \\dependencies: \$(DEPS)
	@echo ==================================

help: test
	@echo $^ =?= $<

env: test
