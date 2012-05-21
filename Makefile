#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include StandardDefinitions.mk

#Site Specific  Flags
ifeq ($(strip $(BOOST_ROOT)),)
	BOOST_ROOT = /usr/local/include
endif
SYSINCLUDES	= -I/usr/include -I$(BOOST_ROOT)
SYSLIBS         = -L/usr/lib
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}

#Generic and Site Specific Flags
CXXFLAGS     += $(SYSINCLUDES) -IAraRootFormat ##$(INC_ARASIM)
LDFLAGS      += -g -I$(BOOST_ROOT) -L${ROOTSYS}/lib $(LD_ARASIM) -L. -lAraSimEvent

# copy from ray_solver_makefile (removed -lAra part)

# added for Fortran to C++


LIBS	= $(ROOTLIBS) -lMinuit $(SYSLIBS)
GLIBS	= $(ROOTGLIBS) $(SYSLIBS)

#ROOT_LIBRARY = libAra.${DLLSUF}

OBJS = Vector.o EarthModel.o IceModel.o Trigger.o Ray.o Tools.o Efficiencies.o Event.o Detector.o Position.o Spectra.o RayTrace.o RayTrace_IceModels.o signal.o secondaries.o Settings.o Primaries.o counting.o RaySolver.o Report.o eventDict.o AraSim.o
CCFILE = Vector.cc EarthModel.cc IceModel.cc Trigger.cc Ray.cc Tools.cc Efficiencies.cc Event.cc Detector.cc Spectra.cc Position.cc RayTrace.cc signal.cc secondaries.cc RayTrace_IceModels.cc Settings.cc Primaries.cc counting.cc RaySolver.cc Report.cc AraSim.cc
CLASS_HEADERS = Trigger.h Detector.h Settings.h Spectra.h IceModel.h Primaries.h Report.h Event.h #need to add headers which added to Tree Branch

PROGRAMS = AraSim

ARAROOTLIB = libAraSimEvent.so

all : $(ARAROOTLIB) $(PROGRAMS) 
	
libAraSimEvent.so : 
	@cd AraRootFormat; make all; make install

AraSim : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAMS)
	@echo "done."

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS)
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(G77LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -dynamiclib -undefined $(UNDEFOPT) $(LDFLAGS) $(G77LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(G77LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif

##-bundle

#%.$(OBJSUF) : %.$(SRCSUF)
#	@echo "<**Compiling**> "$<
#	$(CXX) $(CXXFLAGS) -c $< -o  $@



%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@

%.$(OBJSUF) : %.cc
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@
	@echo "LD_ARASIM=" $(LD_ARASIM)$

# added for fortran code compiling
%.$(OBJSUF) : %.f
	@echo "<**Compiling**> "$<
	$(G77) -c $<


eventDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CLASS_HEADERS) LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
	@cd AraRootFormat; make clean
#############################################################################
