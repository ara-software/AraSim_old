#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
ifeq ($(strip $(BOOST_ROOT)),)
	BOOST_ROOT = /usr/local/include
endif
SYSINCLUDES	= -I/usr/include -I$(BOOST_ROOT)
SYSLIBS         = -L/usr/lib
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}

#Now the bits we're actually compiling


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES)
LDFLAGS      += -g $(ROOTLDFLAGS) 

LIBS          = $(ROOTLIBS) -lMinuit $(SYSLIBS) 
#LIBS          = $(ROOTLIBS) -lgsl -lMinuit $(SYSLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

ROOT_LIBRARY = libAra.${DLLSUF}
#LIB_OBJS = AraSim.o Detector.o Event.o Efficiencies.o Trigger.o IceModel.o EarthModel.o eventDict.o
LIB_OBJS =  Vector.o EarthModel.o IceModel.o Trigger.o Ray.o Tools.o Efficiencies.o Event.o Detector.o Position.o Spectra.o RayTrace.o RayTrace_IceModels.o signal.o eventDict.o Settings.o
CCFILE       =  Vector.cc EarthModel.cc IceModel.cc Trigger.cc Ray.cc Tools.cc Efficiencies.cc Event.cc Detector.cc Spectra.cc Position.cc RayTrace.cc signal.cc RayTrace_IceModels.cc Settings.cc
CLASS_HEADERS = Trigger.h Detector.h

PROGRAMS = araSim

all : $(ROOT_LIBRARY) araSim

araSim : $(ROOT_LIBRARY) AraSim.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD) $(CXXFLAGS) $(LDFLAGS) AraSim.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif

#%.$(OBJSUF) : %.$(SRCSUF)
#	@echo "<**Compiling**> "$<
#	$(CXX) $(CXXFLAGS) -c $< -o  $@



%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@

%.$(OBJSUF) : %.cc
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@



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
#############################################################################



