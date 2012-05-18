#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include StandardDefinitions.mk

#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS)  $(SYSINCLUDES) -I. #$(FFTFLAG)
LDFLAGS      += -g $(ROOTLDFLAGS) 


LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) #$(FFTLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libAraRootSimEvent.${DLLSUF}
LIB_OBJS =  UsefulAraStationEvent.o UsefulIcrrStationEvent.o araRootDict.o
CLASS_HEADERS = UsefulAraStationEvent.h UsefulIcrrStationEvent.h araIcrrDefines.h


all : $(ROOT_LIBRARY)

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,123456),123456)
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

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


araRootDict.C: $(CLASS_HEADERS)
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



