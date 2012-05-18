#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include StandardDefinitions.mk


LD_ARA_UTIL=-L$(ARA_UTIL_LIB_DIR) 
INC_ARA_UTIL=-I$(ARA_UTIL_INC_DIR) 


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(FFTFLAG) $(SYSINCLUDES) -I. 
LDFLAGS      += -g $(ROOTLDFLAGS) 


LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) #$(FFTLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libAraRootSimEvent.${DLLSUF}
LIB_OBJS =  UsefulAraStationEvent.o UsefulIcrrStationEvent.o araRootDict.o
CLASS_HEADERS = UsefulAraStationEvent.h UsefulIcrrStationEvent.h araIcrrDefines.h


all : $(ROOT_LIBRARY)

tests: printStructSizes


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

install: $(ROOT_LIBRARY)
	install -d $(ARA_UTIL_LIB_DIR)
	install -d $(ARA_UTIL_INC_DIR)
ifeq ($(PLATFORM),macosx)
	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) $(ARA_UTIL_LIB_DIR)
else
	install -c -m 755 $(ROOT_LIBRARY) $(ARA_UTIL_LIB_DIR)
endif
	install -c -m 644  $(CLASS_HEADERS) $(ARA_UTIL_INC_DIR)
	install -d $(ARA_UTIL_CALIB_DIR)/ICRR/TestBed
	for file in calib/ICRR/TestBed/*.dat calib/ICRR/TestBed/*.csv calib/ICRR/TestBed/*.txt; do install -c -m 644 "$${file}" $(ARA_UTIL_CALIB_DIR)/ICRR/TestBed; done
	install -d $(ARA_UTIL_CALIB_DIR)/ICRR/Station1
	for file in calib/ICRR/Station1/*.dat calib/ICRR/Station1/*.csv calib/ICRR/Station1/*.txt; do install -c -m 644 "$${file}" $(ARA_UTIL_CALIB_DIR)/ICRR/Station1; done
	install -d $(ARA_UTIL_CALIB_DIR)/ATRI
	for file in calib/ATRI/*.csv calib/ATRI/*.txt; do install -c -m 644 "$${file}" $(ARA_UTIL_CALIB_DIR)/ATRI; done
	install -d $(ARA_UTIL_CALIB_DIR)
	for file in calib/*.csv; do install -c -m 644 "$${file}" $(ARA_UTIL_CALIB_DIR); done




clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
#############################################################################



