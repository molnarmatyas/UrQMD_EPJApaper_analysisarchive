WrkDir := $(shell pwd)
SrcDir = sources
ObjDir = object
DepDir = object/deps
ExeDir = exe
ClsDir = classes

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs) #-lMinuit2
ROOTGLIBS := $(shell root-config --glibs) #-lMinuit2

CXX           = g++
CFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++17
RFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++17 $(ROOTCFLAGS)
LD            = g++
SRCSUF        = cc
LDFLAGS       = -O3 #-m32
RLIBS         = $(ROOTLIBS) $(SYSLIBS)
RGLIBS        = $(ROOTLIBS) $(SYSLIBS)

CLASSES  = Levy_proj_reader
FIT_SOURCES = $(shell ls $(SrcDir)/*.cc)
CLASS_SOURCES = $(addprefix $(ClsDir)/,$(addsuffix .$(SRCSUF), $(CLASSES)))
PROGRAMS = $(notdir $(basename $(FIT_SOURCES)))

define COMPILE_TEMPLATE
-include $(DepDir)/$(notdir $(basename $(1))).d
$(ObjDir)/$(notdir $(basename $(1))).o:
	$(CXX) -c $(RFLAGS) -I$(WrkDir)/$(SrcDir) -I$(WrkDir)/$(ClsDir) -MD -MP -MF $(DepDir)/$(notdir $(basename $(1))).d $(1) -o $$@
endef

define MAKE_TEMPLATE
-include $(DepDir)/$(1).d
$(ObjDir)/$(1).o:
	@echo ""
	$(CXX) $(RFLAGS) -c -I$(WrkDir)/$(SrcDir) -I$(WrkDir)/$(ClsDir) -MD -MP -MF $(DepDir)/$(1).d $(SrcDir)/$(1).cc -o $$@

$(ExeDir)/$(1).exe: $(ObjDir)/$(1).o $(ObjDir)/Levy_proj_reader.o
	@echo ""
	$(LD) $(LDFLAGS) $$^ $(RGLIBS) -o $$@
endef
$(foreach class, $(CLASS_SOURCES), $(eval $(call COMPILE_TEMPLATE,$(class))))
$(foreach program, $(PROGRAMS), $(eval $(call MAKE_TEMPLATE,$(program))))

all: $(addsuffix .exe, $(PROGRAMS))

clean:
	@rm -f $(ExeDir)/*.exe $(ObjDir)/*.o $(DepDir)/*.d

