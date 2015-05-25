BASE   = ../..
NAME   = lib$(shell basename `pwd`)
OUTPUT := $(BASE)/lib/$(NAME).so

GL = 1
LAPACK = 1
OPENCV = 0
#PCL = 1
DEPEND = Core Gui extern_libcolorseg Ors

CXXFLAGS += -D__STDC_CONSTANT_MACROS -DHAVE_LIBAV -DHAVE_X264
LIBS += -lavformat -lavcodec -lavutil -lswscale -lx264 -lz

SRCS = $(shell find . -maxdepth 1 -name '*.cpp' )
OBJS = $(SRCS:%.cpp=%.o)

include $(BASE)/gofMake/generic.mk
