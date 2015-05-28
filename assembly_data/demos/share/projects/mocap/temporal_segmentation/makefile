BASE = ../../..

# ODE = 1
# ROS = 1

# DEPEND = Core Mocap Ors Gui Infer Optim Perception
DEPEND = Core Mocap Ors Gui Optim Perception

CXXFLAGS := -I../src
SRCS := $(shell find -iname '*.cpp') $(shell find ../src/ -iname '*.cpp')
OBJS := $(subst .cpp,.o, ${SRCS})

include $(BASE)/gofMake/generic.mk
