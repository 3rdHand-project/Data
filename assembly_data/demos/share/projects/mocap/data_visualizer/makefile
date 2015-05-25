BASE = ../../..

ODE = 1
# ROS = 1

DEPEND = Core Mocap Ors Gui Perception

SRCS := $(shell find -iname '*.cpp')
OBJS := $(subst .cpp,.o, ${SRCS})

include $(BASE)/gofMake/generic.mk
