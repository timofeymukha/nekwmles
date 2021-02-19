FFLAGS += -I..
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)
FFLAGS += -I$(NEK_SOURCE_ROOT)/core
FFLAGS += -I$(KTH_TOOLBOX_DIR)/driver/frame

LATEST_PFUNIT_DIR := $(PFUNIT_DIR)/PFUNIT-4.2
include $(LATEST_PFUNIT_DIR)/include/PFUNIT.mk

SRCS := spalding.f newton.f
OBJS := $(SRCS:%.f=%.o)

%.o : ../%.f
	$(FC) -c $(FFLAGS) $<

%.o : %.F90
	$(FC) -c $(FFLAGS) $<

libsut.a: $(OBJS)
	$(AR) -r $@ $?


test_all_TESTS := test_spalding.pf test_newton.pf
test_all_OTHER_LIBRARIES := -L. -lsut
$(eval $(call make_pfunit_test,test_all))

all: libsut.a test_all

clean:
	$(RM) *.o *.mod *.a *.inc *.F90 test_all
