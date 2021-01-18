######################################################################
## Include variables specific to the particular IBAMR build.
include $(IBAMR_BUILD_DIR)/config/make.inc
SRC = $(wildcard *.cpp)
SRC_F_2D = $(wildcard fortran/*2d.f)
SRC_F_3D = $(wildcard fortran/*3d.f)

OBJS = $(patsubst %.cpp,%.o,$(SRC))
OBJS_F_2D = $(patsubst %.f,%.o,$(SRC_F_2D))
OBJS_F_3D = $(patsubst %.f,%.o,$(SRC_F_3D))

CPPFLAGS += -I../include/

LSLib2D.a: $(IBAMR_LIB_2D) $(IBTK_LIB_2D) $(OBJS) $(OBJS_F_2D)
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) -DNDIM=$(PDIM) $(OBJS) $(OBJS_F_2D)
	$(AR) rcs LSLib2D.a $(OBJS) $(OBJS_F_2D)

LSLib3D.a: $(IBAMR_LIB_3D) $(IBTK_LIB_3D) $(OBJS) $(OBJS_F_3D)
	$(CXX) -c $(CXXFLAGS) $(LDFLAGS) -DNDIM=$(PDIM) $(OBJS) $(OBJS_F_3D)
	$(AR) rcs LSLib3D.a $(OBJS) $(OBJS_F_3D)

clean:
	$(RM) stamp-{opt,debug}
	$(RM) *.o *.lo *.objs *.ii *.int.c fortran/*.o *.a
	$(RM) -r .libs
	$(RM) LSLib2D.a LSLib3D.a
