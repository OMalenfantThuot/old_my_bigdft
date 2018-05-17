check_PROGRAMS = test_overlap test_laplacian test_coulomb test_electrostatics $(extra_tests)
extra_tests = test_ir_index test_R_1_norm

sources = test_Coulomb.f90 test_Electrostatics.f90 test_Overlap.f90 test_Laplacian.f90 \
	test_ir_index.f90 test_R_1_norm.f90

LIBS += -lGaIn
AM_LDFLAGS= -L$(top_builddir)/src/
AM_FCFLAGS= -I$(top_builddir)/src/

$(sources:.f90=.o): $(top_builddir)/src/libGaIn.a

test_coulomb_SOURCES= test_Coulomb.f90
test_electrostatics_SOURCES=test_Electrostatics.f90
test_overlap_SOURCES= test_Overlap.f90
test_laplacian_SOURCES= test_Laplacian.f90

test_ir_index_SOURCES= test_ir_index.f90
test_R_1_norm_SOURCES= test_R_1_norm.f90

check:
	./test_ir_index 
	./test_R_1_norm 
	./test_overlap
	./test_laplacian
	./test_electrostatics
	./test_coulomb
