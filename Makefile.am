

AM_LDFLAGS = -lm -lgmp -lmpfr -lmpfi -lwcpg $(LAPACK_LIBS) $(BLAS_LIBS) 

bin_PROGRAMS = fxpf_test
fxpf_test_SOURCES = fxpf_test.c \
					fxpf.c \
					fxpf.h \
					filter.h\
					filter.c\
					LPconstructor.c\
					LPconstructor.h\
					matrixalg.c\
					matrixalg.h\
					mpfi_matrixalg.c\
					mpfi_matrixalg.h

lib_LTLIBRARIES = libfxpf.la
libfxpf_la_SOURCES = fxpf.c \
					fxpf.h \
					filter.h\
					filter.c\
					LPconstructor.c\
					LPconstructor.h\
					matrixalg.c\
					matrixalg.h\
					mpfi_matrixalg.c\
					mpfi_matrixalg.h

include_HEADERS = fxpf.h filter.h matrixalg.h mpfi_matrixalg.h

fxpf_test_CFLAGS = $(AM_CFLAGS) 
fxpf_test_CPPFLAGS = $(AM_CPPFLAGS) 
fxpf_test_LDFLAGS = $(AM_LDFLAGS) 

libfxpf_la_LDFLAGS = -version-info 0:9:0 $(AM_LDFLAGS)


ACLOCAL_AMFLAGS = -I m4 

LIBTOOL_DEPS = @LIBTOOL_DEPS@
libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status libtool