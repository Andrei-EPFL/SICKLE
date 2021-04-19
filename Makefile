CC = gcc
CFLAGS = -O3 -Wall -std=c99 -g
OPTS = -DSINGLE_PREC -lfftw3f -lfftw3f_omp -fopenmp -DOMP
#OPTS = -DDOUBLE_PREC -lfftw3 -lfftw3_omp -fopenmp -DOMP

LIBS = -L/global/common/sw/cray/cnl7/haswell/gsl/2.5/gcc/8.2.0/sr445ay/lib -lm -lgsl -lgslcblas
INCL = -I/global/common/sw/cray/cnl7/haswell/gsl/2.5/gcc/8.2.0/sr445ay/include
FFTW_DIR = /opt/cray/pe/fftw/3.3.8.4/$(CRAY_CPU_TARGET)
ifneq ($(FFTW_DIR),)
  LIBS += -L$(FFTW_DIR)/lib
  INCL += -I$(FFTW_DIR)/include
endif
OPTS += $(LIBS) $(INCL)
SRCS = coca.c linhalo.c load_conf.c rand_field.c read_data.c sel_halo.c
EXEC = linhalo

all:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(OPTS)

clean:
	rm $(EXEC)
