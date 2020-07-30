CC = gcc
CFLAGS = -O3 -Wall -std=c99 -g
#OPTS = -lfftw3f
OPTS = -DSINGLE_PREC -lfftw3f -lfftw3f_omp -fopenmp -DOMP
#OPTS = -DDOUBLE_PREC -lfftw3 -lfftw3_omp -fopenmp -DOMP

LIBS = -lm -lgsl -lgslcblas
INCL =
OPTS += $(LIBS) $(INCL)
SRCS = coca.c linhalo.c load_conf.c rand_field.c read_data.c sel_halo.c
EXEC = linhalo

all:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(OPTS)

clean:
	rm $(EXEC)
