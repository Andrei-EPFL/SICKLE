CC = gcc
CXX = g++
CFLAGS = -O3 -Wall
CXXFLAGS = -O3 -Wall
#OPTS = -DSINGLE_PREC -lfftw3f -lfftw3f_omp -fopenmp -DOMP
OPTS = -DDOUBLE_PREC -lfftw3 -lfftw3_omp -fopenmp -DOMP

LIBS = -lm -lgsl -lgslcblas
INCL = -Imake_survey
OPTS += $(LIBS) $(INCL)
SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)
EXEC = linhalo

all: cuboid.o $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $^ $(OPTS) $(INCL)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $*.c $(OPTS) $(INCL)

cuboid.o:
	$(CXX) $(CXXFLAGS) -c make_survey/cuboid.cpp

clean:
	rm $(EXEC) $(OBJS) cuboid.o
