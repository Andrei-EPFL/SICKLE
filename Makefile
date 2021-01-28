CC = gcc
CXX = g++
CFLAGS = -O3 -Wall -g -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
CXXFLAGS = -O3 -Wall -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
#OPTS = -DSINGLE_PREC -lfftw3f -lfftw3f_omp -fopenmp -DOMP
OPTS = -DDOUBLE_PREC -lfftw3 -lfftw3_omp -fopenmp -DOMP
OPTS += -I/Users/czhao/lib/fftw-3.3.8/include -L/Users/czhao/lib/fftw-3.3.8/lib -I/opt/local/include -L/opt/local/lib

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
