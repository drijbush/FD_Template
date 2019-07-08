PROG =	FD

SRCS =	FD_template.f90

OBJS =	FD_template.o

LIBS =	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = gfortran
#F90FLAGS = -O3 -fopenmp -g -std=f2008 -Wall -Wextra -fcheck=all -finit-real=snan
F90FLAGS = -O3 -fopenmp -g -std=f2008 -Wall -Wextra -finit-real=snan
LDFLAGS = -g -fopenmp


all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

FD_template.o: FD_template.o
