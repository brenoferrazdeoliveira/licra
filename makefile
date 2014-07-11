ifeq ($(shell which icc),)
COMPILER = gcc
FLAGS = -O3 -mtune=native
else
COMPILER = icc
FLAGS = -fast
endif

LIB = -lm -lfftw3 -lgsl -lgslcblas

MBBA: dat clean
	@${COMPILER} ${FLAGS} -D LC=0 src/licra.c src/ic.c src/tr.c \
	src/dQ.c src/sd.c src/op.c ${LIB}

5CB: dat clean
	@${COMPILER} ${FLAGS} -D LC=1 src/licra.c src/ic.c src/tr.c \
	src/dQ.c src/sd.c src/op.c ${LIB} 

jm:
	@${COMPILER} ${FLAGS} src/jm.c ${LIB}

clean:
	@rm -f a.out
	@rm -f dat/*.dat

dat:
	@mkdir -p dat
