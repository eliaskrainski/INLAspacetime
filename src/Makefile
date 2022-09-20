toInclude = ${R_LIBRARY_DIR}/INLA/include/

sources = $(wildcard *.c)

all : $(sources:.c=.o) so

%.o : %.c
	@echo gnerating $@
	gcc -c -O3 -Wall -Wextra -fpic -lcblas -I$(toInclude) $^ -o $@

so : $(sources:.c=.o)
	@echo generate .so
	gcc -shared -lcblas *.o -o cgenericModels.so


clean :
	@##-rm -f *.o

.PHONY: all clean
