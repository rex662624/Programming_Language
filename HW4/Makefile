
#define INPUTFILE "input/input1024"

CC ?= gcc
CFLAGS ?= -Wall -std=gnu99 -fopenmp -o

EXEC = basic basic_pthread Strassen_openmp Strassen_basic Strassen_pthread StrassenR_basic2 StrassenR_basic Strassen_sseOpenmp StrassenR_pool

all: $(EXEC)

$(EXEC):% : %.c
	$(CC) $(CFLAGS) $@ $^

test1 : basic
	echo 3 | sudo tee /proc/sys/vm/drop_caches
	perf stat \
		-e cache-misses,cache-references,instructions,cycles \
	./basic 2 
	echo 3 | sudo tee /proc/sys/vm/drop_caches
	perf stat \
		-e cache-misses,cache-references,instructions,cycles \
	./basic 1


clean:	
	rm -f $(EXEC)
