
#define INPUTFILE "input/input1024"

CC ?= gcc
CFLAGS ?= -Wall -std=gnu99 -fopenmp -o

EXEC = basic basic_pthread basic_sse Strassen_openmp Strassen_basic Strassen_pthread StrassenR_basicArray StrassenR_basic Strassen_sseOpenmp StrassenR_pool

all: $(EXEC)

$(EXEC):% : %.c
	$(CC) $(CFLAGS) $@ $^

test1 : basic
	echo 3 | sudo tee /proc/sys/vm/drop_caches
	perf stat \
		-e cache-misses,cache-references,instructions,cycles \
	./basic 1 
	echo 3 | sudo tee /proc/sys/vm/drop_caches
	perf stat \
		-e cache-misses,cache-references,instructions,cycles \
	./basic 2
	rm -f 3

test2 : Strassen_openmp basic_pthread
	./Strassen_openmp 4 1
	./basic_pthread
test3:	Strassen_basic StrassenR_basic
	./Strassen_basic
	./StrassenR_basic
test4:	Strassen_openmp Strassen_pthread
	./Strassen_openmp 4 2
	./Strassen_pthread < ./input/threadNumber
test5:	StrassenR_pool StrassenR_basic
	./StrassenR_basic
	./StrassenR_pool
test6:	basic_sse
	./basic_sse 4
test7:	Strassen_sseOpenmp Strassen_openmp
	./Strassen_openmp 4 2
	./Strassen_sseOpenmp 4

	

clean:	
	rm -f $(EXEC) 3 *.txt
