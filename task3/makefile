# COMPILER = mpixlcxx_r
COMPILER = mpixlC
FLAGS = -O3 -std=c++11

all: task3

task3: task3.cpp GridGenerator.hpp Entities.hpp GridCoordinateDecomposer.hpp BatcherSorter.hpp
	$(COMPILER) $(FLAGS) task3.cpp -o $@

N = 128
SORT = std

submit_test: task3
	for i in {1..8}; do \
		for n in 1 2 4 8 16 32 64 128; do \
			mpisubmit.bg -n $$n --stdout /dev/null --stderr error$$n.txt ./task2 -- $(N) $(N) 0 t result_$(SORT))_$(N).jsonl $(SORT)); \
		done; \
		mpisubmit.bg -n 256 -w 00:10:00 --stdout /dev/null --stderr error256.txt ./task2 -- $(N) $(N) 0 t result_$(SORT))_$(N).jsonl $(SORT)); \
	done


PROCESSOR = 4

submit: task3
	# mpisubmit.bg -n $(PROCESSOR) --stdout output$(PROCESSOR).txt --stderr error$(PROCESSOR).txt ./task3 -- 8 8 4
	mpisubmit.pl -p $(PROCESSOR) --stdout output$(PROCESSOR).txt --stderr error$(PROCESSOR).txt ./task3 -- 512 512 2
	mpisubmit.pl -p 1 --stdout output1.txt --stderr error1.txt ./task3 -- 512 512 2
	mpisubmit.pl -p 64 --stdout output64.txt --stderr error64.txt ./task3 -- 512 512 2
