COMPILER = mpixlcxx_r
# COMPILER = mpixlC
FLAGS = -O3 -std=c++11

all: task3

task3: task3.cpp GridGenerator.hpp Entities.hpp GridCoordinateDecomposer.hpp BatcherSorter.hpp
	$(COMPILER) $(FLAGS) task3.cpp -o task3

N = 4096
GRID = flower

submit_test: task3
	for k in 4 16 64; do \
		for p in 1 2 4 8 16 32; do \
			mpisubmit.bg -n $$p -w 01:59:00 --stdout /dev/null --stderr /dev/null ./task3 -- $(N) $(N) $$k $(GRID) evaluations_$(N)x$(N)_new_batcher_mergesort.jsonl; \
		done; \
		mpisubmit.bg -n 64 -w 00:30:00 --stdout /dev/null --stderr /dev/null ./task3 -- $(N) $(N) $$k $(GRID) evaluations_$(N)x$(N)_new_batcher_mergesort.jsonl; \
		mpisubmit.bg -n 128 -w 00:15:00 --stdout /dev/null --stderr /dev/null ./task3 -- $(N) $(N) $$k $(GRID) evaluations_$(N)x$(N)_new_batcher_mergesort.jsonl; \
		mpisubmit.bg -n 256 -w 00:10:00 --stdout /dev/null --stderr /dev/null ./task3 -- $(N) $(N) $$k $(GRID) evaluations_$(N)x$(N)_new_batcher_mergesort.jsonl; \
	done


submit: task3
	mpisubmit.bg -n 32 --stdout output.txt --stderr error.txt ./task3 -- 80 8 4 flower flower_80x8_4.txt
	
