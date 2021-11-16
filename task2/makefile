COMPILER = mpixlcxx_r
FLAGS = -O3 -std=c++11

all: task2

task2: task2.cpp PointGenertor.hpp Entities.hpp
	$(COMPILER) $(FLAGS) task2.cpp -o $@

N = 128
SORT = std

submit_test: task2
	for i in {1..8}; do \
		for n in 1 2 4 8 16 32 64 128; do \
			mpisubmit.bg -n $$n --stdout /dev/null --stderr error$$n.txt ./task2 -- $(N) $(N) 0 t result_$(SORT))_$(N).jsonl $(SORT)); \
		done; \
		mpisubmit.bg -n 256 -w 00:10:00 --stdout /dev/null --stderr error256.txt ./task2 -- $(N) $(N) 0 t result_$(SORT))_$(N).jsonl $(SORT)); \
	done


PROCESSOR = 8

submit: task2
	mpisubmit.bg -n $(PROCESSOR) --stdout output$(PROCESSOR).txt --stderr error$(PROCESSOR).txt ./task2 -- $N $N 0 t
