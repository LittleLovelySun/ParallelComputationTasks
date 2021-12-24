COMPILER = mpixlcxx_r
FLAGS = -O3 -std=c++11

all: task3

task3: task3.cpp GridGenerator.hpp Entities.hpp GridCoordinateDecomposer.hpp BatcherSorter.hpp
	$(COMPILER) $(FLAGS) task3.cpp -o task3

GRID = flower

submit_test: task3
	for n in 64 128 256; do \
		for k in 4 17 64; do \
			for method in i c; do \
				for p in 1 2 4 8 16 32; do \
					mpisubmit.bg -n $$p -w 01:59:00 --stdout /dev/null --stderr /dev/null ./task3 -- $$n $$n $$k $(GRID) $$method; \
				done; \
				mpisubmit.bg -n 64 -w 00:30:00 --stdout /dev/null --stderr /dev/null ./task3 -- $$n $$n $$k $(GRID) $$method; \
				mpisubmit.bg -n 128 -w 00:15:00 --stdout /dev/null --stderr /dev/null ./task3 -- $$n $$n $$k $(GRID) $$method; \
				mpisubmit.bg -n 256 -w 00:10:00 --stdout /dev/null --stderr /dev/null ./task3 -- $$n $$n $$k $(GRID) $$method; \
			done; \
		done; \
		sleep 1800; \
	done


submit: task3
	mpisubmit.bg -n 32 --stdout output.txt --stderr error.txt ./task3 -- 80 8 4 flower i decompose.txt
