COMPILER = mpixlcxx_r
FLAGS = -O3 -std=c++11

PROCESSORS = 4

all: task2

task2: task2.cpp PointGenertor.hpp Entities.hpp
	$(COMPILER) $(FLAGS) task2.cpp -o $@

submit: task2
	mpisubmit.bg -n $(PROCESSORS) --stdout output$(PROCESSORS).txt --stderr error$(PROCESSORS).txt ./task2 -- 3 5
