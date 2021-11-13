COMPILER = mpixlcxx_r
FLAGS = -O3 -Wall

PROCESSORS = 4

all: task2

task2: task2.cpp
	$(COMPILER) $(FLAGS) $^ -o $@

submit: task2
	mpisubmit.bg -n $(PROCESSORS) --stdout output$(PROCESSORS).txt --stderr error$(PROCESSORS).txt ./task2
