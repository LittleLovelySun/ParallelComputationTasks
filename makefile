COMPILER = g++
FLAGS = -O3 -Wall

all: task1 task2

task1: task1.cpp
	$(COMPILER) $(FLAGS) $^ -o $@

task2: task2.cpp
	$(COMPILER) $(FLAGS) $^ -o $@
