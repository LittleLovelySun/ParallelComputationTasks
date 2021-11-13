#include <iostream>
#include <cstring>
#include "SortingScheduler.hpp"

using namespace std;


void Test() {
	for (int i = 1; i <= 24; i++) {
		SortingSheduler sheduler(i);
		cout << "Test for " << i << ": " << (sheduler.Test() ? "OK" : "FAILED") << endl;
	}
}

void Help() {
	cout << "Usage: ./bsort [n|--test]" << endl;
	cout << "Options: " << endl;
	cout << "    n        size of sorting network, natural" << endl;
	cout << "    --test   test mode for n from 1 to 24" << endl;
	cout << "    --help   print this message" << endl;
}


int main(int argc, char const* argv[]) {
	if (argc != 2) {
		cout << "Invalid arguments. Check --help for usage" << endl;
		return -1;
	}

	if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
		Help();
		return 0;
	}

	if (!strcmp(argv[1], "--test")) {
		Test();
		return 0;
	}

	int n = atoi(argv[1]);
	if (n < 1) {
		cout << "Invalid n value (" << argv[1] << ")" << endl;
		return -1;
	}

	SortingSheduler sheduler(n);
	sheduler.Shedule();
	return 0;
}