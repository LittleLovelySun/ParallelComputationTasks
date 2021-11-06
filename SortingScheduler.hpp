#include <iostream>
#include <vector>

using namespace std;

struct Comparator {
	int first;
	int second;

	Comparator(int first, int second) {
		this->first = first;
		this->second = second;
	}
};

class SortingSheduler {
	int n;
	vector<Comparator> comparators;

	void Merge(int start1, int start2, int n, int m, int step = 1);
	void Sort(int start, int size);

	void PrintComparators() const;
	int GetTacts() const;

	bool IsSorted(const vector<int>& array) const;
public:
	SortingSheduler(int n);
	void Shedule();
	bool Test();
};

void SortingSheduler::Merge(int leftPos, int rightPos, int leftSize, int rightSize, int offset) {
	if (leftSize == 0 || rightSize == 0)
		return;

	if (leftSize == 1 && rightSize == 1) {
		comparators.push_back(Comparator(leftPos, rightPos));
		return;
	}

	Merge(leftPos, rightPos, (leftSize + 1) / 2, (rightSize + 1) / 2, offset * 2);
	Merge(leftPos + offset, rightPos + offset, leftSize / 2, rightSize / 2, offset * 2);

	for (int i = 1; i < leftSize - 1; i += 2)
		comparators.push_back(Comparator(leftPos + i * offset, leftPos + (i + 1) * offset));

	if (leftSize % 2 == 0)
		comparators.push_back(Comparator(leftPos + (leftSize - 1) * offset, rightPos));

	for (int i = 1 - leftSize % 2; i < rightSize - 1; i += 2)
		comparators.push_back(Comparator(rightPos + i * offset, rightPos + (i + 1) * offset));
}


void SortingSheduler::Sort(int position, int size) {
	if (size < 2)
		return;

	int middle = size / 2;
	Sort(position, middle);
	Sort(position + middle, size - middle);
	Merge(position, position + middle, middle, size - middle);
}

void SortingSheduler::PrintComparators() const {
	for (auto comparator = comparators.begin(); comparator != comparators.end(); comparator++)
		cout << comparator->first << " " << comparator->second << endl;
}

int SortingSheduler::GetTacts() const {
	vector<int> used(n, 0);

	for (auto comparator = comparators.begin(); comparator != comparators.end(); comparator++) {
		int next = max(used[comparator->first], used[comparator->second]) + 1;
		used[comparator->first] = next;
		used[comparator->second] = next;
	}

	int maxTact = 0;

	for (int i = 0; i < n; i++)
		maxTact = max(maxTact, used[i]);

	return maxTact;
}

SortingSheduler::SortingSheduler(int n) {
	this->n = n;
}

void SortingSheduler::Shedule() {
	Sort(0, n);

	cout << n << " 0 0" << endl;
	PrintComparators();

	cout << comparators.size() << endl;
	cout << GetTacts() << endl;
}


bool SortingSheduler::IsSorted(const vector<int>& array) const {
	for (int i = 1; i < array.size(); i++)
		if (array[i - 1] > array[i])
			return false;

	return true;
}

bool SortingSheduler::Test() {
	Sort(0, n);

	int total = 1 << n;

	for (int i = 0; i < total; i++) {
		vector<int> array(n);
		for (int j = 0; j < n; j++)
			array[j] = (i >> j) & 1;

		for (auto comparator = comparators.begin(); comparator != comparators.end(); comparator++)
			if (array[comparator->first] > array[comparator->second])
				swap(array[comparator->first], array[comparator->second]);

		if (!IsSorted(array))
			return false;
	}

	return true;
}