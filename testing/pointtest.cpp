#include <iostream>

int main() {
	int x = 42;
	int* pointer = &x;
	
	int a = *pointer;

	std::cout << a << std::endl;

	return 0;
}
