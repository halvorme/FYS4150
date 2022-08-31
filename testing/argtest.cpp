#include <iostream>
#include <string>

int main(int argc, char** argv){
	if (argc != 4) {
		std::string exec_name = argv[0];
	
		std::cerr << "Error: Wrong number of args." << std::endl;
		std::cerr << "Usage " << exec_name << " <int> <double> <string>" << std::endl;
		
		return 1;
	}
	
	int heltall = atoi(argv[1]);
	double desimaltall = atof(argv[2]);
	std::string streng = argv[3];
	
	std::cout << "The input was: Heltallet " << heltall << ", desimaltallet " 
		<< desimaltall << " og strengen " << streng << std::endl;
	
	return 0;
}
