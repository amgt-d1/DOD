#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <time.h>


// dataset identifier
unsigned int dataset_id = 0;

// data dimensionality
unsigned int dimensionality = 2;

// kNN graph average degree
unsigned int degree = 1;

// system parameter
unsigned int core_no = 1;

// distance function type (0: L2, 1: L1, 2: jaccard distance, 3: edit distance, 4: angular distance, 5: L4)
unsigned int type = 0;


// get current time
void get_current_time() {

	time_t t = time(NULL);
	printf(" %s\n", ctime(&t));
}

// parameter input
void input_parameter() {

	std::string prefix = "../parameter/";

	std::ifstream ifs_degree(prefix + "degree.txt");
	std::ifstream ifs_dataset_id(prefix + "dataset_id.txt");
	std::ifstream ifs_core_no(prefix + "thread_num.txt");

	if (ifs_degree.fail()) {
		std::cout << " degree.txt does not exist." << std::endl;
		std::exit(0);
	}
	else if (ifs_dataset_id.fail()) {
		std::cout << " dataset_id.txt does not exist." << std::endl;
		std::exit(0);
	}
	else if (ifs_core_no.fail()) {
		std::cout << " core_no.txt does not exist." << std::endl;
		std::exit(0);
	}

	while (!ifs_degree.eof()) { ifs_degree >> degree; }
	while (!ifs_dataset_id.eof()) { ifs_dataset_id >> dataset_id; }
	while (!ifs_core_no.eof()) { ifs_core_no >> core_no; }

	// set dimensionality
	if (dataset_id == 0) {
		dimensionality = 128;
	}
}