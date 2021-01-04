#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <random>
#include <time.h>


// dataset identifier
unsigned int dataset_id = 0;

// data dimensionality
unsigned int dimensionality = 2;

// system parameter
unsigned int core_no = 1;

// range radius parameter
float radius = 0;

// minPts parameter
unsigned int k = 50;

// distance function type
unsigned int type = 0;

// graph degree
unsigned int degree = 20;


// get current time
void get_current_time() {

	time_t t = time(NULL);
	printf(" %s\n", ctime(&t));
}

/// parameter input
void input_parameter() {

	std::string prefix = "../parameter/";

	std::ifstream ifs_radius(prefix + "radius.txt");
	std::ifstream ifs_k(prefix + "k.txt");
	std::ifstream ifs_dataset_id(prefix + "dataset_id.txt");
	std::ifstream ifs_core_no(prefix + "thread_num.txt");
	std::ifstream ifs_degree(prefix + "degree.txt");

	if (ifs_radius.fail()) {
		std::cout << " radius.txt does not exist." << std::endl;
		std::exit(0);
	}
	else if (ifs_k.fail()) {
		std::cout << " k.txt does not exist." << std::endl;
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
	else if (ifs_degree.fail()) {
		std::cout << " degree.txt does not exist." << std::endl;
		std::exit(0);
	}

	while (!ifs_radius.eof()) { ifs_radius >> radius; }
	while (!ifs_k.eof()) { ifs_k >> k; }
	while (!ifs_dataset_id.eof()) { ifs_dataset_id >> dataset_id; }
	while (!ifs_core_no.eof()) { ifs_core_no >> core_no; }
	while (!ifs_degree.eof()) { ifs_degree >> degree; }

	// set dimensionality
	if (dataset_id == 0) {
		dimensionality = 128;
	}
}