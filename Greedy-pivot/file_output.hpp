#include "data.hpp"
#include <chrono>
#include <unistd.h>


// variable for time measure
std::chrono::system_clock::time_point start, end;

// result
double outlier_count = 0;			// #outliers
double false_positive_count = 0;
double cpu_approximate = 0;
double cpu_exact = 0;
unsigned long int distance_comp_count = 0;
unsigned long int visit_count = 0;
double memory_usage = 0;


// compute memory usage
double process_mem_usage() {

    double resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    resident_set = rss * page_size_kb;

	return resident_set / 1000;
}

// result output
void output_result() {

	std::string f_name = "result/data-id(" + std::to_string(dataset_id) + ")_radius(" + std::to_string(radius) + ")_k(" + std::to_string(k) + ")_thread_no(" + std::to_string(core_no) + ")_greedy.csv";
	
	std::ofstream file;
	file.open(f_name.c_str(), std::ios::out | std::ios::app);

	if (file.fail()) {
		std::cerr << " cannot open the output file." << std::endl;
		file.clear();
		return;
	}

	file << 
	"filtering time [millisec]" << "," << 
	"verification time [millisec]" << "," << 
	"total time [millisec]" << "," << 
	"#outliers" << "," << 
	"outlier rate" << "," <<
	"false positives" << "," <<
	"#distance comp." << "," <<
	"visit count" << "," <<
	"memory [MB]" << "," <<
	"degree" << "\n";
	file <<
	cpu_approximate << "," <<
	cpu_exact << "," <<
	cpu_approximate + cpu_exact << "," <<
	outlier_count << "," <<
	(outlier_count / dataset.size()) * 100 << "," <<
	false_positive_count << "," <<
	distance_comp_count << "," <<
	visit_count << "," <<
	memory_usage << "," <<
	degree << "\n\n";

	file.close();
}
