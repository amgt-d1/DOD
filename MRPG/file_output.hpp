#include "data.hpp"
#include <chrono>


// variable for time measure
std::chrono::system_clock::time_point start, end;

// result
double cpu_initial_edge = 0;
double cpu_update_edge = 0;
double cpu_exact_knn = 0;
double cpu_reverse = 0;
double cpu_connect = 0;
double cpu_detour = 0;
double cpu_remove = 0;


// result output
void output_result() {

	std::string f_name = "result/graph-construct_time_data-id(" + std::to_string(dataset_id) + ")_degree(" + std::to_string(degree) + ")_thread_no(" + std::to_string(core_no) + ").csv";
	std::ofstream file;
	file.open(f_name.c_str(), std::ios::out | std::ios::app);

	if (file.fail()) {
		std::cerr << " cannot open the output file." << std::endl;
		file.clear();
		return;
	}

	file <<
		"init time [millisec]" << "," <<
		"update time [millisec]" << "," <<
		"exact compu. time [millisec]" << "," <<
		"reverse edge time [millisec]" << "," <<
		"connect time [millisec]" << "," <<
		"detour remove time [millisec]" << "," <<
		"edge remove time [millisec]" << "," <<
		"total time [millisec]" << "\n";
	file << 
		cpu_initial_edge << "," << 
		cpu_update_edge << "," << 
		cpu_exact_knn << "," << 
		cpu_reverse << "," <<
		cpu_connect << "," << 
		cpu_detour << "," <<
		cpu_remove << "," <<
		cpu_initial_edge + cpu_update_edge + cpu_exact_knn + cpu_reverse + cpu_connect + cpu_detour + cpu_remove << "\n\n";

	file.close();
}

// graph output
void output_graph() {

	std::string f_name = "result/graph/data-id(" + std::to_string(dataset_id) + ")_degree(" + std::to_string(degree) + ")_graph.txt";
	std::ofstream file;
	file.open(f_name.c_str(), std::ios::out | std::ios::app);

	if (file.fail()) {
		std::cerr << " cannot open the graph output file." << std::endl;
		file.clear();
		return;
	}

	for (unsigned int i = 0; i < dataset.size(); ++i) {
		for (unsigned int j = 0; j < dataset[i].edges.size(); ++j) file << i << "\t" << dataset[i].edges[j] << "\t" << dataset[i].pivot_flag << "\t" << dataset[i].exact_flag << "\n";
	}
}
