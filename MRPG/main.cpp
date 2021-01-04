#include "mrpg.hpp"


int main() {

	// show the current time
	get_current_time();

	// file input
	input_parameter();

	// data input
	input_data();

	std::cout << " --------------------\n";
	std::cout << " dataset id: " << dataset_id << "\n";
	std::cout << " dataset cardinality: " << dataset.size() << "\n";
	std::cout << " dataset dimensionality: " << dimensionality << "\n";
	std::cout << " degree: " << degree << "\n";
	std::cout << " #threads: " << core_no << "\n";
	std::cout << " --------------------\n\n";

	// MRPG construction
	build_mrpg();

	// result output
	output_result();

	// graph output
	output_graph();

	return 0;
}