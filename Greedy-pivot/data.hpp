#include <vector>
#include "file_input.hpp"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>


// definition of data object
class data {

public:

	// for numeric data
	std::vector<float> pt;

	// for string data
	std::vector<char> str;

	// for set data
	std::unordered_set<unsigned int> set;

	// identifier
	unsigned int identifier = 0;

	// #distance computation
	unsigned long int distance_comp_count = 0;

	// #visits
	unsigned long int visit_count = 0;
	unsigned long int visit_count_first = 0;

	// edge set
	std::vector<unsigned int> edges;

	// threshold
	float threshold = 0;

	// flag (0: if approximate counting < k)
	bool flag = 1;

	// pivot flag
	bool pivot_flag = 0;

	// exact flag
	bool exact_flag = 0;
	

	/***************/
	/* constructor */
	/***************/

	// standard
	data() {}

	// with identifier
	data(unsigned int id) { identifier = id; }

	/*******************/
	/* member function */
	/*******************/

	// pt update
	void update_pt(std::vector<float>& point) {

		for (unsigned int i = 0; i < point.size(); ++i) pt.push_back(point[i]);
	}

	// word update
	void update_str(std::string s) {

		for (unsigned int i = 0; i < s.size(); ++i) str.push_back(s[i]);
	}

	// set update
	void update_set(std::unordered_set<unsigned int>& s) { set = s; }

	bool operator>(const data &d) const { return threshold > d.threshold; }
};


// definition of dataset
std::vector<data> dataset;


// input data
void input_data() {

	std::string prefix = "../dataset/";

	// id variable
	unsigned int id = 0;

	std::cout << " data input begins...\n";

	// data input
	if (dataset_id == 0) {

		// input deep10M
		std::string f_name = prefix + "sift_base.fvecs";

		// file input
		std::ifstream ifs_file(f_name, std::ios::in | std::ios::binary);

		// error check
		if (ifs_file.fail()) {
			std::cout << " SIFT does not exist." << std::endl;
			std::exit(0);
		}

		float t;
		unsigned int cnt = 0;
		std::vector<float> temp;

		while (!ifs_file.eof()) {

			// read
			ifs_file.read((char*)& t, sizeof(float));

			if (cnt > 0) {
				
				temp.push_back(t);
			
				if (cnt % 129 == 0) {

					temp.pop_back();

					// make a data object
					data d(id);
					d.update_pt(temp);

					// insert into dataset
					dataset.push_back(d);

					// increment identifier
					++id;

					temp.clear();
					cnt = 0;
				}
			}
			++cnt;
		}
	}

	std::cout << " data input ends...\n\n";
}

// input graph
void input_graph() {

	std::string prefix = "../MRPG/result/graph/";

	// input graph
	std::string f_name = prefix + "data-id(" + std::to_string(dataset_id) + ")_degree(" + std::to_string(degree) + ")_graph.txt";

	// file input
	std::ifstream ifs_file(f_name);

	// error check
	if (ifs_file.fail()) {
		std::cout << " graph file does not exist." << std::endl;
		std::exit(0);
	}

	std::cout << " graph input begins...\n";

	unsigned int id = 0;
	unsigned int id_ = 0;
	bool flag = 0;
	bool flag_ = 0;

	// file read
	while (!ifs_file.eof()) {

		// input
		ifs_file >> id >> id_ >> flag >> flag_;

		// edge insertion
		if (id != id_) dataset[id].edges.push_back(id_);

		// pivot flag
		dataset[id].pivot_flag = flag;

		// exact flag
		dataset[id].exact_flag = flag_;
	}

	std::cout << " graph input ends...\n\n";
}
