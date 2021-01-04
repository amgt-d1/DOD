#include <vector>
#include "file_input.hpp"
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cfloat>
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

	// neighbor object list
	std::multimap<float, unsigned int> neighbor_list;

	// reverse kNN list
	std::unordered_set<unsigned int> reverse_knn_list;

	// temporal neighbors
	std::unordered_set<unsigned int> neighbors;

	// threshold
	float threshold = FLT_MAX;

	// update flag
	bool update_flag = 0;

	// pivot flag
	bool pivot_flag = 0;

	// exact flag
	bool exact_flag = 0;

	// edge (vector ver.)
	std::vector<unsigned int> edges;


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


	// get threshold
	void get_threshold() {

		if (neighbor_list.size() >= degree) {
			
			auto itr = neighbor_list.end();
			--itr;
			threshold = itr->first;
		}
	}

	// get accessed object
	void get_accessed(std::unordered_set<unsigned int> &accessed) {

		auto it = neighbor_list.begin();
		while (it != neighbor_list.end()) {
			accessed.insert(it->second);
			++it;
		}

		// insert its own id
		accessed.insert(identifier);
	}

	// neighbor list update
	bool update_neighbor_list(float distance, unsigned int id) {

		bool flag = 0;

		if (distance < threshold){

			// check duplication
			bool f = 1;
			auto p = neighbor_list.equal_range(distance);
			for (auto it = p.first; it != p.second; ++it) {
				if (it->second == id) {
					f = 0;
					break;
				}
			}

			if (f) {

				flag = 1;

				// insert new entry
				neighbor_list.insert({distance, id});

				// remove
				while (neighbor_list.size() > degree) {
					auto it = neighbor_list.end();
					--it;
					neighbor_list.erase(it);
				}

				// threshold update
				get_threshold();

				// update flag
				update_flag = 1;
			}
		}

		return flag;
	}

	// neighbor list update for exact kNN computation
	void update_neighbor_list_(const float distance, const unsigned int id) {

		const unsigned int coefficient = 4;

		if (distance < threshold || neighbor_list.size() < degree * coefficient){

			// check duplication
			bool f = 1;
			auto p = neighbor_list.equal_range(distance);
			for (auto it = p.first; it != p.second; ++it) {
				if (it->second == id) {
					f = 0;
					break;
				}
			}

			if (f) {

				// insert new entry
				neighbor_list.insert({distance, id});

				// remove
				while (neighbor_list.size() > degree * coefficient) {
					auto it = neighbor_list.end();
					--it;
					neighbor_list.erase(it);
				}

				// threshold update
				if (neighbor_list.size() == degree * coefficient) {
					auto it = neighbor_list.end();
					--it;
					threshold = it->first;
				}
			}
		}
	}

	// copy neighbors
	void copy_neighbors() {

		// clear neighbors
        neighbors.clear();

		// copy
		auto it = neighbor_list.begin();
		while (it != neighbor_list.end()) {
			neighbors.insert(it->second);
			++it;
		}
	}

};


// definition of dataset
std::vector<data> dataset;


// input data
void input_data() {

	std::string prefix = "../dataset/";

	// id variable
	unsigned int id = 0;

	// input
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
}

// distance computation
float compute_distance(const data& d1, const data& d2) {

	float distance = 0;

	if (type == 0) {

		// L2 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += powf(d1.pt[i] - d2.pt[i], 2.0);
		distance = sqrt(distance);
	}
	else if (type == 1) {

		// L1 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += fabsf(d1.pt[i] - d2.pt[i]);
	}
	else if (type == 2) {

		// jaccard distance
		double overlap = 0;

		if (d1.set.size() <= d2.set.size()) {

			auto itr = d1.set.begin();
			while (itr != d1.set.end()) {

				if (d2.set.find(*itr) != d2.set.end()) ++overlap;

				++itr;
			}
		}
		else {

			auto itr = d2.set.begin();
			while (itr != d2.set.end()) {

				if (d1.set.find(*itr) != d1.set.end()) ++overlap;

				++itr;
			}
		}

		float jaccard = (float)(overlap / (d1.set.size() + d2.set.size() - overlap));
		distance = 1.0 - jaccard;
	}
	else if (type == 3) {

		// edit distance
		const size_t len1 = d1.str.size();
		const size_t len2 = d2.str.size();
		std::vector<std::vector<size_t>> d(len1 + 1, std::vector<size_t>(len2 + 1));

		d[0][0] = 0;
		for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
		for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

		for (unsigned int i = 1; i <= len1; ++i) {
			for (unsigned int j = 1; j <= len2; ++j) d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (d1.str[i - 1] == d2.str[j - 1] ? 0 : 1) });
		}

		distance = (float)d[len1][len2];
	}
	else if (type == 4) {

		// angular distance

		// inner product
		float ip = 0;

		// L2 norm
		float l2norm_1 = 0;
		float l2norm_2 = 0;

		for (unsigned int i = 0; i < dimensionality; ++i) {

			ip += d1.pt[i] * d2.pt[i];

			l2norm_1 += powf(d1.pt[i], 2.0);
			l2norm_2 += powf(d2.pt[i], 2.0);
		}

		l2norm_1 = sqrt(l2norm_1);
		l2norm_2 = sqrt(l2norm_2);

		float cosine_sim = ip / (l2norm_1 * l2norm_2);

		distance = acos(cosine_sim) / 3.1416;
	}
	else if (type == 5) {

		// L4 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += powf(d1.pt[i] - d2.pt[i], 4.0);
		distance = powf(distance, 0.25);
	}

	return distance;
}

// distance computation with early termination
float compute_distance(const data& d1, const data& d2, const float &threshold) {

	float distance = 0;

	if (type == 0) {

		// L2 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += powf(d1.pt[i] - d2.pt[i], 2.0);
		distance = sqrt(distance);
	}
	else if (type == 1) {

		// L1 norm
		for (unsigned int i = 0; i < dimensionality; ++i) distance += fabsf(d1.pt[i] - d2.pt[i]);
	}
	else if (type == 2) {

		// jaccard distance
		double overlap = 0;

		if (d1.set.size() <= d2.set.size()) {

			auto itr = d1.set.begin();
			while (itr != d1.set.end()) {

				if (d2.set.find(*itr) != d2.set.end()) ++overlap;

				++itr;
			}
		}
		else {

			auto itr = d2.set.begin();
			while (itr != d2.set.end()) {

				if (d1.set.find(*itr) != d1.set.end()) ++overlap;

				++itr;
			}
		}

		float jaccard = (float)(overlap / (d1.set.size() + d2.set.size() - overlap));
		distance = 1.0 - jaccard;
	}
	else if (type == 3) {

		// edit distance
		const size_t len1 = d1.str.size();
		const size_t len2 = d2.str.size();
		std::vector<std::vector<size_t>> d(len1 + 1, std::vector<size_t>(len2 + 1));

		d[0][0] = 0;
		for (unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
		for (unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

		for (unsigned int i = 1; i <= len1; ++i) {
			for (unsigned int j = 1; j <= len2; ++j) d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (d1.str[i - 1] == d2.str[j - 1] ? 0 : 1) });
		}

		distance = (float)d[len1][len2];
	}
	else if (type == 4) {

		// angular distance

		// inner product
		float ip = 0;

		// L2 norm
		float l2norm_1 = 0;
		float l2norm_2 = 0;

		for (unsigned int i = 0; i < dimensionality; ++i) {

			ip += d1.pt[i] * d2.pt[i];

			l2norm_1 += powf(d1.pt[i], 2.0);
			l2norm_2 += powf(d2.pt[i], 2.0);
		}

		l2norm_1 = sqrt(l2norm_1);
		l2norm_2 = sqrt(l2norm_2);

		float cosine_sim = ip / (l2norm_1 * l2norm_2);

		distance = acos(cosine_sim) / 3.1416;
	}
	else if (type == 5) {

		// L4 norm
		for (unsigned int i = 0; i < dimensionality; ++i) {
			distance += powf(d1.pt[i] - d2.pt[i], 4.0);

			if (i % 30 == 0) {
				if (powf(threshold, 4.0) < distance) break;
			}
		}
		distance = powf(distance, 0.25);
	}

	return distance;
}

// avoid integer
float add_noise(const data& d1, const data& d2) {

	float noise = 0;

	if (type == 3 || type == 4) {

        std::default_random_engine engine(d1.identifier + d2.identifier);
        std::uniform_real_distribution<> rn(0, 0.9999);
        noise = rn(engine);
    }

	return noise;
}