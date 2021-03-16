#include "vptree.hpp"


// Greedy couting
unsigned int greedy_counting(const unsigned int id) {

	// max hop-count
	unsigned int counter = 0;
	float size = 0;
	while (1) {
		size = powf(degree, counter);
		if (size >= k) break;
		++counter;
	}
	const unsigned int max_hop = counter;

	// probability
	std::mt19937 mt(id);
    std::uniform_real_distribution<> prob(0, 1.0);

	// init count
	unsigned int count = 0;

	// prepare dist
	float dist = 0;

	// prepare hash table
	std::unordered_set<unsigned int> visited;
	visited.reserve(k);
	visited.insert(id);

	// init queue
	std::deque<std::pair<unsigned int, unsigned int>> queue;
	queue.push_back({id, 0});

	while(queue.size() > 0) {

		// pop the front
		const unsigned int identifier = queue[0].first;
		const unsigned int hop_count = queue[0].second;
		queue.pop_front();

		for (unsigned int i = 0; i < dataset[identifier].edges.size(); ++i) {

			// get id
			const unsigned int id_ = dataset[identifier].edges[i];

			// increment visit
			++dataset[id].visit_count;

			// check visited
			if (visited.find(id_) == visited.end()) {

				// insert into hash table
				visited.insert(id_);

				// distance computation
				dist = compute_distance(dataset[id], dataset[id_]);
				++dataset[id].distance_comp_count;

				// update queue
				if (dist <= radius) {
					++count;
					if (count == k) break;
					queue.push_back({id_, hop_count + 1});
				}
				else {
					if (dataset[id_].pivot_flag && hop_count < max_hop) {
						if (prob(mt) <= (float)(k / size)) queue.push_back({id_, hop_count + 1});
					}
				}
			}
		}

		if (count == k) break;
	}

	return count;	
}

// 1-hop couting
unsigned int onehop_counting(const unsigned int id) {

	// init count
	unsigned int count = 0;

	// prepare dist
	float dist = 0;

	for (unsigned int i = 0; i < dataset[id].edges.size(); ++i) {

		// get id
		const unsigned int id_ = dataset[id].edges[i];

		// compute distance
		dist = compute_distance(dataset[id], dataset[id_], radius);
		++dataset[id].distance_comp_count;

		if (dist <= radius) ++count;
		if (count == k) break;
	}

	return count;	
}

// Approximate counting
void detect_outlier_approximate() {

	start = std::chrono::system_clock::now();

	#pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(dynamic) reduction(+:outlier_count_approx) reduction(+:outlier_count)
		for (unsigned int i = 0; i < dataset.size(); ++i) {

			if (dataset[i].exact_flag == 0) {

				// approximate counting
				const unsigned int count = greedy_counting(i);

				if (count < k) dataset[i].flag = 0;
			}
			else {

				const unsigned int count = onehop_counting(i);

				if (count < k) outlier_count = outlier_count + 1;
			}
		}
	}

	end = std::chrono::system_clock::now();
	cpu_approximate = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " filtering time: " << cpu_approximate << "[msec]\n\n";
}

// Exact counting
void detect_outlier_exact(vp_tree &vpt, const unsigned int &mode) {

    start = std::chrono::system_clock::now();

    #pragma omp parallel num_threads(core_no)
	{
		#pragma omp for schedule(dynamic) reduction(+:outlier_count) reduction(+:false_positive_count)
		for (unsigned int i = 0; i < dataset.size(); ++i) {

			if (dataset[i].flag == 0) {

				unsigned int count = 0;

				// exact counting
				if (mode == 1) {

					// range search on VP-tree
					count = vpt.range_search(dataset[i]);
				}
				else {

					float distance = 0;

					// scan
					for (unsigned int j = 0; j < dataset.size(); ++j) {

						if (i != j) {
							distance = compute_distance(dataset[i], dataset[j], radius);
							++dataset[i].distance_comp_count;
							if (distance <= radius) ++count;
							if (count == k) break;
						}
					}
				}

				if (count < k) {
					outlier_count = outlier_count + 1;
				}
				else{
					false_positive_count = false_positive_count + 1;
				}
			}
		}
	}

    end = std::chrono::system_clock::now();
	cpu_exact = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << " verification time: " << cpu_exact << "[msec]\n\n";
}


int main () {

    // show the current time
	get_current_time();

	// file input
	input_parameter();

	// data input
	input_data();

	// graph input
	memory_usage = process_mem_usage();
	input_graph();

	std::cout << " --------------------------\n";
	std::cout << " dataset id: " << dataset_id << "\n";
	std::cout << " dataset cardinality: " << dataset.size() << "\n";
	std::cout << " dataset dimensionality: " << dimensionality << "\n";
	std::cout << " radius: " << radius << "\n";
	std::cout << " k: " << k << "\n";
	std::cout << " degree: " << degree << "\n";
	std::cout << " #threads: " << core_no << "\n";
	std::cout << " --------------------------\n\n";

    // build vp-tree if necessary
	unsigned int mode = 0;
	vp_tree vpt;
	
	//if (dataset_id == 0) {
	//	vpt.build();
	//	mode = 1;
	//}

	memory_usage = process_mem_usage() - memory_usage;

	// approximate counting
	detect_outlier_approximate();

    // exact counting
    detect_outlier_exact(vpt, mode);

    // output result
    output_result();

    return 0;
}
